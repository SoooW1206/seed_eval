# seed-eval: details

Full documentation, moved from the top-level `README.md`. Paths are relative to the repository root.

A parameterized evaluator for **seed selection schemes in `vg giraffe`** on long reads.
For each *condition* (one seed-index configuration) it builds the index with
`vg minimizer`, maps simulated reads with known truth, compares the alignments to
the truth, and collects accuracy, speed, memory, index size and seed density into
one table and plot.

Conditions are rows in `conditions.tsv`. Adding a scheme or parameter setting
(another `k`/`w`/`s`, weighted minimizers, a different vg build with a new seed
type) means adding a row; the driver does not need to change as long as the new
scheme is built by `vg minimizer` and mapped by `vg giraffe`.

Current setup (Paten lab rotation, Oct 2026):

- Graph: HPRC v2.0 Minigraph-Cactus CHM13 graph with HG002 held out, haplotype-sampled for HG002
  (`hprc-v2.0-mc-chm13-eval-sampled16o_fragmentlinked-for-real-hifi-HG002-full.gbz`);
  CHM13 is its `REFERENCE` sample.
- Reads: simulated HG002 HiFi reads (pbsim, real Revio error model, >= 10 kb) with
  truth positions on CHM13; 1k, 10k, 100k and 1M read sets.
- Mapping: `vg giraffe -b hifi`, unchanged across conditions, so only the seed
  index differs between conditions.

## Requirements

- A built `vg` (tested with `v1.77.0-26-gc91cc1d91`). Keep a fixed copy for
  experiments (e.g. `~/bin/vg_baseline`); do not replace a binary that running jobs use.
- Slurm (`sbatch`, `sacct`) for everything heavy. Most jobs use the `medium` partition
  with 16 threads. `filters` and the two small `prepare` jobs use `short`, `fingerprint`
  uses 2 CPUs, and `timing --real` uses `long` with 64 threads.
- Python 3.9+. matplotlib for plots (`collect` and the prototype `--plot` / `--deciles`
  options); numpy for `minimaxmer.py --counts-mmc` and `plot_reads_count.py`. `jq` on the
  compute nodes (`density`, `filters` and `fingerprint` check for it).
- A graph `.gbz`; its distance index is built next to it as `<graph>.dist` by `prepare`.
- Truth reads as GAM files with positions on a `REFERENCE` path of the graph.

## Quick start

1. Copy the example config and fill in your paths (tab-separated `key`, `value`):

   ```bash
   cp config.example.tsv config.tsv
   ```

   | Key | Meaning |
   |---|---|
   | `data_dir` | directory for large outputs (indexes, results, jobs, logs) |
   | `graph` | the graph `.gbz`; the distance index is `<graph>.dist` next to it |
   | `reads_dir` | directory with the truth GAMs; the default read set is `reads_dir/HG002-sim-hifi-1k.gam` |
   | `vg` | vg executable used by rows whose `vg_binary` is `default`, and by `prepare` |
   | `exclude_nodes` | optional: Slurm node list for `sbatch --exclude`, to keep runs on one node type. Used by `--fixed-nodes`, and required by `timing` unless `--no-exclusive` |

   `config.tsv` is ignored by git. Another file can be passed with `--config FILE`.

2. Edit `conditions.tsv` (one row per seed configuration, see below).

3. Run the stages in order. Each prints the sbatch scripts it writes; add `--submit`
   to submit, and wait for the jobs to finish before the next step.

   ```bash
   python3 seed_eval.py prepare --submit                    # once per graph
   python3 seed_eval.py index --submit                      # once per condition
   python3 seed_eval.py map --reads READS.gam --submit      # for each read set
   python3 seed_eval.py density --reads READS.gam --submit  # for each read set
   python3 seed_eval.py filters --reads READS.gam --submit  # for each read set
   python3 seed_eval.py collect --reads READS.gam --snapshot
   ```

Other settings stay at the top of `seed_eval.py`: `PARTITION`, `THREADS`,
`GIRAFFE_PRESET`, `COMPARE_RANGE`, `MAP_MEM` and the time limits.

## Layout

```
seed-eval/                      this repo (code + small summaries)
  seed_eval.py                  the driver
  conditions.tsv                one row per condition
  config.example.tsv            template for config.tsv (your paths; not in git)
  prototype/                    minimaxmer prototype, count-table tool, analysis scripts, tests
  patches/vg/, patches/gbwtgraph/   minimaxmer patches for vg and GBWTGraph, with their licenses
  results_summary/<reads>/      committed snapshots of the summary tables and plots
  LICENSE                       MIT

data_dir/                       large files, not in git
  graph.paths_meta.tsv, graph.length.txt   outputs of prepare
  indexes/<name>.min|.zipcodes|.params|.vg_sha256|.counts_sha256
  counts/                       KMC k-mer counts (.kff) and minimaxmer count tables (.mmc)
  annotation/                   CHM13 annotation BEDs (T2T censat, segmental duplications)
  results/<reads>/<name>/       per-condition map outputs, logs, density/ subdirectory
  results/<reads>/              results.tsv, comparison.png, filter/minimizer summaries
  results/<reads>/timing/       exclusive-node timing jobs (collect reads these)
  results/<reads>/timing_repeats_runs/   repeated shared-node timing jobs (timing_repeats.py)
  jobs/                         generated sbatch scripts
  logs/                         Slurm logs
```

## One-time preparation per graph (`prepare`)

`python3 seed_eval.py prepare [--submit]` writes three jobs to `data_dir/jobs/`. They
do not depend on the reads; run them once for a new graph.

| Job | What it does | Used by |
|---|---|---|
| `prepare_check_paths` | `vg paths -M` into `graph.paths_meta.tsv`, and prints the `REFERENCE` paths so you can check that the truth's reference (CHM13 here) is one | sanity check for `gamcompare` |
| `prepare_build_dist` | `vg index -j` distance index `<graph>.dist` (~97 GiB peak, 44 min on this graph) | every index/map job |
| `prepare_graph_length` | `vg stats -l` total sequence length into `graph.length.txt` | denominator of `density_graph` |

## conditions.tsv

Tab-separated, one row per condition. All columns must be present; optional ones may be empty.

| Column | Required | Meaning |
|---|---|---|
| `name` | yes | unique name; used for index and result file names |
| `vg_binary` | yes | absolute path of the vg used for every step of this condition, or `default` for the `vg` in the config |
| `index_args` | yes | arguments passed verbatim to `vg minimizer` (e.g. `-k 31 -w 20`, `-c -k 31 -s 11`, `-W -k 31 -w 20`) |
| `index_mem` | yes | `--mem` of the index job |
| `giraffe_args` | no | extra arguments passed verbatim to `vg giraffe` (empty in the current sweep) |
| `map_mem` | no | `--mem` of the map and density jobs; empty means `MAP_MEM` (64G) |
| `density_group` | no | groups conditions for tables and plots, e.g. `2/21`; `reference` marks the baseline |
| `label` | no | the seed method, e.g. `minimizer`, `closed syncmer`, `weighted minimizer` |

**Adding a condition:** add a row with a new `name`. If you change `index_args` of an
existing row, also change its `name` (or delete its old index); otherwise the
`.params` check below stops the run.

Memory measured on this graph at k = 31 (density 2/21): index build ~56 GiB plain,
~150 GiB with `-W` (dominated by k-mer counting, so similar for any `w`); giraffe
~41 GB. The table uses 80G / 256G for index jobs (128G for plain 2/11 rows) and
96G for mapping at density 2/11.

Theoretical density: minimizers 2/(w+1), closed syncmers 2/(k-s+1), so `w = k - s`
pairs them (GBWTGraph's `include/gbwtgraph/minimizer.h`, line 1358: "Minimizers and
closed syncmers should have roughly the same seed density when w = k - s.").
`vg minimizer --help` gives `k` a maximum of 31. GBWTGraph's
`MinimizerHeader::sanitize` (`src/minimizer.cpp`) adjusts other values rather than
rejecting them: `k` below 2 becomes 2, `w` or `s` of 0 becomes 1, and `s >= k`
becomes `k - 1`. Nothing requires `w >= k`; the vg default is k29 w11.

**`minimax_k31w20` is the exception.** It uses `vg_minimax` (vg + gbwtgraph on the
`minimaxmer` branches) with `--minimaxmer-counts` in both `index_args` and
`giraffe_args`, and its `w` is 20, not the formula value for its group. Its
`density_group` is `2/11` because its measured density matches that group: on the 1k
simulated reads minimaxmers at w = 20 select 0.186 k-mers per bp (1.99 times
`mini_k31w20`), against 0.181-0.182 for the 2/11 rows on 10k reads. Compare it with
the 2/11 rows, with the 2/21 rows as the same-`w` reference. Its index build peaked at
120 GiB and its giraffe needs the 2/11 memory plus ~34 GB for the counts table, hence
160G for both.

## Stages

Run from the repo directory. Every stage except `collect`, `checksums` and `determinism` only **writes** sbatch
scripts to `data_dir/jobs/` and prints them; add `--submit` to submit.
`--only NAME ...` restricts to some conditions, `--conditions FILE` uses another table,
`--reads GAM` picks the truth read set (default: 1k). On `--submit`, `--fixed-nodes`
adds `--exclude=<exclude_nodes>` and `--sbatch-args "..."` adds any other arguments to
the `sbatch` command line; the generated scripts do not change.

| Stage | Command | Output |
|---|---|---|
| index | `python3 seed_eval.py index [--force] --submit` | `indexes/<name>.min`, `.zipcodes`, `.params`, `.vg_sha256` (and `.counts_sha256` with minimaxmer counts); skipped if already built with the same parameters |
| checksums | `python3 seed_eval.py checksums [--only NAME ...] [--force]` | runs locally: writes missing `.vg_sha256` / `.counts_sha256` records for indexes built before they existed |
| map | `python3 seed_eval.py map --reads GAM [--reps N] --submit` | `results/<reads>/<name>/`: `mapped.gam`, `annot.gam`, `compare.tsv`, `gamcompare.txt`, `giraffe.log`, `stats.txt`, `*.time`, `report.tsv`, `run_info.tsv` |
| density | `python3 seed_eval.py density --reads GAM --submit` | `results/<reads>/<name>/density/`: `prov.gam`, `per_read.tsv`, `density_summary.tsv` |
| filters | `python3 seed_eval.py filters --reads GAM --submit` | one job for all conditions: `results/<reads>/filter_summary.tsv`, `minimizer_summary.tsv` (+ per-condition copies in `density/`) |
| collect | `python3 seed_eval.py collect --reads GAM [--snapshot]` | runs locally: `results/<reads>/results.tsv`, `mcnemar.tsv`, `comparison.png`, `density_sweep.png`, and `timing.tsv`, `timing_runs.tsv`, `placement.tsv` when timing jobs exist; `--snapshot` copies these and the filter/minimizer summaries to `results_summary/<reads>/` |
| fingerprint | `python3 seed_eval.py fingerprint --reads GAM [--run-tag TAG] --submit` | one job per run that has a `mapped.gam` (`<name>/rep*`, or `<name>/TAG/rep*`): `fingerprint.tsv` in the run directory |
| determinism | `python3 seed_eval.py determinism --reads GAM --run-tag TAG` | runs locally: `results/<reads>/determinism/summary.tsv`, `reads.tsv` |
| timing | `python3 seed_eval.py timing --reads GAM --only A B C REF [--real [--sampled]] [--no-exclusive] [--jobs N] [--rounds R] --submit` | one exclusive job, submitted twice (`--real`: two jobs in opposite orders; `--jobs N`: N scripts, submitted by hand): `results/<reads>/timing/job<ID>/runs.tsv` (`--no-exclusive`: `timing_repeats_runs/job<ID>/`), `run_info.tsv`, per run `round<r>_<name>/giraffe.log`, `giraffe.time` (`--real`: also `fingerprint.tsv`, statistics, `placement.tsv`) |

What each job runs:

- **index**: `vg minimizer <index_args> -p -t 16 -d <dist> -o <name>.min -z <name>.zipcodes <gbz>`.
- **map**: `vg giraffe -p -b hifi -t 16 <giraffe_args> -Z <gbz> -d <dist> -m <min> -z <zipcodes> -G <reads>`,
  then `vg annotate -m --search-limit=-1` (positions on reference paths),
  `vg gamcompare --range 200 -T` against the truth GAM, and `vg stats -a`.
  With `--reps N` it runs N times into `rep1..N/` (replicates, e.g. for the determinism
  test); replicates 2..N
  delete their GAMs after their statistics are written. `--run-tag TAG` writes the
  runs to `<name>/TAG/rep1..N/` instead (kept apart from the default runs, which
  `collect` reads), and `--keep-gams` keeps the GAMs of every replicate.
- **density**: giraffe with `--track-provenance` on the read set; for each read,
  `stage.minimizer.results` (seeds the scheme picks in the read, before giraffe's
  filters) and `stage.seed.results` (graph hits of the minimizers giraffe keeps).
  Run it on every read set: `collect` only uses the density run of the same read set.
- **filters**: per giraffe funnel filter, how many minimizers pass and fail, plus
  seeds per read (median, mean, 90th percentile).
- **fingerprint**: `vg view -aj mapped.gam | jq` keeps, per read, the name, MAPQ,
  score, identity, number of mappings, start (`node:offset:strand`) and a hash of the
  path (nodes, offsets, strands and edit lengths; no read bases). The job checks that
  there is one row per alignment and no duplicate read names.
- **determinism**: compares the fingerprints of the `TAG` replicates **by read name**
  (threaded Giraffe writes reads in any order), and separately together with the
  earlier `rep1/` run if it has a fingerprint. Each read is `identical`, `score_only`
  (same path, different MAPQ/score/identity), `path_diff_same_start`, `start_diff` or
  `mapped_flip`; `correct_flip` counts reads whose gamcompare verdict differs between
  runs. `summary.tsv` has one row per comparison; `reads.tsv` lists only reads that
  differ, one row per run, with the node it ran on.
- **timing**: speed on fixed hardware. One job reserves a whole node (`--exclusive`,
  `--exclude=<exclude_nodes>` written into the script) and runs only giraffe (16
  threads) for the `--only` conditions, 3 rounds, reloading the indexes every run.
  Each round starts one place further along the `--only` order (A B C REF, then
  B C REF A, then C REF A B), so no condition is always first or last. Each GAM is
  deleted once written. After every run a row is appended to `runs.tsv`: `job_id`, node,
  round, order, condition, `load_seconds`, `map_seconds`, `reads_per_sec_per_thread`,
  `reads_per_cpu_second`, `cpu_seconds`, `memory_gb` (all from `giraffe.log`).
  `--submit` submits the same script twice; two jobs running at the same time are on
  different nodes, one queued behind the other may land on the same node (see `node`).
- **timing --real** (`--reads` is a real-read FASTQ, mapped with `-f`; the lab's practice
  for real-read speed is the full read set): 64 threads, partition `long`, 300000 MB,
  30 h, one round per job; two scripts, `_order1` in the `--only` order and `_order2`
  reversed, each submitted once. After each giraffe run and before its GAM is deleted,
  `vg filter -T "name;is_aligned;mapping_quality;score;identity;nodes"` writes
  `fingerprint.tsv` (per read: aligned, MAPQ, score, identity, first and last node with
  strand, number of nodes; no bases), and `real_stats.tsv`, `mapq_hist.tsv` and
  `identity_hist.tsv` are computed from it (MAPQ and identity over aligned reads). The job
  stops if the fingerprint does not have exactly as many reads as giraffe reports.
  `runs.tsv` gets the extra columns `peak_rss_gib`, `reads`, `aligned_fraction`,
  `mapq0_fraction`, `mapq60_fraction`, `identity_mean`, `identity_median`. At the end the
  job compares its runs read by read into `placement.tsv` (every method against the
  `minimizer` label of its density group). Read set name: the file name without
  `.fastq.gz`, `.fq.gz`, `.fastq` or `.fq`. Real-read (64 threads) and simulated-read (16 threads) speeds are
  reported separately and are not comparable.
- **timing, repeated on shared nodes**: `--no-exclusive` leaves out `--exclusive` and
  `--exclude` (and then needs no `exclude_nodes` in the config); `--rounds R` sets the
  rounds per job (default 3); `--jobs N` writes N scripts `timing_<reads>_job1..N.sbatch`,
  job k starting its rounds k-1 places further along the `--only` order, so over the jobs
  every condition runs first and last. `--jobs` refuses `--submit`: submit the scripts by
  hand in `--dependency=afterany` chains to limit how many run at once. Without these
  options the generated scripts are unchanged. On shared nodes wall-clock times depend on
  the neighbours, so only ratios within a job are compared, over many jobs
  (`prototype/timing_repeats.py`). These jobs write to `results/<reads>/timing_repeats_runs/`,
  not `timing/`, so `collect`'s timing tables hold the exclusive jobs only.
- **timing --real --sampled**: `--reads` is a sample of real reads (e.g. 100k from
  `prototype/sample_fastq.py`), mapped like simulated reads (16 threads, the rotated
  rounds, map memory, partition `medium`) with `--real`'s fingerprints, statistics and
  `placement.tsv` (from round 1). Works with `--no-exclusive`, `--jobs` and `--rounds`;
  full `--real` runs do not take `--jobs` or `--rounds`.

Typical order: index → for each read set: map and density → filters → collect.
To test mapping determinism: `map --run-tag det --reps 3 --keep-gams --fixed-nodes`
→ `fingerprint --run-tag det` (and `fingerprint` for an earlier `rep1/`) → `determinism --run-tag det`.

## Outputs

### results.tsv (one row per condition with map results for that read set)

| Column | Meaning |
|---|---|
| `name`, `label`, `density_group` | from `conditions.tsv` |
| `reads_set` | read set name (truth GAM file stem) |
| `index_args` | from `indexes/<name>.params`, i.e. what the index was actually built with |
| `giraffe_args`, `git_commit`, `run_date` | from the map run's `run_info.tsv` (commit with `-dirty` if `seed_eval.py` or `conditions.tsv` had uncommitted edits); runs from before `run_info.tsv` existed have commit `NA` and the gamcompare file time as an approximate date |
| `total_reads`, `aligned` | `vg stats -a` |
| `eligible` | reads whose truth has a position on a reference path (others cannot be scored) |
| `correct`, `wrong`, `accuracy` | `gamcompare`: mapped within 200 bp of the truth on the same path; `accuracy = correct / eligible` |
| `reads_per_sec_per_thread` | giraffe's `Mapping speed` (excludes index loading); median over replicates |
| `speed_min`, `speed_max`, `reps`, `reps_consistent` | spread and number of replicates; `reps_consistent` = all replicates got the same number of reads right |
| `giraffe_mem_gb` | giraffe's reported memory footprint |
| `giraffe_wall_s` | wall time of the giraffe command, mostly index loading for small read sets |
| `index_min_gb`, `index_zipcodes_gb`, `index_build_min`, `index_peak_gib` | index sizes, build time and peak memory |
| `occurrences` | `Minimizer occurrences` reported by `vg minimizer` |
| `density_graph` | `occurrences / graph length` (graph length from `prepare_graph_length`) |
| `density_reads` | **main density measure**: seeds picked per bp of read (`stage.minimizer.results / length`) |
| `seeds_per_read_median`, `seeds_per_read_mean`, `seeds_per_read_p90` | seeds per read (`stage.seed.results`, graph hits of the kept minimizers): median, mean and nearest-rank 90th percentile |
| `kept_minimizers_per_read` | minimizers left after giraffe's filters, per read |
| `mean_hits_per_kept_minimizer` | `stage.seed.results / kept minimizers`: how many graph hits each kept minimizer has on average |
| `timed_jobs` | number of exclusive-node `timing` jobs (`results/<reads>/timing/`) this condition ran in (`NA` without one) |
| `timed_speed_median` | median over those jobs of the per-job median speed (reads/s/thread); descriptive only |
| `timed_ratio_min`, `timed_ratio_max` | range over those jobs of the reads/s/thread ratio against the `reference` condition (giraffe default) **in the same job**; `NA` in jobs without it |
| `seeds_per_read_p99`, `seeds_per_read_max` | high end of seeds per read: nearest-rank 99th percentile and maximum (the maximum is one read, so it is noisy); appended after the older columns, which keep their positions |

### timing.tsv (exclusive-node speed, written when timing runs exist)

One row per condition and timing job: `node`, `runs`, `speed_median`/`speed_min`/`speed_max`
(reads/s/thread), `ratio_vs_ref` (this median ÷ the reference condition's median in the
same job), and medians of reads per CPU-second, load seconds and map seconds, plus the
largest giraffe memory footprint. Runs that did not finish are left out. Compare speeds
within a job; `--snapshot` copies this file too.

These exclusive-node numbers (`timing.tsv` and the `timed_*` columns of `results.tsv`) use
reads/s/thread against the giraffe default, a different definition from
`prototype/timing_repeats.py` (load, map and reads per CPU-second against the density-matched
minimizer, over repeated shared-node jobs). Speed conclusions are based on `timing_repeats`.

`timing_runs.tsv` concatenates every job's `runs.tsv` (one row per run). For real reads,
`placement.tsv` concatenates the jobs' placement tables: per job and pair, `reads`,
`aligned_both`, `aligned_only_condition`, `aligned_only_against`, `different_start` (reads
aligned by both whose first node or strand differs), and among those
`higher_score_condition`, `higher_score_against`, `equal_score`. A real read set has no
truth, so `collect` writes no `results.tsv` or `mcnemar.tsv` for it (only the timing tables).

### Plots

`comparison.png` has one bar per condition: mismapped reads as % of reads with a truth
position on CHM13 (`eligible`) and seed hits per read (graph positions)
on the same read set (solid bar: median, light bar: mean). There are no error bars on
accuracy: every condition maps the same reads, so differences are judged with the
paired test in `mcnemar.tsv`. With at most 6 conditions each bar is labelled
`wrong/eligible` and, for non-reference conditions, the difference in mismapped reads
against the reference with its Bonferroni-corrected p, e.g. `−12 wrong vs ref (p=1)`
(negative: fewer mismapped reads than the reference); with more conditions the bars
carry no labels. Conditions are grouped by method (`label`), densest first within a
method, with the `reference` group last; each method has its own color and the
reference is gray. There is no speed panel: speeds are compared with the repeated timing
runs (see "Note on speed"); `results.tsv` keeps the single-run speed columns for reference. Seed columns are only filled from a density run on
the same read set; without one they are `NA`, and if no condition has one the seeds
panels are left out. A third panel shows the high end of seeds per read: a bar for the
99th percentile and the maximum as a number above it (vertical when the bars carry no
other labels). It is a separate panel because p99 and the maximum are one to two orders
of magnitude above the median and would flatten the median and mean bars.

`density_sweep.png` is drawn only when the read set has results at two or more
density groups. Its x axis is the measured density on reads (`density_reads`); the
four panels show mismapped reads (% of reads with a truth position on CHM13), seed
hits per read (mean solid, median dashed), the 99th percentile of seed hits per read
(its own panel, for the same scale reason) and giraffe memory. Each method is one line
and reference conditions are gray points.

### mcnemar.tsv (paired accuracy comparison)

`collect` compares conditions read by read on the same eligible reads, using each
condition's `rep1/compare.tsv` (or its single run). Tested pairs: every condition
against the `reference` group, and within each density group every other method
against the `minimizer` label. A `minimaxmer` condition is also tested against every other
method of its group. The `family` column says which kind of pair a row is (`vs_reference`
or `within_<group>`). Only discordant reads matter: `only_condition_correct`
(b) and `only_against_correct` (c). `p_exact` is the two-sided exact McNemar test
(b ~ Binomial(b + c, 1/2) if the two conditions are equally accurate);
`p_bonferroni = min(1, p_exact × m_tests)`, with `m_tests` the number of tests in this
read set's table. `net_reads = b − c` (positive: the condition gets more reads right
than `against`), `net_percent` is the same as % of `eligible_both`.
`wrong_condition` and `wrong_against` are the mismapped counts on the shared reads.
`collect` stops if a condition's `compare.tsv` and `results.tsv` disagree on the number
of mismapped reads.

### Filter summaries

`filter_summary.tsv` lists, per condition and giraffe funnel filter, the number of
items passing and failing; the `minimizer` stage filters act on minimizers.
`minimizer_summary.tsv` gives picked, kept and hit counts and seeds per read
(median, mean, 90th percentile) per condition.

### Note on speed

Mapping speed depends on how busy the compute node is, so it varies between runs. Speed
comparisons are taken from the repeated timing runs (`timing --no-exclusive --jobs N`,
summarised by `prototype/timing_repeats.py`: ratios within each job against the
density-matched minimizer, distributed over jobs). The speed columns of `results.tsv`
(`reads_per_sec_per_thread`, `speed_min`/`speed_max`) come from single map runs and are for
reference only. The exclusive-node timing runs (`timing.tsv`, `timing_runs.tsv`) were all made
on one node, `phoenix-17`.

To measure speed this way:

```bash
python3 seed_eval.py timing --reads READS.gam --only A B REF --no-exclusive --jobs 10
python3 prototype/timing_repeats.py data_dir/results/<reads>/timing_repeats_runs JOB_ID ... --ref REF
```

The first command writes 10 jobs (submit them yourself). Each job maps with every condition
in turn. The second command gives each condition's time divided by REF's time in the same job.

## Minimaxmer prototype (`prototype/`)

`prototype/minimaxmer.py` is a Python prototype of minimaxmer
selection (Yoshi Suzuki's [himapper_rp](https://github.com/yoshihikosuzuki/himapper_rp)):
in each window of `w` k-mers it picks the minimum of
`(min(genome_count, 65534), 65535 - reads_count, hash >> 32)`, i.e. a k-mer rare in the
genome and common in the reads (step 1), and counts the windows each one wins (step 2).
`--rspace` adds Yoshi's step 3 (sparse thinning). The k-mer encoding, hash (`wang_hash_64`),
canonical orientation and tie rule follow gbwtgraph, so with constant counts the selection
is vg's `(w, k)`-minimizers (up to ties in the high 32 bits of the hash). As in Yoshi's
code, k-mers absent from the genome (count 0) are not candidates.

With `--counts-mmc FILE` (a `minimax_counts` table), both counts come from that file, and
k-mers not in it are not candidates. `--reads-mode bins` groups `reads_count` into 0 / 1-9 /
10+, and `--reads-mode constant` ignores it. Without `--counts-mmc`, `genome_count` is 1 for
every k-mer and `reads_count` is counted from the input reads (`run_on_reads.sh` uses this).

`prototype/test_minimaxmer.py` tests it on synthetic sequences. `prototype/vg_oracle.cpp`
is a small program that prints gbwtgraph's own k-mer codes, hashes and minimizers; build
it against the vg tree (command in its header) and pass it with `--oracle` to compare
the prototype with gbwtgraph directly:

```bash
python3 -I prototype/test_minimaxmer.py --oracle /path/to/vg_oracle
```

On real reads, `prototype/run_on_reads.sh` (an sbatch script, submitted from the repo
directory) feeds `name<TAB>sequence` from a GAM through `vg view -aj | jq` into the
prototype twice: with constant counts (must equal vg's minimizers) and with reads
counts from those reads. Sequences only pass through pipes; outputs are per-read
counts (read name, length, selected) and summaries. `prototype/compare_with_vg.py`
matches these per read **by name** with a density run's `per_read.tsv`, whose
`minimizers` column is Giraffe's `stage.minimizer.results`:

```bash
sbatch -o LOG prototype/run_on_reads.sh TRUTH.gam OUT_DIR /path/to/vg
python3 -I prototype/compare_with_vg.py OUT_DIR/per_read_constant.tsv \
    data_dir/results/<reads>/<condition>/density/per_read.tsv OUT_DIR/per_read_counted.tsv --out comparison.tsv
```

### Design: no step 3 in vg

Yoshi runs all three steps (window minimum, windows won, uniformly sparse thinning)
independently on the reference and on the reads, and keeps the k-mers selected on
both. Our vg version (`vg_minimax`) implements step 1 only, on both sides: the
index stores every window minimum (vg's rule, every tied minimum is selected), and
giraffe selects the read's window minima the same way, then applies its own filters
(window downsampling and the `--max-min` / `--num-bp-per-min` cap with `-b hifi`)
instead of step 3.

Step 3 cannot be done on the graph side the way Yoshi does it on a linear reference:
gbwtgraph does not select minimizers along whole haplotypes but on short haplotype
windows around each node (a node plus `k + w - 2` bases of each haplotype
continuation). That gives the same result as whole haplotypes for any rule that depends
on one window only, but Yoshi's thinning accepts seeds greedily by priority and blocks
their neighbours, so a decision can depend on seeds far away and cannot be made
consistently on those short windows. A possible later experiment is step 3 on the read
side only: it would only pick among the read's step-1 window winners (one of the tied
minima, by Yoshi's robust winnowing), which are winners of the same window in the
graph, so in an error-free window that exists on a graph haplotype every kept seed would
still be in the index. The Python prototype implements step 3 (`--rspace`) for such
comparisons.

### Repeated speed measurements

- `python3 -I prototype/sample_fastq.py --total N --n 100000 --seed S --names NAMES < FASTQ > OUT`:
  a uniform sample of exactly n reads without replacement (`random.Random(S).sample(range(N), n)`),
  kept in input order; stops unless the input has exactly N four-line records. Writes the
  sampled read names (no sequences) to NAMES.
- `python3 -I prototype/timing_repeats.py results/<reads>/timing_repeats_runs JOB_ID ... [--ref mini_k31w20] [--out PREFIX] [--plot PREFIX]`:
  for the given jobs only (those Slurm reports as completed), each condition's median over rounds of `load_seconds`,
  `map_seconds` and `reads_per_cpu_second`, divided by the reference's in the same job;
  across jobs the median, quartiles, range, and the fraction of jobs with a ratio above 1.
  For real reads also aligned fraction, MAPQ and identity over all runs, the pooled MAPQ
  histogram, and `placement.tsv` summed over jobs. The same ratios against each condition's density-matched
  minimizer (`group_ratio_*`), `memory_gb` median and range (`PREFIX_memory.tsv`), and with `--plot` one
  figure per density group (`PREFIX_d2-11.png`, ...).

### Analyses of minimaxmers against the 2/11 group

These scripts compare methods read by read (numbers only):

- `prototype/stage_work.py RESULTS_DIR METHOD ...`: work per read in each giraffe stage
  (seeds, zipcode trees, chains, alignments, aligner invocations and bases), from
  `density/stage_counts.jsonl`, which `sbatch -o LOG prototype/stage_counts.sh RESULTS_DIR VG
  CONDITION ...` extracts from the density provenance with `vg view -aj | jq` (stage `results`
  counts and the aligner's per-read totals; no timings). Median, mean, p99, total, and an exact sign test against
  the first method. It helps explain speed differences but does not replace measured
  speed: minimaxmers give fewer seeds yet map more slowly, because the cost of looking up
  the count table is not part of these counts.
- `prototype/hard_reads.py RESULTS_DIR METHOD ... --selector mini_k31w20`: accuracy and exact
  McNemar tests on the hard reads, the 10% with the most seed hits under a condition that
  is not compared (so the choice favours none), and, where every method has a density run,
  the 10% with the most seed hits averaged over the methods after the first (`group_top10`),
  and the reads whose seed hits differ between the first method and any other (`hits_differ`); `--deciles PLOT` plots the error rate per decile of
  the selector's seed hits.
- `prototype/mapq_pr.py RESULTS_DIR METHOD ...`: MAPQ distribution, reads kept and wrong at
  MAPQ thresholds, a paired test on reads wrong at MAPQ 60, and the precision-recall curve
  of the lab's `plot-pr.R`.
- `prototype/mapq60_reads.py MAPQ60_DIR ANNOTATION_DIR METHOD ...`: one row per read wrong at MAPQ 60 in any
  method, from re-mapping just those reads with `--track-provenance` (data dir
  `jobs/mapq60_provenance_*.sbatch`): truth span, each method's mapped span, distance, MAPQ, seeds, chains,
  whether the read went to the other copy of a segmental duplication, and the T2T censat classes
  (`chm13v2.0_censat_v2.1.bed`) and SD bases (`chm13v2.0_SD.bed`) under the truth span; `caveat` marks
  truth on chrY and alignments without a CHM13 position.

Their outputs are in `results_summary/HG002-sim-hifi-10k/` (`stage_work.tsv`,
`hard_reads.tsv`, `mapq.tsv`, `pr.png`) and `results_summary/HG002-sim-hifi-100k/`
(`hard_reads.tsv`, `mapq.tsv`, `pr.png`, `error_by_seed_decile.png`, `mapq60_reads.tsv`, `mapq60_summary.tsv`).

### Development checks

Three small scripts were used only while developing and are not part of any stage or analysis
(numbers only):

- `prototype/compare_gam_json.py A.json B.json`: compares two `vg view -aj` outputs read by
  read, ignoring timing fields (used to check that two runs or builds give the same alignments).
- `prototype/density_table.py --baseline DENSITY_SUMMARY.tsv MODE=PER_READ.tsv:SUMMARY.txt ...`:
  a table of the prototype's measured density per counting mode against a baseline condition's
  density run (the 1k density comparison of minimaxmer count modes).
- `prototype/kff_header.py FILE.kff`: prints the header fields of a KFF file (k, max, min,
  data_size), to check the KMC output before building the count table.

### Count table (`prototype/minimax_counts.cpp`)

The counts minimaxmer selection needs, for every canonical k-mer in the graph:
`genome_count` = number of distinct graph positions over all haplotypes (the count
`vg minimizer --threshold` uses for weighted minimizers, and the graph analogue of
Yoshi's occurrences in the genome), and `reads_count` = the k-mer's count in a KFF file
(KMC on the sample's reads; 0 if absent). Both saturate at 65535; canonical means the
smaller 2-bit key of the two strands, as in gbwtgraph. Build it against the vg tree
(command in its header) and run:

```bash
minimax_counts -g GRAPH.gbz -f READS.kff -o PREFIX -c COVERAGE -t THREADS
```

It writes `PREFIX.mmc` (24-byte header `MMXCNT01`, k, n; then n sorted 64-bit keys and n
32-bit values `genome_count << 16 | reads_count`, little-endian), `PREFIX.summary.tsv`
(first: distinct graph k-mers, fraction with `genome_count` 1, the `reads_count` peak
among `genome_count` 1 k-mers with `reads_count` ≥ 5 and its ratio to `COVERAGE`, and the
counts file size), and histograms of each count plus a joint table (numbers only).
The k of the KFF sets k. For the reads, count with
`kmc -k31 -ci2 -cs65535 -okff` (16-bit counts like Yoshi's; k-mers seen once dropped).
`prototype/test_minimax_counts.py` checks the output on a small graph against a brute
force over the haplotype walks (genome counts) and a `kmc_dump` of the same reads
(reads counts), and reads the counts file back.

`results_summary/HG002Revio_full_k31_counts/` holds the outputs for our graph and the full
HG002 Revio HiFi reads (12.8 M reads, 66x haploid coverage), without the 31.8 GB counts
file itself. They were made with
`kmc -k31 -m128 -ci2 -cs65535 -okff -t16 -hp READS.fastq.gz OUT TMP` and
`minimax_counts -g GRAPH.gbz -f OUT.kff -o PREFIX -c 66.0 -t 32`:

| File | Contents |
|---|---|
| `summary.tsv` | `minimax_counts` summary (headline numbers first) |
| `genome_hist.tsv` | number of graph k-mers per `genome_count` |
| `reads_hist.tsv` | number of graph k-mers per `reads_count` (all graph k-mers) |
| `joint.tsv` | graph k-mers by binned `genome_count` × binned `reads_count` |
| `reads_hist_genome1.tsv` | `reads_count` histogram of the k-mers with `genome_count` 1, from `python3 -I prototype/plot_reads_count.py hist PREFIX.mmc OUT.tsv` |
| `reads_count_genome1.png` | that histogram, with the coverage and half of it marked: `python3 -I prototype/plot_reads_count.py plot reads_hist_genome1.tsv OUT.png --coverage 66` |

## Minimaxmer patches for vg (`patches/`)

`vg_minimax` is vg with two patches: `patches/vg/` against vg commit `c91cc1d91`
(`vg minimizer` and `vg giraffe --minimaxmer-counts`), and `patches/gbwtgraph/` against
GBWTGraph commit `7fac0af` (minimaxmer selection, `include/gbwtgraph/minimax.h`), which is
the GBWTGraph commit vg `c91cc1d91` pins as its `deps/gbwtgraph` submodule. Each directory
has the license of the project it patches (`LICENSE`). To build:

```bash
git clone https://github.com/vgteam/vg.git vg_minimax && cd vg_minimax
git checkout c91cc1d91
git submodule update --init --recursive        # deps/gbwtgraph is now at 7fac0af
git apply /path/to/seed-eval/patches/vg/0001-*.patch
git -C deps/gbwtgraph apply /path/to/seed-eval/patches/gbwtgraph/0001-*.patch
make -j 16                                     # vg's own build instructions apply (see vg's README)
```

`git am` instead of `git apply` keeps the patches' commit messages. The vg patch leaves the
submodule pointer unchanged, so `git status` in the vg checkout shows `deps/gbwtgraph` as
modified. The counts file is written by `prototype/minimax_counts.cpp` (see "Count table").

## Safeguards

- **Index provenance.** Each index job writes `indexes/<name>.params` (vg binary and
  version, `index_args`, graph and distance index). `index`, `map`, `density` and
  `timing` refuse to generate scripts, and the jobs refuse to start, if it disagrees with the
  row. `index` deletes the old `.params` before rebuilding, so a half-built index is
  never used.
- **Program and counts checksums.** `vg_minimax` reports the same `vg version` as
  `vg_baseline`, so the version line cannot tell them apart. Each index job also
  writes `indexes/<name>.vg_sha256` (path, size and sha256 of the vg binary) and, for a
  row whose `index_args` has `--minimaxmer-counts FILE`, `indexes/<name>.counts_sha256`
  for that counts file. These are separate files so that `.params`, which queued jobs
  compare word for word, never changes. `map`, `density` and `timing` refuse to
  generate scripts without them, and the jobs recompute the sha256 before mapping and
  stop if the binary or counts file is not the recorded one. For indexes built before
  these records existed, `python3 seed_eval.py checksums [--only NAME ...]` writes
  them, but only for files last modified before the index's `.params` was written
  (otherwise rebuild the index); it keeps existing records unless `--force`.
- **Atomic outputs.** The main outputs (index files, GAMs, `compare.tsv`, `per_read.tsv`,
  the `filters` summaries, fingerprints) are written as `.tmp` and renamed when complete.
  Small tables are written directly; `.params` is written last.
- **Run provenance.** Map and density jobs record commit, date, node, job id and
  giraffe arguments in `run_info.tsv`. Commit before submitting so the commit is not `-dirty`.

## Notes and known caveats

- **Giraffe caps the minimizers it uses.** With `-b hifi`, giraffe downsamples
  minimizers in windows of `min(read length / 15, 227)` bp and then keeps at most
  `max(79, read length / 152)` (`--max-min`, `--num-bp-per-min`). The current
  sweep keeps the preset fixed; these are command-line options and could be varied
  through `giraffe_args`.
- Small read sets give few errors; use 100k+ reads to compare accuracy.
- **The `hifi` preset turns giraffe's explored MAPQ cap off** (`explored-cap false` in
  vg's `giraffe_main.cpp`). That cap is the only place giraffe treats a minimizer's
  hash as a uniform random number (`get_prob_of_disruption_in_column`), so a scheme that
  changes the hash values (weighted minimizers, minimaxmers) breaks no
  assumption of giraffe's MAPQ under `-b hifi`; with presets that use the cap, that
  MAPQ cap would be computed from non-uniform hashes.
- **Mapping was deterministic in our test** (`results_summary/HG002-sim-hifi-100k/determinism_summary.tsv`):
  for `sync_k31s11` and `default_k31w50W` on the 100k set, three runs plus one run five
  days earlier gave identical fingerprints for all 100,000 reads. Scope: 2 conditions,
  100k reads, 384-core nodes only, 16 threads; other conditions, node types or thread
  counts were not tested.
- Indexes are large (~24 GB `.min` at density 2/21, ~12 GB at w50).
- Changing the graph means pointing `graph` in the config at it, re-running `prepare`
  and rebuilding every index (the `.params` check refuses indexes built on another
  graph path); changing only the reads needs `map`/`density`/`filters`/`collect` with
  `--reads`.

## Known limitations

- **chrY truth is not excluded.** The lab's `gamcompare` runs use `--ignore 'CHM13#0#chrY'`;
  ours do not. On the 100k set, 354 reads have truth positions only on `CHM13#0#chrY`; every
  method gets ~19-20% of them wrong (~10-11% of all wrong reads). Without them the wrong
  reads are minimaxmer 578 against 579 / 598 / 580 for `mini_k31w10` / `sync_k31s21` /
  `miniW_k31w10`, and every McNemar conclusion stays the same (no test below p = 0.11).
  The read list is `results/HG002-sim-hifi-100k/truth_chrY_reads.tsv` in the data directory.
- **Real reads overlap the reads the minimaxmer counts were made from.** The count table
  (`HG002Revio_full_k31_graph.mmc`, `reads_count`) is built from all 12,800,930 HG002 Revio
  reads. The random real-read sample `HG002Revio-random100k-s20261027` is 100,000 of those
  reads (about 0.8%), so a small part of `reads_count` comes from the reads being mapped; the
  effect is expected to be small. The full real-read timing (`timing --real`) maps all of
  them. The simulated reads are not affected: they were simulated with pbsim2 from HG002's
  haplotype sequences (extracted from the graph `hprc-chm-hg002-v2.0.full` with HG002 as the
  reference sample), with read lengths and quality strings taken from 100k HG002 Revio reads
  (`--sample-fastq`, `--depth 4 --length-min 10000 --difference-ratio 6:50:54`, per
  `reads/sim/hifi/HG002/README.txt` and vg's `scripts/make_pbsim_reads.sh`). The count table
  does not contain their sequences; that counts and simulated reads come from the same person
  is the design (one index per sample), not a leak.
- **No centromere-specific accuracy.** There is no list of centromeric HG002 simulated reads.
  The lab's centromeric list is for CHM13 simulated reads
  (`reads/sim/hifi/CHM13/CHM13-sim-hifi.centromeric.txt`, 192,312 names), but those reads
  are no longer where their README puts them. Even with them, minimaxmers could not be run
  on CHM13 reads, because the counts are HG002's.
- **Minimaxmers lack Yoshi's step 3 (uniformly sparse thinning).** `vg_minimax` implements
  step 1 on both sides (see "Design: no step 3 in vg"). On the graph, minimizers are
  chosen on short node-local haplotype windows, while step 3 accepts seeds greedily and
  blocks their neighbours, so its decisions can depend on seeds far away. They cannot be made
  consistently on those windows. Step 3 on the read side only is a possible later
  experiment. As a result, `minimax_k31w20` has about twice the density of `mini_k31w20`
  (about 2/11), so it is compared with the 2/11 group and not with Yoshi's sparser selection.
- **Minimaxmers leave about 0.5% more real reads unaligned.** On the random 100k real reads,
  501 reads are aligned by `mini_k31w10` but not by `minimax_k31w20` (28 the other way). All 495
  of them with a CHM13 position lie in centromeric satellite arrays (`chm13v2.0_censat_v2.1.bed`:
  mostly HOR, then hsat2, hsat1A, hsat3), and `mini_k31w10` aligns them mostly with low MAPQ
  (3 of 501 at MAPQ 60) and identity around 0.95. 372 of them are aligned by both
  `sync_k31s21` and `miniW_k31w10`, so the loss is mostly specific to minimaxmers. A guess, not
  tested: in these arrays minimaxmers pick k-mers that occur in only one copy in the graph and
  not in HG002's copy. On the simulated reads, the extra unaligned reads are almost all
  reads without a truth position on CHM13, so accuracy is not affected. Tables:
  `results_summary/HG002Revio-random100k-s20261027/unaligned_reads/`.
