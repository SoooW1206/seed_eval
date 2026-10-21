# seed-eval

Author: Sophia Tianxiang Wang (rotation in the Paten lab)

This is work in progress; the code and results may change.

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
- Slurm (`sbatch`, `sacct`) for everything heavy. Jobs use the `medium` partition
  with 16 threads; the `filters` job uses `short`.
- Python 3.8+ with matplotlib (only `collect` plots); `jq` on the compute nodes
  (`density` and `filters` check for it).
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
   | `exclude_nodes` | optional: Slurm node list that `--fixed-nodes` passes to `sbatch --exclude`, to keep runs on one node type |

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
  results_summary/<reads>/      committed snapshots of the summary tables and plots

data_dir/                       large files, not in git
  graph.paths_meta.tsv, graph.length.txt   outputs of prepare
  indexes/<name>.min|.zipcodes|.params
  results/<reads>/<name>/       per-condition map outputs, logs, density/ subdirectory
  results/<reads>/              results.tsv, comparison.png, filter/minimizer summaries
  jobs/                         generated sbatch scripts
  logs/                         Slurm logs
```

## One-time preparation per graph (`prepare`)

`python3 seed_eval.py prepare [--submit]` writes three jobs to `data_dir/jobs/`. They
do not depend on the reads; run them once for a new graph.

| Job | What it does | Used by |
|---|---|---|
| `prepare_check_paths` | `vg paths -M` into `graph.paths_meta.tsv`; confirms the truth's reference (CHM13 here) is a `REFERENCE` path | sanity check for `gamcompare` |
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

## Stages

Run from the repo directory. Every stage except `collect` only **writes** sbatch
scripts to `data_dir/jobs/` and prints them; add `--submit` to submit.
`--only NAME ...` restricts to some conditions, `--conditions FILE` uses another table,
`--reads GAM` picks the truth read set (default: 1k). On `--submit`, `--fixed-nodes`
adds `--exclude=<exclude_nodes>` and `--sbatch-args "..."` adds any other arguments to
the `sbatch` command line; the generated scripts do not change.

| Stage | Command | Output |
|---|---|---|
| index | `python3 seed_eval.py index [--force] --submit` | `indexes/<name>.min`, `.zipcodes`, `.params`; skipped if already built with the same parameters |
| map | `python3 seed_eval.py map --reads GAM [--reps N] --submit` | `results/<reads>/<name>/`: `mapped.gam`, `annot.gam`, `compare.tsv`, `gamcompare.txt`, `giraffe.log`, `stats.txt`, `*.time`, `report.tsv`, `run_info.tsv` |
| density | `python3 seed_eval.py density --reads GAM --submit` | `results/<reads>/<name>/density/`: `prov.gam`, `per_read.tsv`, `density_summary.tsv` |
| filters | `python3 seed_eval.py filters --reads GAM --submit` | one job for all conditions: `results/<reads>/filter_summary.tsv`, `minimizer_summary.tsv` (+ per-condition copies in `density/`) |
| collect | `python3 seed_eval.py collect --reads GAM [--snapshot]` | runs locally: `results/<reads>/results.tsv`, `comparison.png`, `density_sweep.png`; `--snapshot` copies them and the filter/minimizer summaries to `results_summary/<reads>/` |
| fingerprint | `python3 seed_eval.py fingerprint --reads GAM [--run-tag TAG] --submit` | one job per run that has a `mapped.gam` (`<name>/rep*`, or `<name>/TAG/rep*`): `fingerprint.tsv` in the run directory |
| determinism | `python3 seed_eval.py determinism --reads GAM --run-tag TAG` | runs locally: `results/<reads>/determinism/summary.tsv`, `reads.tsv` |

What each job runs:

- **index**: `vg minimizer <index_args> -p -t 16 -d <dist> -o <name>.min -z <name>.zipcodes <gbz>`.
- **map**: `vg giraffe -p -b hifi -t 16 <giraffe_args> -Z <gbz> -d <dist> -m <min> -z <zipcodes> -G <reads>`,
  then `vg annotate -m --search-limit=-1` (positions on reference paths),
  `vg gamcompare --range 200 -T` against the truth GAM, and `vg stats -a`.
  With `--reps N` it runs N times into `rep1..N/` (for timing); replicates 2..N
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
| `giraffe_args`, `git_commit`, `run_date` | from the map run's `run_info.tsv` (commit with `-dirty` if the code had uncommitted edits); runs from before `run_info.tsv` existed have commit `NA` and the gamcompare file time as an approximate date |
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

### Plots

`comparison.png` has one bar per condition: mismapped reads as % of eligible with a
95% Wilson interval, mapping speed, and seeds per read on the same read set (solid
bar: median, light bar: mean). Conditions are grouped by method (`label`), densest
first within a method, with the `reference` group last; each method has its own
color and the reference is gray. With more than 6 conditions the bars carry no value
labels. If the read set has no replicates, the speed panel is titled "single run, not
for comparison". Seed columns are only filled from a density run on the same read
set; without one they are `NA`, and if no condition has one the seeds panel is left out.

`density_sweep.png` is drawn only when the read set has results at two or more
density groups. Its x axis is the measured density on reads (`density_reads`); the
three panels show mismapped reads (% of eligible, 95% Wilson interval), seeds per read
(mean solid, median dashed) and giraffe memory. Each method is one line and reference
conditions are gray points. In the error panel the methods are shifted slightly apart
along x so their intervals do not overlap; the shift is for display only.

### Filter summaries

`filter_summary.tsv` lists, per condition and giraffe funnel filter, the number of
items passing and failing; the `minimizer` stage filters act on minimizers.
`minimizer_summary.tsv` gives picked, kept and hit counts and seeds per read
(median, mean, 90th percentile) per condition.

### Note on speed

Mapping speed depends on how busy the compute node is, so it can vary between runs.
Use replicates (`--reps`) on a large read set and treat small differences with caution.

## Minimaxmer prototype (`prototype/`)

Work in progress. `prototype/minimaxmer.py` is a Python prototype of minimaxmer
selection (Yoshi Suzuki's [himapper_rp](https://github.com/yoshihikosuzuki/himapper_rp)):
in each window of `w` k-mers it picks the minimum of
`(genome_count, -reads_count, hash)`, i.e. a k-mer rare in the genome and common in the
reads (step 1), and counts the windows each one wins (step 2). Yoshi's step 3 (uniformly
sparse thinning) is not implemented. The k-mer encoding, hash (`wang_hash_64`), canonical
orientation and tie rule follow gbwtgraph, so with constant counts the selection is
exactly vg's `(w, k)`-minimizers; as in Yoshi's code, k-mers absent from the genome
(count 0) are not candidates. Genome counts are a placeholder (all 1) for now.

`prototype/test_minimaxmer.py` tests it on synthetic sequences. `prototype/vg_oracle.cpp`
is a small program that prints gbwtgraph's own k-mer codes, hashes and minimizers; build
it against the vg tree (command in its header) and pass it with `--oracle` to compare
the prototype with gbwtgraph directly:

```bash
python3 -I prototype/test_minimaxmer.py --oracle /path/to/vg_oracle
```

## Safeguards

- **Index provenance.** Each index job writes `indexes/<name>.params` (vg binary and
  version, `index_args`, graph and distance index). `index`, `map` and `density`
  refuse to generate scripts, and the jobs refuse to start, if it disagrees with the
  row. `index` deletes the old `.params` before rebuilding, so a half-built index is
  never used.
- **Atomic outputs.** Every file is written as `.tmp` and renamed when complete.
- **Run provenance.** Map and density jobs record commit, date, node, job id and
  giraffe arguments in `run_info.tsv`. Commit before submitting so the commit is not `-dirty`.

## Notes and known caveats

- **Giraffe caps the minimizers it uses.** With `-b hifi`, giraffe downsamples
  minimizers in windows of `min(read length / 15, 227)` bp and then keeps at most
  `max(79, read length / 152)` (`--max-min`, `--num-bp-per-min`). The current
  sweep keeps the preset fixed; these are command-line options and could be varied
  through `giraffe_args`.
- Small read sets give few errors; use 100k+ reads to compare accuracy.
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
