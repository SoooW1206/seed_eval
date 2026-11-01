# seed-eval

Author: Sophia Tianxiang Wang (rotation in the Paten lab, October-November 2026)

seed-eval compares seed selection schemes in `vg giraffe` on long reads: minimizers, closed
syncmers, weighted minimizers and minimaxmers. For each scheme it builds a seed index, maps
simulated reads with known positions, and reports accuracy, speed, memory and seeds per read.
Details of every stage, column and output: [docs/DETAILS.md](docs/DETAILS.md).

## Using seed-eval

Needs `vg`, Slurm, Python 3.9+ with matplotlib, and `jq`. Copy `config.example.tsv` to
`config.tsv` and fill in the paths. Put one row per seed index in `conditions.tsv`. Then run
the stages in order, waiting for each stage's jobs to finish:

```bash
python3 seed_eval.py prepare --submit                    # distance index, once per graph
python3 seed_eval.py index --submit                      # seed index per condition
python3 seed_eval.py map --reads READS.gam --submit      # map, then score against the truth
python3 seed_eval.py density --reads READS.gam --submit  # seeds per read
python3 seed_eval.py filters --reads READS.gam --submit  # giraffe filter statistics
python3 seed_eval.py collect --reads READS.gam --snapshot  # tables and plots
```

## Minimaxmers

Minimaxmers come from Yoshihiko Suzuki's
[himapper_rp](https://github.com/yoshihikosuzuki/himapper_rp). They pick seeds by k-mer counts
instead of by a random hash. Each k-mer has a `genome_count` (how often it occurs in the
genome) and a `reads_count` (how often it occurs in the sample's reads). The algorithm has
three steps:

1. **Window minimum.** In each window of `w` k-mers, pick the k-mer that is rarest in the
   genome; among those, the one most common in the reads; the hash breaks any tie. K-mers
   that are not in the genome are never picked, so most k-mers with read errors are skipped.
2. **Windows won.** Count how many windows each picked k-mer wins.
3. **Thinning.** Keep the picked k-mers, most windows won first, and drop the others near
   each kept one, so the seeds are spread out evenly.

Yoshi runs these steps on the reference and on the reads, and keeps k-mers picked on both.
In vg (`vg_minimax`, built from `patches/`), the graph and the reads use step 1 only. The
genome is the graph; `genome_count` is the number of graph positions of the k-mer. vg picks
minimizers in short windows around each node, but thinning can depend on seeds far away, so
vg has no thinning. `-w 20` minimaxmers then select about as many k-mers per base as `-w 10`
minimizers (density about 2/11). `prototype/minimaxmer.py` has all three steps.

To run them: count the reads' k-mers with KMC, build the count table with
`prototype/minimax_counts.cpp`, build `vg_minimax`, and give the same counts file to both
`vg minimizer` and `vg giraffe` with `--minimaxmer-counts` in `conditions.tsv` (details: [count
table](docs/DETAILS.md#count-table-prototypeminimax_countscpp), [patches](docs/DETAILS.md#minimaxmer-patches-for-vg-patches)).

## Example results

100k simulated HG002 HiFi reads on the HPRC v2.0 CHM13 graph, haplotype-sampled for HG002
with HG002 held out, `-b hifi`.
Wrong reads out of 93,392 scored: minimizer `mini_k31w10` 645, minimaxmer `minimax_k31w20` 650,
closed syncmer `sync_k31s21` 668, weighted minimizer `miniW_k31w10` 645, lab default
`default_k31w50W` 653. Minimaxmers took about 1.35 times as long to map as `mini_k31w10`
and used 99 GB of memory instead of 72 GB. All tables: [results_summary/](results_summary/).

## Acknowledgments

- Minimaxmers were proposed by Yoshihiko Suzuki
  ([himapper_rp](https://github.com/yoshihikosuzuki/himapper_rp)). The vg Giraffe
  implementation (`patches/`) implements the method of himapper_rp (commit 2d4e0fa).
  `window_orientation`, `robust_window_counts` and `sparse_thinning` in `prototype/minimaxmer.py`
  are translated from himapper_rp's `minimaxmer.c` and carry its BSD license notice
  (`prototype/LICENSE.himapper_rp`).
- vg ([vgteam/vg](https://github.com/vgteam/vg)) and GBWTGraph
  ([jltsiren/gbwtgraph](https://github.com/jltsiren/gbwtgraph)) are MIT-licensed; the patches
  are against vg commit `c91cc1d91` and GBWTGraph commit `7fac0af`.
- k-mers in the reads were counted with [KMC](https://github.com/refresh-bio/KMC).
- The evaluation follows the conventions of the Paten lab's long-read Giraffe experiments
  (read sets, `gamcompare` scoring, the precision-recall plot of `plot-pr.R`).

## License

seed-eval is released under the MIT license (`LICENSE`). Two parts follow other licenses:
`window_orientation`, `robust_window_counts` and `sparse_thinning` in `prototype/minimaxmer.py`
are translated from himapper_rp and follow its BSD license notice
(`prototype/LICENSE.himapper_rp`); the patches in `patches/` follow the MIT licenses of vg and
GBWTGraph (`patches/vg/LICENSE`, `patches/gbwtgraph/LICENSE`).
