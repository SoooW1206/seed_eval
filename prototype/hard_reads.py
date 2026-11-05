#!/usr/bin/env python3
"""Accuracy on hard reads: the reads with the most seed hits.

Usage: python3 -I hard_reads.py READS_DIR METHOD ... --selector NAME [--out TABLE.tsv] [--deciles PLOT.png]

READS_DIR: data_dir/results/<reads>. Correctness comes from each method's compare.tsv (rep1/ when present),
seed hits from density/per_read.tsv. Only reads eligible for every method are counted. Subsets:
  all                 every read
  selector_top10      the 10% of reads with the most seed hits under --selector, a condition that is not one of
                      the compared methods, so the choice favours none of them
  group_top10         the 10% with the most seed hits averaged over the methods after the first (only if they all
                      have density runs)
  hits_differ         reads whose seed hits under the first method differ from any other method (same condition)
For each subset: wrong reads per method, and an exact McNemar test of the first method against each other one
(only_first_correct / only_other_correct). --deciles plots the error rate per decile of the selector's seed
hits, one line per method. Prints numbers only.
"""
import argparse
import csv
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import seed_eval  # noqa: E402  (plot style and method colors)


def compare(reads_dir, c):
    for p in (Path(reads_dir) / c / "rep1" / "compare.tsv", Path(reads_dir) / c / "compare.tsv"):
        if p.is_file():
            with open(p, newline="") as f:
                return {r["read"]: (r["eligible"] == "1", r["correct"] == "1") for r in csv.DictReader(f, delimiter="\t")}
    sys.exit(f"no compare.tsv for {c}")


def seeds(reads_dir, c):
    p = Path(reads_dir) / c / "density" / "per_read.tsv"
    if not p.is_file():
        return None
    with open(p, newline="") as f:
        return {r["read"]: int(r["seeds"]) for r in csv.DictReader(f, delimiter="\t") if r["seeds"] != "NA"}


def mcnemar(b, c):
    """Two-sided exact McNemar p, summed in log space."""
    n = b + c
    if n == 0:
        return 1.0
    logs = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) - n * math.log(2) for i in range(min(b, c) + 1)]
    top = max(logs)
    return min(1.0, 2 * math.exp(top) * sum(math.exp(x - top) for x in logs))


def top_fraction(values, fraction):
    """Reads with the highest values; ties broken by read name so the choice is reproducible."""
    ranked = sorted(values, key=lambda r: (-values[r], r))
    return ranked[:math.ceil(fraction * len(ranked))]


def deciles_plot(path, reads_set, methods, cmp, eligible, hits, selector):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ranked = sorted((r for r in eligible if r in hits), key=lambda r: (hits[r], r))
    groups = [ranked[i * len(ranked) // 10:(i + 1) * len(ranked) // 10] for i in range(10)]
    labels = {"mini": "minimizer", "sync": "closed syncmer", "miniW": "weighted minimizer", "minimax": "minimaxmer"}
    colors = dict(zip(["minimizer", "closed syncmer", "weighted minimizer", "minimaxmer"], seed_eval.METHOD_COLORS))
    fig, ax = plt.subplots(figsize=(9, 5))
    x = list(range(1, 11))
    for m in methods:
        rate = [100 * sum(not cmp[m][r][1] for r in g) / len(g) for g in groups]
        ax.plot(x, rate, color=colors.get(labels.get(m.split("_")[0]), seed_eval.REFERENCE_COLOR), lw=2, marker="o",
                ms=6, markeredgecolor="#fcfcfb", label=m)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{i}\n{hits[g[0]]}-{hits[g[-1]]}\nn={len(g)}" for i, g in zip(x, groups)], fontsize=8)
    ax.set_xlabel(f"Decile of seed hits per read under {selector} (1 = fewest); hit range; eligible reads", color="#52514e")
    ax.set_ylabel("Mismapped reads (% of eligible reads)", color="#52514e")
    ax.set_title(f"Errors concentrate on the reads with the most seed hits\n{reads_set}", loc="left", fontsize=11)
    ax.set_ylim(bottom=0)
    seed_eval.style_axes(ax)
    ax.legend(frameon=False, fontsize=9)
    fig.patch.set_facecolor("#fcfcfb")
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor="#fcfcfb")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reads_dir")
    parser.add_argument("methods", nargs="+")
    parser.add_argument("--selector", required=True)
    parser.add_argument("--out")
    parser.add_argument("--deciles")
    args = parser.parse_args()
    reads_set = Path(args.reads_dir).name
    cmp = {m: compare(args.reads_dir, m) for m in args.methods}
    eligible = [r for r in cmp[args.methods[0]] if all(r in cmp[m] and cmp[m][r][0] for m in args.methods)]
    subsets = {"all": eligible}
    selector_hits = seeds(args.reads_dir, args.selector)
    if selector_hits is not None:
        subsets["selector_top10"] = top_fraction(selector_hits, 0.1)
    method_hits = {m: seeds(args.reads_dir, m) for m in args.methods}
    if all(h is not None for h in method_hits.values()):
        others = args.methods[1:]
        mean_hits = {r: statistics.fmean(method_hits[m][r] for m in others) for r in method_hits[others[0]]}
        subsets["group_top10"] = top_fraction(mean_hits, 0.1)
        first = method_hits[args.methods[0]]
        subsets["hits_differ"] = [r for r in first if any(first[r] != method_hits[m].get(r) for m in others)]
    first = args.methods[0]
    columns = ["reads_set", "subset", "eligible"] + [f"wrong_{m}" for m in args.methods]
    for m in args.methods[1:]:
        columns += [f"only_{first}_correct_vs_{m}", f"only_{m}_correct", f"p_{first}_vs_{m}"]
    rows = []
    eligible_set = set(eligible)
    for name, subset in subsets.items():
        reads = [r for r in subset if r in eligible_set]
        row = {"reads_set": reads_set, "subset": name, "eligible": len(reads)}
        for m in args.methods:
            row[f"wrong_{m}"] = sum(not cmp[m][r][1] for r in reads)
        for m in args.methods[1:]:
            b = sum(cmp[first][r][1] and not cmp[m][r][1] for r in reads)
            c = sum(cmp[m][r][1] and not cmp[first][r][1] for r in reads)
            row.update({f"only_{first}_correct_vs_{m}": b, f"only_{m}_correct": c, f"p_{first}_vs_{m}": f"{mcnemar(b, c):.2g}"})
        rows.append(row)
    lines = ["\t".join(columns)] + ["\t".join(str(r[c]) for c in columns) for r in rows]
    print("\n".join(lines))
    if args.out:
        with open(args.out, "w") as f:
            f.write("\n".join(lines) + "\n")
    if args.deciles:
        if selector_hits is None:
            sys.exit(f"--deciles needs a density run for {args.selector}")
        deciles_plot(args.deciles, reads_set, args.methods, cmp, eligible, selector_hits, args.selector)


if __name__ == "__main__":
    main()
