#!/usr/bin/env python3
"""MAPQ distribution and precision-recall by MAPQ threshold, as in the lab's plot-pr.R.

Usage: python3 -I mapq_pr.py READS_DIR METHOD ... [--out TABLE.tsv] [--plot PR.png]

READS_DIR: data_dir/results/<reads>. Uses each method's compare.tsv (rep1/ when present), eligible reads only.
Table: reads and wrong reads per MAPQ bin (0, 1-9, 10-29, 30-59, 60) and with MAPQ >= 60, 30, 1, plus a paired
exact test of the first method against each other one on "wrong at MAPQ 60" (reads only the first / only the
other gets wrong at MAPQ 60). Plot: as long-read-giraffe-experiments' plot-pr.R: going down from MAPQ 60,
TP and FP are the correct and wrong reads at or above the threshold, Precision = TP / (TP + FP),
Recall = TP / eligible reads; axes -log10(1 - Recall) and -log10(1 - Precision); point size: reads at that
MAPQ; MAPQ 0 and 60 labelled. Prints numbers only.
"""
import argparse
import collections
import csv
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import seed_eval  # noqa: E402
from hard_reads import mcnemar  # noqa: E402

BINS = [(0, 0), (1, 9), (10, 29), (30, 59), (60, 60)]
THRESHOLDS = [60, 30, 1]
LABELS = {"mini": "minimizer", "sync": "closed syncmer", "miniW": "weighted minimizer", "minimax": "minimaxmer"}


def mapq(reads_dir, c):
    """read -> (MAPQ, correct) for eligible reads; compare.tsv's mq column."""
    for p in (Path(reads_dir) / c / "rep1" / "compare.tsv", Path(reads_dir) / c / "compare.tsv"):
        if p.is_file():
            with open(p, newline="") as f:
                return {r["read"]: (int(r["mq"]), r["correct"] == "1") for r in csv.DictReader(f, delimiter="\t")
                        if r["eligible"] == "1"}
    sys.exit(f"no compare.tsv for {c}")


def color(m):
    colors = dict(zip(["minimizer", "closed syncmer", "weighted minimizer", "minimaxmer"], seed_eval.METHOD_COLORS))
    return colors.get(LABELS.get(m.split("_")[0]), seed_eval.REFERENCE_COLOR)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reads_dir")
    parser.add_argument("methods", nargs="+")
    parser.add_argument("--out")
    parser.add_argument("--plot")
    parser.add_argument("--reference", help="a condition drawn in gray as a reference, not as a compared method")
    args = parser.parse_args()
    reads_set = Path(args.reads_dir).name
    data = {m: mapq(args.reads_dir, m) for m in args.methods}
    first = args.methods[0]
    columns = ["reads_set", "condition", "eligible"]
    columns += [f"mapq_{a}" if a == b else f"mapq_{a}_{b}" for a, b in BINS]
    columns += [f"wrong_mapq_{a}" if a == b else f"wrong_mapq_{a}_{b}" for a, b in BINS]
    columns += [f"kept_mapq_ge{t}" for t in THRESHOLDS] + [f"wrong_mapq_ge{t}" for t in THRESHOLDS]
    columns += ["wrong_at_60_only_" + first, "wrong_at_60_only_this", "wrong_at_60_both", "p_wrong_at_60"]
    rows = []
    for m in args.methods:
        reads = data[m]
        row = {"reads_set": reads_set, "condition": m, "eligible": len(reads)}
        for a, b in BINS:
            sub = [ok for q, ok in reads.values() if a <= q <= b]
            key = f"{a}" if a == b else f"{a}_{b}"
            row[f"mapq_{key}"], row[f"wrong_mapq_{key}"] = len(sub), sum(not ok for ok in sub)
        for t in THRESHOLDS:
            row[f"kept_mapq_ge{t}"] = sum(q >= t for q, _ in reads.values())
            row[f"wrong_mapq_ge{t}"] = sum(q >= t and not ok for q, ok in reads.values())
        if m != first:
            bad = lambda d, r: r in d and d[r][0] >= 60 and not d[r][1]
            shared = [r for r in data[first] if r in reads]
            only_first = sum(bad(data[first], r) and not bad(reads, r) for r in shared)
            only_this = sum(bad(reads, r) and not bad(data[first], r) for r in shared)
            row.update({"wrong_at_60_only_" + first: only_first, "wrong_at_60_only_this": only_this,
                        "wrong_at_60_both": sum(bad(data[first], r) and bad(reads, r) for r in shared),
                        "p_wrong_at_60": f"{mcnemar(only_first, only_this):.2g}"})
        rows.append(row)
    lines = ["\t".join(columns)] + ["\t".join(str(r.get(c, "")) for c in columns) for r in rows]
    print("\n".join(lines))
    if args.out:
        with open(args.out, "w") as f:
            f.write("\n".join(lines) + "\n")
    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 5.5))
        # Axis range as plot-pr.R: 1e-4, extended to 1e-6 above 10,000 reads per condition and to 1e-9 above
        # 1,000,000. As there (ggplot scale limits), points beyond the range are not drawn.
        reads_per_condition = statistics.fmean(len(data[m]) for m in args.methods)
        limit = 9 if reads_per_condition > 1000000 else 6 if reads_per_condition > 10000 else 4
        labels = []
        for index, m in enumerate(args.methods):
            line_color = seed_eval.REFERENCE_COLOR if m == args.reference else color(m)
            reads = data[m]
            by = collections.defaultdict(lambda: [0, 0])
            for q, ok in reads.values():
                by[q][0 if ok else 1] += 1
            tp = fp = 0
            xs, ys, sizes, qs = [], [], [], []
            for q in sorted(by, reverse=True):
                tp, fp = tp + by[q][0], fp + by[q][1]
                precision, recall = tp / (tp + fp), tp / len(reads)
                if precision < 1 and recall < 1:   # -log10(0) is undefined; plot-pr.R drops these too
                    x, y = -math.log10(1 - recall), -math.log10(1 - precision)
                    if x <= limit and y <= limit:
                        xs.append(x)
                        ys.append(y)
                        sizes.append(sum(by[q]))
                        qs.append(q)
            ax.plot(xs, ys, color=line_color, lw=2, label=m + (" (reference)" if m == args.reference else ""))
            ax.scatter(xs, ys, s=[10 + 60 * s / len(reads) for s in sizes], color=line_color, edgecolor="#fcfcfb", lw=1, zorder=3)
            labels += [(q, x, y, index, line_color) for x, y, q in zip(xs, ys, qs) if q in (0, 60)]
        # The MAPQ 0 and 60 points of the methods sit close together. As geom_text_repel does in plot-pr.R,
        # their labels are moved apart: one column per MAPQ (left of the 60 points, right of the 0 points),
        # stacked at a fixed spacing from the highest point, down when that point is near the top edge and up
        # otherwise, each tied to its point by a line in the method's color.
        step = 0.035 * limit
        for q in (0, 60):
            group = sorted((label for label in labels if label[0] == q), key=lambda label: label[3])
            if not group:
                continue
            top = max(y for _, _, y, _, _ in group)
            down = top > limit - 0.6
            column = (min(x for _, x, _, _, _ in group) - 0.06 * limit if q == 60
                      else max(x for _, x, _, _, _ in group) + 0.06 * limit)
            for rank, (_, x, y, _, line_color) in enumerate(group):
                text_y = top - step * (rank + 1) if down else top + step * (rank + 1)
                ax.annotate(str(q), (x, y), xytext=(column, text_y), textcoords="data", fontsize=8, color="#52514e",
                            ha="center", va="center",
                            arrowprops=dict(arrowstyle="-", color=line_color, lw=0.8, shrinkA=4, shrinkB=3))
        ticks = list(range(0, limit + 1))
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"1e-{t}" for t in ticks])
        ax.set_yticks(ticks)
        ax.set_yticklabels([f"1e-{t}" for t in ticks])
        ax.set_xlim(0, limit)
        ax.set_ylim(0, limit)
        ax.set_xlabel("1 - Recall", color="#52514e")
        ax.set_ylabel("1 - Precision", color="#52514e")
        ax.set_title(f"Precision-recall by MAPQ threshold, {reads_set}\n"
                     "(as plot-pr.R in long-read-giraffe-experiments; point size: reads at that MAPQ)",
                     loc="left", fontsize=10)
        seed_eval.style_axes(ax)
        ax.legend(frameon=False, fontsize=9)
        fig.patch.set_facecolor("#fcfcfb")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150, facecolor="#fcfcfb")
        plt.close(fig)


if __name__ == "__main__":
    main()
