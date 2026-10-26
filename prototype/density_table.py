#!/usr/bin/env python3
"""Side-by-side density table for minimaxmer prototype runs on the same reads.

Usage: python3 -I density_table.py --baseline DENSITY_SUMMARY.tsv MODE=PER_READ.tsv:SUMMARY.txt ... [--out TABLE.tsv]

DENSITY_SUMMARY.tsv: a seed_eval density run's density_summary.tsv (giraffe's own minimizers on these
reads), the reference density. PER_READ.tsv / SUMMARY.txt: minimaxmer.py --per-read output and the
summary it prints (with the "windows_won  selected_kmers" histogram). Prints numbers only.
"""
import argparse
import csv
import statistics


def quantile(values, q):
    values = sorted(values)
    return values[min(len(values) - 1, int(q * len(values)))]


def windows_won(summary):
    hist, inside = {}, False
    for line in open(summary):
        fields = line.rstrip("\n").split("\t")
        if fields[0] == "windows_won":
            inside = True
        elif inside and len(fields) == 2:
            hist[int(fields[0])] = int(fields[1])
    return hist


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("runs", nargs="+")
    parser.add_argument("--out")
    args = parser.parse_args()
    with open(args.baseline, newline="") as f:
        base = next(csv.DictReader(f, delimiter="\t"))
    base_density = float(base["density_pooled"])
    columns = ["mode", "reads", "total_bp", "selected", "density", "density_vs_" + base["name"],
               "per_read_selected_p10", "per_read_selected_p50", "per_read_selected_p90",
               "per_read_density_p10", "per_read_density_p50", "per_read_density_p90",
               "windows_won_mean", "windows_won_1_fraction", "windows_won_p90", "windows_won_max"]
    rows = [{"mode": base["name"] + " (giraffe)", "reads": base["reads"], "total_bp": base["total_bp"],
             "selected": base["minimizers"], "density": f"{base_density:.6g}", "density_vs_" + base["name"]: "1"}]
    for run in args.runs:
        mode, files = run.split("=", 1)
        per_read_path, summary = files.split(":")
        with open(per_read_path, newline="") as f:
            per_read = [(int(r["length"]), int(r["selected"])) for r in csv.DictReader(f, delimiter="\t")]
        bp, selected = sum(n for n, _ in per_read), sum(s for _, s in per_read)
        won = windows_won(summary)
        won_list = [w for w, n in won.items() for _ in range(n)]
        sel, dens = [s for _, s in per_read], [s / n for n, s in per_read if n]
        rows.append({"mode": mode, "reads": len(per_read), "total_bp": bp, "selected": selected,
                     "density": f"{selected / bp:.6g}", "density_vs_" + base["name"]: f"{selected / bp / base_density:.4g}",
                     "per_read_selected_p10": quantile(sel, 0.1), "per_read_selected_p50": quantile(sel, 0.5),
                     "per_read_selected_p90": quantile(sel, 0.9),
                     "per_read_density_p10": f"{quantile(dens, 0.1):.4g}", "per_read_density_p50": f"{quantile(dens, 0.5):.4g}",
                     "per_read_density_p90": f"{quantile(dens, 0.9):.4g}",
                     "windows_won_mean": f"{statistics.fmean(won_list):.4g}",
                     "windows_won_1_fraction": f"{won.get(1, 0) / len(won_list):.4g}",
                     "windows_won_p90": quantile(won_list, 0.9), "windows_won_max": max(won_list)})
    lines = ["\t".join(columns)] + ["\t".join(str(r.get(c, "")) for c in columns) for r in rows]
    print("\n".join(lines))
    if args.out:
        with open(args.out, "w") as f:
            f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
