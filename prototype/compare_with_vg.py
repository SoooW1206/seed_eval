#!/usr/bin/env python3
"""Compare the prototype's per-read selection counts with Giraffe's.

Usage: python3 -I compare_with_vg.py CONSTANT.tsv VG_PER_READ.tsv [COUNTED.tsv] [--out comparison.tsv]

CONSTANT.tsv / COUNTED.tsv: the prototype's --per-read output (read, length, selected) with
constant counts and with reads counts. VG_PER_READ.tsv: per_read.tsv of a seed_eval density run
(read, length, minimizers, seeds), where `minimizers` is Giraffe's stage.minimizer.results.
Reads are matched by name, never by line order; names and lengths must agree on both sides.
Prints and writes numbers only.
"""
import argparse
import collections
import csv
import sys


def load(path, count_column):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    table = {r["read"]: (int(r["length"]), int(r[count_column])) for r in rows}
    if len(table) != len(rows):
        sys.exit(f"{path}: duplicate read names")
    return table


def compare(a, b, name_a, name_b):
    """Rows (section, metric, value) comparing per-read counts a - b after checking names and lengths."""
    if set(a) != set(b):
        sys.exit(f"read names differ: {len(set(a) - set(b))} only in {name_a}, {len(set(b) - set(a))} only in {name_b}")
    bad_length = sum(1 for read in a if a[read][0] != b[read][0])
    if bad_length:
        sys.exit(f"{bad_length} reads have different lengths in {name_a} and {name_b}")
    section = f"{name_a}_vs_{name_b}"
    diffs = collections.Counter(a[read][1] - b[read][1] for read in a)
    total_bp = sum(length for length, _ in a.values())
    rows = [(section, "reads", len(a)), (section, "lengths_equal", len(a)),
            (section, "reads_equal", diffs.get(0, 0)), (section, "reads_different", len(a) - diffs.get(0, 0)),
            (section, f"total_{name_a}", sum(c for _, c in a.values())),
            (section, f"total_{name_b}", sum(c for _, c in b.values())),
            (section, f"density_{name_a}", round(sum(c for _, c in a.values()) / total_bp, 6)),
            (section, f"density_{name_b}", round(sum(c for _, c in b.values()) / total_bp, 6))]
    rows += [(section, f"difference_{d:+d}", n) for d, n in sorted(diffs.items()) if d != 0]
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("constant")
    parser.add_argument("vg_per_read")
    parser.add_argument("counted", nargs="?")
    parser.add_argument("--out")
    args = parser.parse_args()
    constant = load(args.constant, "selected")
    rows = compare(constant, load(args.vg_per_read, "minimizers"), "prototype_constant", "vg")
    if args.counted:
        rows += compare(load(args.counted, "selected"), constant, "prototype_counted", "prototype_constant")
    for section, metric, value in rows:
        print(f"{section}\t{metric}\t{value}")
    if args.out:
        with open(args.out, "w") as f:
            f.write("section\tmetric\tvalue\n")
            f.writelines(f"{s}\t{m}\t{v}\n" for s, m, v in rows)


if __name__ == "__main__":
    main()
