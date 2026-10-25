#!/usr/bin/env python3
"""reads_count distribution of the k-mers that occur once in the graph, from a minimax_counts counts file.

Usage:
  python3 -I plot_reads_count.py hist PREFIX.mmc OUT.tsv
      reads_count histogram of the k-mers with genome_count 1 (reads_count, kmers); reads the counts file
      in chunks, numbers only
  python3 -I plot_reads_count.py plot HIST.tsv OUT.png --coverage 66 [--xmax 150]
      k-mers per reads_count, with the coverage and half of it marked
"""
import argparse
import csv

MAGIC = b"MMXCNT01"
CHUNK = 1 << 27  # entries per read: 512 MiB of values


def hist(mmc, out):
    import numpy as np
    with open(mmc, "rb") as f:
        header = f.read(24)
    if header[:8] != MAGIC:
        raise SystemExit(f"{mmc}: not a minimax_counts file")
    n = int.from_bytes(header[16:24], "little")
    values = np.memmap(mmc, dtype="<u4", mode="r", offset=24 + 8 * n, shape=(n,))
    counts = np.zeros(65536, dtype=np.int64)
    for start in range(0, n, CHUNK):
        chunk = np.asarray(values[start:start + CHUNK])
        counts += np.bincount(chunk[(chunk >> 16) == 1] & 0xFFFF, minlength=65536)
    with open(out, "w") as f:
        f.write("reads_count\tkmers\n")
        f.writelines(f"{r}\t{c}\n" for r, c in enumerate(counts) if c)
    print(f"entries\t{n}\ngenome_count_1\t{counts.sum()}\npeak_from_5\t{int(np.argmax(counts[5:65535])) + 5}")


def plot(tsv, out, coverage, xmax):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with open(tsv, newline="") as f:
        rows = {int(r["reads_count"]): int(r["kmers"]) for r in csv.DictReader(f, delimiter="\t")}
    total = sum(rows.values())
    xs = list(range(0, xmax + 1))
    ys = [rows.get(x, 0) / 1e6 for x in xs]
    beyond = sum(c for r, c in rows.items() if r > xmax)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    fig.patch.set_facecolor("#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    ax.grid(axis="y", color="#e6e5e0", linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#c3c2b7")
    ax.tick_params(colors="#52514e")

    ax.fill_between(xs, ys, step="mid", color="#2a78d6", alpha=0.18, linewidth=0)
    ax.step(xs, ys, where="mid", color="#2a78d6", linewidth=2)
    top = max(ys) * 1.15
    for x, label in ((coverage / 2, f"{coverage / 2:g}\n(coverage / 2)"), (coverage, f"{coverage:g}\n(coverage)")):
        ax.axvline(x, color="#8a8984", linewidth=1.2, linestyle="--")
        ax.text(x + 1.5, top * 0.97, label, color="#52514e", fontsize=9, va="top")
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, top)
    ax.set_xlabel("reads_count (occurrences of the 31-mer in the reads)", color="#52514e")
    ax.set_ylabel("31-mers (millions)", color="#52514e")
    ax.set_title("Read counts of 31-mers that occur once in the graph", loc="left", color="#2b2a27", fontsize=12)
    peaks = [max(range(lo, hi), key=lambda r: rows.get(r, 0)) for lo, hi in ((5, int(coverage * 0.75)),
                                                                             (int(coverage * 0.75), xmax))]
    note = (f"{total / 1e9:.2f} billion 31-mers with genome_count 1; peaks at {peaks[0]} and {peaks[1]}.\n"
            f"{rows.get(0, 0) / 1e6:.1f} million have reads_count 0 (absent from the reads or seen once); "
            f"{beyond / 1e6:.1f} million beyond {xmax} not shown.")
    fig.text(0.01, 0.01, note, color="#8a8984", fontsize=8)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(out, dpi=150, facecolor="#fcfcfb")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="mode", required=True)
    h = sub.add_parser("hist")
    h.add_argument("mmc")
    h.add_argument("out")
    p = sub.add_parser("plot")
    p.add_argument("tsv")
    p.add_argument("out")
    p.add_argument("--coverage", type=float, required=True)
    p.add_argument("--xmax", type=int, default=150)
    args = parser.parse_args()
    if args.mode == "hist":
        hist(args.mmc, args.out)
    else:
        plot(args.tsv, args.out, args.coverage, args.xmax)


if __name__ == "__main__":
    main()
