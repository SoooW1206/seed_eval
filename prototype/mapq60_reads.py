#!/usr/bin/env python3
"""Read-by-read table of the reads wrong at MAPQ 60, with the region their truth position falls in.

Usage: python3 -I mapq60_reads.py MAPQ60_DIR ANNOTATION_DIR METHOD ... [--out TABLE.tsv] [--summary SUMMARY.tsv]

MAPQ60_DIR: results/<reads>/mapq60 from jobs/mapq60_provenance_*.sbatch: truth_positions.tsv and, per method,
positions.tsv (annotated: name, MAPQ, score, refpos list), prov_positions.tsv (stage counts) and compare.tsv.
ANNOTATION_DIR: chm13v2.0_censat_v2.1.bed and chm13v2.0_SD.bed (T2T, chromosome names without "CHM13#0#").
A refpos list has one entry per node; each read's span is min..max offset on the path with the most entries.
distance: bp between the mapped and the truth span (0 when they overlap), "other_chrom", or "unmapped".
Regions: censat classes (name up to the first "_" or "(") overlapping the truth span, with bp; SD: bp of the
truth span inside segmental duplications. sd_partner: for a wrong read, the mapped span lies in an SD whose
other copy (the SD name) overlaps the truth span, i.e. the read went to the duplicated copy.
caveat: truth_on_chrY (the lab ignores chrY truth) and mapped_without_CHM13_position (some method's alignment
has no CHM13 refpos after annotate, so gamcompare calls it wrong though it may not be).
Prints numbers and names only.
"""
import argparse
import collections
import csv
import re
import sys
from pathlib import Path

PREFIX = "CHM13#0#"


def span(refpos):
    """refpos list 'path:offset:strand,...' -> (chrom, start, end, strand, other paths) of the main path."""
    by = collections.defaultdict(list)
    for item in filter(None, refpos.split(",")):
        path, offset, strand = item.rsplit(":", 2)
        by[path].append((int(offset), strand))
    if not by:
        return None
    path = max(by, key=lambda p: len(by[p]))
    offsets = [o for o, _ in by[path]]
    strand = collections.Counter(s for _, s in by[path]).most_common(1)[0][0]
    return path.removeprefix(PREFIX), min(offsets), max(offsets) + 1, strand, len(by) - 1


def distance(a, b):
    if a is None or b is None:
        return "unmapped"
    if a[0] != b[0]:
        return "other_chrom"
    return max(0, max(a[1], b[1]) - min(a[2], b[2]))


def load_bed(path):
    bed = collections.defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            bed[fields[0]].append((int(fields[1]), int(fields[2]), fields[3], fields[8] if len(fields) > 8 else ""))
    return bed


def overlaps(bed, s):
    if s is None:
        return []
    return [(a, b, name, rgb, min(b, s[2]) - max(a, s[1])) for a, b, name, rgb in bed.get(s[0], [])
            if a < s[2] and s[1] < b]


def covered(intervals, s):
    """bp of span s covered by the union of intervals."""
    points = sorted((max(a, s[1]), min(b, s[2])) for a, b, *_ in intervals)
    total, end = 0, s[1]
    for a, b in points:
        a = max(a, end)
        if b > a:
            total, end = total + b - a, b
    return total


def stages(text):
    return dict((k.strip(), int(v)) for k, v in (x.split("=") for x in text.split(";") if x))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mapq60_dir")
    parser.add_argument("annotation_dir")
    parser.add_argument("methods", nargs="+")
    parser.add_argument("--out")
    parser.add_argument("--summary")
    args = parser.parse_args()
    d, a = Path(args.mapq60_dir), Path(args.annotation_dir)
    censat, sd = load_bed(a / "chm13v2.0_censat_v2.1.bed"), load_bed(a / "chm13v2.0_SD.bed")
    with open(d / "truth_positions.tsv") as f:
        truth = {r[0]: span(r[1]) for r in csv.reader(f, delimiter="\t")}
    mapped, verdict = {}, {}
    for m in args.methods:
        # refpos only exists after annotate (positions.tsv); stage counts come from the provenance GAM
        with open(d / m / "prov_positions.tsv") as f:
            counts = {r[0]: stages(r[4]) for r in csv.reader(f, delimiter="\t")}
        with open(d / m / "positions.tsv") as f:
            mapped[m] = {r[0]: (int(r[1]), span(r[3]), counts[r[0]]) for r in csv.reader(f, delimiter="\t")}
        with open(d / m / "compare.tsv", newline="") as f:
            verdict[m] = {r["read"]: r["correct"] == "1" for r in csv.DictReader(f, delimiter="\t")}

    columns = ["read", "caveat", "truth", "truth_len", "truth_other_paths", "censat", "sd_bp", "sd_fraction",
               "wrong_at_60", "wrong_below_60", "correct"]
    for m in args.methods:
        columns += [f"{m}_{x}" for x in ("correct", "mapq", "mapped", "distance", "sd_partner", "seeds", "chains")]
    rows = []
    for read in sorted(truth, key=lambda r: (truth[r][0], truth[r][1]) if truth[r] else ("", 0)):
        t = truth[read]
        classes = collections.Counter()
        for _, _, name, _, bp in overlaps(censat, t):
            classes[re.split(r"[_(]", name)[0]] += bp
        sd_hits = overlaps(sd, t)
        sd_bp = covered(sd_hits, t) if t else 0
        row = {"read": read, "truth": f"{t[0]}:{t[1]}-{t[2]}:{t[3]}" if t else "NA", "truth_len": t[2] - t[1] if t else "NA",
               "truth_other_paths": t[4] if t else "NA",
               "censat": ",".join(f"{k}:{v}" for k, v in classes.most_common()) or "none",
               "sd_bp": sd_bp, "sd_fraction": f"{sd_bp / (t[2] - t[1]):.2f}" if t else "NA"}
        groups = collections.defaultdict(list)
        for m in args.methods:
            q, s, st = mapped[m][read]
            ok = verdict[m][read]
            partner = "NA"
            if not ok and s is not None and t is not None:
                partner = "no"
                for _, _, name, _, _ in overlaps(sd, s):
                    hit = re.fullmatch(r"(\w+):(\d+)-(\d+)", name)
                    if hit and hit[1] == t[0] and int(hit[2]) < t[2] and t[1] < int(hit[3]):
                        partner = "yes"
            groups["correct" if ok else "wrong_at_60" if q >= 60 else "wrong_below_60"].append(m)
            row.update({f"{m}_correct": int(ok), f"{m}_mapq": q,
                        f"{m}_mapped": f"{s[0]}:{s[1]}-{s[2]}:{s[3]}" if s else "unmapped",
                        f"{m}_distance": distance(s, t), f"{m}_sd_partner": partner,
                        f"{m}_seeds": st.get("seed", "NA"), f"{m}_chains": st.get("chain", "NA")})
        for g in ("wrong_at_60", "wrong_below_60", "correct"):
            row[g] = ",".join(groups[g]) or "-"
        caveats = []
        if t and t[0] == "chrY":
            caveats.append("truth_on_chrY")   # the lab's gamcompare ignores CHM13#0#chrY truth positions
        if any(mapped[m][read][1] is None for m in args.methods):
            caveats.append("mapped_without_CHM13_position")   # gamcompare calls it wrong; unclear if it is
        row["caveat"] = ",".join(caveats) or "-"
        rows.append(row)
    lines = ["\t".join(columns)] + ["\t".join(str(r[c]) for c in columns) for r in rows]
    if args.out:
        Path(args.out).write_text("\n".join(lines) + "\n")

    # Summary: per read, how the methods split, and the region
    first = args.methods[0]
    summary = ["\t".join(["read", "caveat", "truth", "censat", "sd_fraction", f"{first}_mapq", f"{first}_distance",
                          f"{first}_sd_partner", "wrong_at_60", "wrong_below_60", "correct"])]
    for r in rows:
        summary.append("\t".join(str(r[c]) for c in ["read", "caveat", "truth", "censat", "sd_fraction", f"{first}_mapq",
                                                      f"{first}_distance", f"{first}_sd_partner", "wrong_at_60",
                                                      "wrong_below_60", "correct"]))
    print("\n".join(summary))
    if args.summary:
        Path(args.summary).write_text("\n".join(summary) + "\n")


if __name__ == "__main__":
    main()
