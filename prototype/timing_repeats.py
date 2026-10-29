#!/usr/bin/env python3
"""Speed ratios over many repeated timing jobs (seed_eval.py timing --no-exclusive --jobs N).

Usage: python3 -I timing_repeats.py TIMING_DIR JOB_ID ... [--ref mini_k31w20] [--out PREFIX] [--plot PREFIX]

TIMING_DIR: data_dir/results/<reads>/timing_repeats_runs (where timing --no-exclusive writes; exclusive jobs are
in timing/); each JOB_ID has job<ID>/runs.tsv. Only the given jobs are used. Within each job, every condition's load_seconds,
map_seconds and reads_per_cpu_second are the median over its rounds and are divided by the reference condition's
value in the same job. Across jobs: median, quartiles, range of these ratios, and the fraction of jobs with a
ratio above 1 (for seconds: slower than the reference; for reads_per_cpu_second: faster).
The same ratios are also taken against each condition's density-matched reference (group_ref: the
minimizer-labelled condition of its density group in conditions.tsv; the reference group keeps --ref), in the
group_ratio_* columns. --plot PREFIX draws one figure per density group (PREFIX_d2-11.png, ...), every
ratio against that group's minimizer, the giraffe default added for context, same y ranges across figures.
Real reads (runs.tsv with aligned_fraction): also aligned fraction, MAPQ 0/60 fractions and identity over all
runs, the MAPQ histogram pooled over each condition's runs, and placement.tsv summed over jobs.
Only jobs that sacct reports as COMPLETED are used; the others are skipped (named on stderr), and --plot
titles the figure "PRELIMINARY: n of N jobs finished" until all given jobs are.
Writes PREFIX_ratios.tsv, PREFIX_memory.tsv (memory_gb median and range over all runs) (and PREFIX_real.tsv, PREFIX_mapq.tsv, PREFIX_placement.tsv for real reads).
Prints numbers only.
"""
import argparse
import collections
import csv
import re
import statistics
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import seed_eval  # noqa: E402

METRICS = ["load_seconds", "map_seconds", "reads_per_cpu_second"]
REAL_METRICS = ["aligned_fraction", "mapq0_fraction", "mapq60_fraction", "identity_mean"]
MAPQ_BINS = [(0, 0), (1, 9), (10, 29), (30, 59), (60, 60)]


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def spread(values):
    """median, Q1, Q3, min, max (quartiles as statistics.quantiles, inclusive method)."""
    values = sorted(values)
    q1, _, q3 = statistics.quantiles(values, n=4, method="inclusive") if len(values) > 1 else values * 3
    return statistics.median(values), q1, q3, values[0], values[-1]


def write(path, columns, rows):
    lines = ["\t".join(columns)] + ["\t".join(str(r[c]) for c in columns) for r in rows]
    print("\n".join(lines) + "\n")
    if path:
        Path(path).write_text("\n".join(lines) + "\n")


def plot(args, conditions, rows, baseline, ratios_against, n_finished):
    """One figure per density group, PREFIX_d<group>.png (2/11 -> d2-11), one panel per metric. Every ratio in a
    figure, of the group's conditions and of the reference group's (giraffe default, last, for context), is
    ratios_against the group's BASELINE_LABEL condition; each point is one job; a line at the median and at 1.
    A panel has the same y range in every figure of the data set, so the figures can be put side by side."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = seed_eval.method_colors(list(rows.values()))
    reads_set = Path(args.timing_dir).resolve().parent.name
    status = (f" — PRELIMINARY: {n_finished} of {len(args.job_ids)} jobs finished"
              if n_finished < len(args.job_ids) else "")
    file_order = list(rows)
    context = [c for c in conditions if c in rows and rows[c]["density_group"] == seed_eval.REFERENCE_GROUP]
    figures = []   # (group, its reference, its conditions with the density-matched methods first)
    for g, ref in baseline.items():
        if g != seed_eval.REFERENCE_GROUP:
            own = sorted((c for c in conditions if c in rows and rows[c]["density_group"] == g),
                         key=lambda c: (rows[c]["label"] not in seed_eval.DENSITY_MATCHED_LABELS, file_order.index(c)))
            figures.append((g, ref, own))
    if not figures:
        print("no density group with a baseline condition in the runs; nothing plotted", file=sys.stderr)
        return
    data = {(g, c, m): ratios_against(c, ref, m) for g, ref, own in figures for c in own + context for m in METRICS}
    limits = {}
    for m in METRICS:   # union over the figures, with a margin
        values = [x for (_, _, metric), xs in data.items() if metric == m for x in xs] + [1]
        pad = 0.06 * ((max(values) - min(values)) or 0.1)
        limits[m] = (min(values) - pad, max(values) + pad)

    def context_label(c):
        w = re.search(r"-w (\d+)", rows[c]["index_args"])
        return "giraffe default" + (f" (density ≈2/{int(w[1]) + 1})" if w else "") + ",\nfor context"

    for g, ref, own in figures:
        shown = own + context
        fig, axes = plt.subplots(1, len(METRICS), figsize=(4.4 * len(METRICS), 4.9))
        for ax, m in zip(axes, METRICS):
            for i, c in enumerate(shown):
                ratios = data[(g, c, m)]
                color = seed_eval.row_color(rows[c], colors)
                # small fixed horizontal offsets so equal values stay visible; no randomness
                xs = [i + 0.12 * (k - (len(ratios) - 1) / 2) / max(1, len(ratios) / 4) for k in range(len(ratios))]
                ax.scatter(xs, ratios, s=28, color=color, edgecolor="#fcfcfb", lw=0.8, zorder=3)
                if ratios:
                    med = statistics.median(ratios)
                    ax.plot([i - 0.3, i + 0.3], [med, med], color=color, lw=2, zorder=2)
            if context:
                ax.axvline(len(own) - 0.5, color="#c3c2b7", lw=1, zorder=1)
                ax.text(len(own) + (len(context) - 1) / 2, 1.01, context_label(context[0]),
                        transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7, color="#52514e")
            ax.axhline(1, color="#52514e", lw=1, ls="--", zorder=1)
            ax.set_xlim(-0.5, len(shown) - 0.5)
            ax.set_ylim(*limits[m])
            ax.set_xticks(range(len(shown)))
            ax.set_xticklabels([f"{c}\nn={len(data[(g, c, m)])}" for c in shown], rotation=45, ha="right", fontsize=7.5)
            what = "higher = faster" if m == "reads_per_cpu_second" else "higher = slower"
            ax.set_title(f"{m}\n({what}; one point per job)", fontsize=9, loc="left", pad=26)
            seed_eval.style_axes(ax)
            ax.grid(axis="x", visible=False)
        fig.suptitle(f"{reads_set} — density {g} — ratio vs {ref}{status}", x=0.01, ha="left", fontsize=11)
        fig.patch.set_facecolor("#fcfcfb")
        fig.tight_layout()
        path = f"{args.plot}_d{g.replace('/', '-')}.png"
        fig.savefig(path, dpi=150, facecolor="#fcfcfb")
        plt.close(fig)
        print(f"plot {path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("timing_dir")
    parser.add_argument("job_ids", nargs="+")
    parser.add_argument("--ref", default="mini_k31w20")
    parser.add_argument("--out", help="output prefix")
    parser.add_argument("--plot", metavar="PREFIX",
                        help="strip plots of the ratios, one PNG per density group: PREFIX_d<group>.png")
    parser.add_argument("--conditions", default=str(Path(__file__).resolve().parent.parent / "conditions.tsv"),
                        help="conditions.tsv, for the plot colors as in collect")
    args = parser.parse_args()
    out = lambda suffix: f"{args.out}_{suffix}.tsv" if args.out else None

    # Only jobs Slurm reports as COMPLETED; running, pending or failed ones are skipped
    states = subprocess.run(["sacct", "-X", "-n", "-P", "-o", "JobID,State", "-j", ",".join(args.job_ids)],
                            check=True, capture_output=True, text=True).stdout
    states = dict(line.split("|")[:2] for line in states.splitlines() if line)
    finished = [job for job in args.job_ids if states.get(job) == "COMPLETED"]
    for job in args.job_ids:
        if job not in finished:
            print(f"skipping job {job}: {states.get(job, 'unknown to sacct')}", file=sys.stderr)
    if not finished:   # not an error: nothing to summarize yet, and no empty table or plot is written
        print(f"no finished jobs yet (0 of {len(args.job_ids)}); nothing written", file=sys.stderr)
        sys.exit(0)

    runs = {}
    for job in finished:
        path = Path(args.timing_dir) / f"job{job}" / "runs.tsv"
        if not path.is_file():
            sys.exit(f"no {path}")
        runs[job] = read_tsv(path)
    conditions = list(dict.fromkeys(r["condition"] for rows in runs.values() for r in rows))
    if args.ref not in conditions:
        sys.exit(f"reference {args.ref} not in the runs")

    # Per job: condition -> metric -> median over rounds
    medians = {}
    for job, rows in runs.items():
        by = collections.defaultdict(lambda: collections.defaultdict(list))
        for r in rows:
            for m in METRICS:
                if r[m] not in ("", "NA"):
                    by[r["condition"]][m].append(float(r[m]))
        medians[job] = {c: {m: statistics.median(v) for m, v in ms.items()} for c, ms in by.items()}
    print("jobs:", len(runs), "nodes:", " ".join(sorted({r["node"] for rows in runs.values() for r in rows})),
          "runs per condition:", {c: sum(r["condition"] == c for rows in runs.values() for r in rows) for c in conditions})
    print()

    # Second reference: the BASELINE_LABEL condition of the same density group (as placement.tsv pairs them);
    # the reference group and conditions without such a baseline keep --ref
    cond_rows = {c["name"]: c for c in seed_eval.load_conditions(args.conditions)}
    baseline = {r["density_group"]: name for name, r in cond_rows.items()
                if r["label"] == seed_eval.BASELINE_LABEL and name in conditions}
    group_ref = {c: baseline.get(cond_rows[c]["density_group"], args.ref)
                 if c in cond_rows and cond_rows[c]["density_group"] != seed_eval.REFERENCE_GROUP else args.ref
                 for c in conditions}

    def ratios_against(c, ref, m):
        """Per job: c's median over rounds divided by ref's median over rounds in the same job."""
        return [medians[j][c][m] / medians[j][ref][m] for j in runs
                if m in medians[j].get(c, {}) and m in medians[j].get(ref, {})]

    def described(prefix, ratios):
        med, q1, q3, lo, hi = spread(ratios)
        return {f"{prefix}_median": f"{med:.4f}", f"{prefix}_q1": f"{q1:.4f}", f"{prefix}_q3": f"{q3:.4f}",
                f"{prefix}_min": f"{lo:.4f}", f"{prefix}_max": f"{hi:.4f}",
                f"jobs_{prefix}_above_1": f"{sum(x > 1 for x in ratios)}/{len(ratios)}"}

    stats = ["median", "q1", "q3", "min", "max"]
    columns = (["condition", "metric", "jobs"] + [f"ratio_{s}" for s in stats] + ["jobs_ratio_above_1", "value_median",
               "group_ref"] + [f"group_ratio_{s}" for s in stats] + ["jobs_group_ratio_above_1"])
    rows, group_ratios = [], {}
    for m in METRICS:
        for c in conditions:
            ratios = ratios_against(c, args.ref, m)
            group_ratios[(c, m)] = ratios_against(c, group_ref[c], m)
            values = [medians[j][c][m] for j in runs if m in medians[j].get(c, {})]
            if not ratios:
                continue
            rows.append({"condition": c, "metric": m, "jobs": len(ratios), **described("ratio", ratios),
                         "value_median": f"{statistics.median(values):.4g}", "group_ref": group_ref[c],
                         **described("group_ratio", group_ratios[(c, m)])})
    write(out("ratios"), columns, rows)

    # Giraffe's reported memory footprint over every run (not a ratio: it does not depend on the node's load)
    rows = []
    for c in conditions:
        values = [float(r["memory_gb"]) for rows_ in runs.values() for r in rows_
                  if r["condition"] == c and r["memory_gb"] not in ("", "NA")]
        if values:
            rows.append({"condition": c, "runs": len(values), "memory_gb_median": f"{statistics.median(values):.2f}",
                         "memory_gb_min": f"{min(values):.2f}", "memory_gb_max": f"{max(values):.2f}"})
    write(out("memory"), ["condition", "runs", "memory_gb_median", "memory_gb_min", "memory_gb_max"], rows)
    if args.plot:
        plot(args, conditions, cond_rows, baseline, ratios_against, len(finished))

    if "aligned_fraction" not in next(iter(runs.values()))[0]:
        return
    # Real reads: statistics over every run, MAPQ histograms pooled per condition, placement summed over jobs
    columns = ["condition", "runs"] + [f"{m}_{s}" for m in REAL_METRICS for s in ("median", "min", "max")]
    rows = []
    for c in conditions:
        row = {"condition": c}
        for m in REAL_METRICS:
            values = [float(r[m]) for rows_ in runs.values() for r in rows_ if r["condition"] == c and r[m] not in ("", "NA")]
            row["runs"] = len(values)
            row[f"{m}_median"], row[f"{m}_min"], row[f"{m}_max"] = (f"{statistics.median(values):.6g}",
                                                                    f"{min(values):.6g}", f"{max(values):.6g}")
        rows.append(row)
    write(out("real"), columns, rows)

    columns = ["condition", "runs", "aligned"] + [f"mapq_{a}" if a == b else f"mapq_{a}_{b}" for a, b in MAPQ_BINS]
    rows = []
    for c in conditions:
        hist, n = collections.Counter(), 0
        for job, rows_ in runs.items():
            for r in rows_:
                if r["condition"] == c:
                    n += 1
                    for h in read_tsv(Path(args.timing_dir) / f"job{job}" / f"round{r['round']}_{c}" / "mapq_hist.tsv"):
                        hist[int(h["mapq"])] += int(h["count"])
        total = sum(hist.values())
        row = {"condition": c, "runs": n, "aligned": total}
        for a, b in MAPQ_BINS:
            key = f"mapq_{a}" if a == b else f"mapq_{a}_{b}"
            row[key] = f"{sum(v for q, v in hist.items() if a <= q <= b) / total:.5f}"
        rows.append(row)
    write(out("mapq"), columns, rows)

    numeric = ["reads", "aligned_both", "aligned_only_condition", "aligned_only_against", "different_start",
               "higher_score_condition", "higher_score_against", "equal_score"]
    sums = collections.defaultdict(lambda: collections.Counter())
    jobs = collections.Counter()
    for job in runs:
        path = Path(args.timing_dir) / f"job{job}" / "placement.tsv"
        if not path.is_file():
            print(f"no {path}", file=sys.stderr)
            continue
        for r in read_tsv(path):
            key = (r["condition"], r["against"])
            jobs[key] += 1
            sums[key].update({k: int(r[k]) for k in numeric})
    columns = ["condition", "against", "jobs"] + numeric
    write(out("placement"), columns, [{"condition": a, "against": b, "jobs": jobs[(a, b)], **sums[(a, b)]}
                                      for a, b in sums])


if __name__ == "__main__":
    main()
