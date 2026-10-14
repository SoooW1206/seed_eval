#!/usr/bin/env python3
"""Seed-method evaluator for vg giraffe.

Reads a conditions table (one seed-index configuration per row) and generates
one sbatch job per condition for the chosen stage. Scripts are only written and
printed by default; pass --submit to actually submit them.

The driver is data-driven: every condition-specific piece of a command comes
from its row in the table, so adding a condition means adding a row.

Columns:
  name           unique condition name; used for output file names
  vg_binary      absolute path to the vg executable used for every step
  index_args     arguments passed verbatim to `vg minimizer`
  index_mem      sbatch --mem for the index job (e.g. 80G). Measured at k=31 on
                 the HG002 sampled graph: plain ~56 GiB -> 80G, -W ~150 GiB -> 256G
  giraffe_args   extra arguments passed verbatim to `vg giraffe` (may be empty)
  density_group  only used for result tables and plots
  label          only used for result tables and plots

Stages:
  index    vg minimizer -> indexes/<name>.min + .zipcodes
  map      giraffe -> annotate -> gamcompare -> stats, on --reads
  density  giraffe --track-provenance on --reads; per-read seed counts
  filters  one small job: per-filter pass/fail totals and kept minimizers from the density run
  collect  (runs locally, no sbatch) parse logs of all stages into
           results/<reads>/results.tsv and results/<reads>/comparison.png
"""
import argparse
import csv
import math
import re
import subprocess
import sys
from pathlib import Path

DATA_DIR = Path("/private/groups/patenlab/twang285/seed-eval")
GRAPH_NAME = "hprc-v2.0-mc-chm13-eval-sampled16o_fragmentlinked-for-real-hifi-HG002-full"
GBZ = DATA_DIR / f"{GRAPH_NAME}.gbz"
DIST = DATA_DIR / f"{GRAPH_NAME}.dist"
READS_DIR = Path("/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/reads/sim/hifi/HG002")
DEFAULT_READS = READS_DIR / "HG002-sim-hifi-1k.gam"
PARTITION = "medium"
THREADS = 16  # vg minimizer: more than 16 threads rarely helps
GIRAFFE_PRESET = "hifi"
COMPARE_RANGE = 200  # bp window for gamcompare to call a read correct
MAP_MEM = "64G"  # measured on 1k reads: giraffe footprint ~41 GB, job peak ~36 GiB
MAP_TIME = "02:00:00"
INDEX_TIME = "12:00:00"

REQUIRED_COLUMNS = ["name", "vg_binary", "index_args", "index_mem"]
OPTIONAL_COLUMNS = ["giraffe_args", "density_group", "label"]


def load_conditions(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = set(REQUIRED_COLUMNS + OPTIONAL_COLUMNS) - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"{path}: missing columns {sorted(missing)}")
        conditions = [{k: (v or "").strip() for k, v in row.items()} for row in reader]
    validate(conditions, path)
    return conditions


def validate(conditions, path):
    for i, c in enumerate(conditions, start=2):  # line 1 is the header
        empty = [col for col in REQUIRED_COLUMNS if not c[col]]
        if empty:
            sys.exit(f"{path} line {i}: empty required columns {empty}")
    names = [c["name"] for c in conditions]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        sys.exit(f"{path}: duplicate condition names {duplicates}")


def sbatch_header(job_name, mem, time, log, partition=PARTITION, cpus=THREADS):
    return f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH -c {cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --output={log}
set -euo pipefail
"""


def index_prefix(c):
    return DATA_DIR / "indexes" / c["name"]


def require_index(c):
    out = index_prefix(c)
    return f"""for f in {out}.min {out}.zipcodes; do
    [ -s "$f" ] || {{ echo "missing index file $f; run the index stage first" >&2; exit 1; }}
done
"""


def giraffe_cmd(c, reads, extra=""):
    out = index_prefix(c)
    return (f"{c['vg_binary']} giraffe -p -b {GIRAFFE_PRESET} -t {THREADS} {c['giraffe_args']} {extra}"
            f" \\\n    -Z {GBZ} -d {DIST} -m {out}.min -z {out}.zipcodes -G {reads}")


def index_script(c, reads):
    out = index_prefix(c)
    log = DATA_DIR / "logs" / f"index_{c['name']}.%j.log"
    return sbatch_header(f"idx-{c['name']}", c["index_mem"], INDEX_TIME, log) + f"""mkdir -p {out.parent}
/usr/bin/time -v {c['vg_binary']} minimizer {c['index_args']} -p -t {THREADS} \\
    -d {DIST} -o {out}.min.tmp -z {out}.zipcodes.tmp {GBZ}
mv {out}.min.tmp {out}.min
mv {out}.zipcodes.tmp {out}.zipcodes
echo "min size (bytes): $(stat -c %s {out}.min)"
echo "zipcodes size (bytes): $(stat -c %s {out}.zipcodes)"
"""


def results_dir(c, reads):
    return DATA_DIR / "results" / Path(reads).stem / c["name"]


def map_script(c, reads):
    o = results_dir(c, reads)
    log = DATA_DIR / "logs" / f"map_{Path(reads).stem}_{c['name']}.%j.log"
    vg = c["vg_binary"]
    return sbatch_header(f"map-{c['name']}", MAP_MEM, MAP_TIME, log) + require_index(c) + f"""O={o}
mkdir -p $O
/usr/bin/time -v -o $O/giraffe.time {giraffe_cmd(c, reads)} \\
    --report-name $O/report.tsv > $O/mapped.gam.tmp 2> $O/giraffe.log
mv $O/mapped.gam.tmp $O/mapped.gam
/usr/bin/time -v -o $O/annotate.time {vg} annotate -t {THREADS} -a $O/mapped.gam -x {GBZ} -m --search-limit=-1 \\
    > $O/annot.gam.tmp
mv $O/annot.gam.tmp $O/annot.gam
/usr/bin/time -v -o $O/gamcompare.time {vg} gamcompare --threads {THREADS} --range {COMPARE_RANGE} \\
    $O/annot.gam {reads} -T -a {c['name']} > $O/compare.tsv.tmp 2> $O/gamcompare.txt
mv $O/compare.tsv.tmp $O/compare.tsv
{vg} stats -a $O/mapped.gam > $O/stats.txt.tmp
mv $O/stats.txt.tmp $O/stats.txt
cat $O/gamcompare.txt
"""


def density_script(c, reads):
    o = results_dir(c, reads) / "density"
    log = DATA_DIR / "logs" / f"density_{Path(reads).stem}_{c['name']}.%j.log"
    # stage.minimizer.results: seeds the scheme picks in the read (before hit-cap filtering)
    # stage.seed.results: seeds giraffe actually uses after filtering
    jq_filter = ('[.name, (.sequence|length), (.annotation.stage.minimizer.results // "NA"), '
                 '(.annotation.stage.seed.results // "NA")] | @tsv')
    return sbatch_header(f"dens-{c['name']}", MAP_MEM, MAP_TIME, log) + require_index(c) + f"""command -v jq > /dev/null || {{ echo "jq not found on $(hostname)" >&2; exit 1; }}
O={o}
mkdir -p $O
{giraffe_cmd(c, reads, "--track-provenance")} > $O/prov.gam.tmp 2> $O/giraffe.log
mv $O/prov.gam.tmp $O/prov.gam
# Only names, lengths and counts leave jq; read sequences are never printed.
{{ printf 'read\\tlength\\tminimizers\\tseeds\\n'
  {c['vg_binary']} view -aj $O/prov.gam | jq -r '{jq_filter}'; }} > $O/per_read.tsv.tmp
mv $O/per_read.tsv.tmp $O/per_read.tsv
awk -F'\\t' -v name={c['name']} 'NR > 1 && $3 != "NA" {{
    n++; len += $2; m += $3; s += $4; d += $3 / $2 }}
    END {{ OFS = "\\t"
        print "name", "reads", "total_bp", "minimizers", "seeds", "density_pooled", "density_mean_per_read", "kept_minimizer_hits_per_bp"
        print name, n, len, m, s, m / len, d / n, s / len }}' $O/per_read.tsv > $O/density_summary.tsv
cat $O/density_summary.tsv
"""


def filters_script(conditions, reads):
    """One small job for all conditions: per-filter totals from the density run's provenance GAM."""
    reads_set = Path(reads).stem
    out_dir = DATA_DIR / "results" / reads_set
    log = DATA_DIR / "logs" / f"filters_{reads_set}.%j.log"
    # Minimizer filters are the funnel filters of stage "minimizer", numbered in the order they run;
    # the last one's passed count is the number of minimizers giraffe keeps (NA if the read has no
    # filter record). stage.seed.results counts seeds, i.e. graph hits of the kept minimizers.
    jq_filter = (
        '(.annotation.filter // {} | to_entries | map(.value + {order: (.key | tonumber)})) as $f'
        ' | (["R", (.annotation.stage.minimizer.results // 0),'
        ' ([$f[] | select(.stage == "minimizer")] | sort_by(.order) | last | .passed.count_total // "NA"),'
        ' (.annotation.stage.seed.results // 0)] | @tsv),'
        ' ($f[] | ["F", .order, .name, .stage, .passed.count_total, .failed.count_total] | @tsv)')
    body = ""
    for c in conditions:
        d = results_dir(c, reads) / "density"
        body += f"""
D={d}
if [ -s $D/prov.gam ]; then
    # Only counts leave jq; read sequences are never printed.
    {c['vg_binary']} view -aj $D/prov.gam | jq -r '{jq_filter}' > $D/filter_lines.tmp
    awk -F'\\t' -v OFS='\\t' -v name={c['name']} '$1 == "F" {{
            key = $2 SUBSEP $3 SUBSEP $4; if (!(key in passed)) order[++n] = key
            passed[key] += $5; failed[key] += $6 }}
        END {{ for (i = 1; i <= n; i++) {{ split(order[i], k, SUBSEP)
            print name, k[1], k[2], k[3], passed[order[i]], failed[order[i]] }} }}' $D/filter_lines.tmp \\
        | sort -t$'\\t' -k2,2n > $D/filters.tsv
    awk -F'\\t' -v OFS='\\t' -v name={c['name']} '$1 == "R" && $3 != "NA" {{ reads++; picked += $2; kept += $3; hits += $4 }}
        END {{ print "name", "reads", "picked_minimizers", "kept_minimizers", "kept_minimizer_hits",
                     "kept_minimizers_per_read", "mean_hits_per_kept_minimizer"
               print name, reads, picked, kept, hits, kept / reads, hits / kept }}' $D/filter_lines.tmp > $D/minimizers.tsv
    rm $D/filter_lines.tmp
    cat $D/filters.tsv >> $OUT/filter_summary.tsv.tmp
    tail -n +2 $D/minimizers.tsv >> $OUT/minimizer_summary.tsv.tmp
else
    echo "skipping {c['name']}: no $D/prov.gam (run the density stage first)" >&2
fi
"""
    return sbatch_header(f"filters-{reads_set}", "4G", "00:30:00", log, partition="short", cpus=1) + f"""command -v jq > /dev/null || {{ echo "jq not found on $(hostname)" >&2; exit 1; }}
OUT={out_dir}
printf 'name\\torder\\tfilter\\tstage\\tpassed\\tfailed\\n' > $OUT/filter_summary.tsv.tmp
printf 'name\\treads\\tpicked_minimizers\\tkept_minimizers\\tkept_minimizer_hits\\tkept_minimizers_per_read\\tmean_hits_per_kept_minimizer\\n' \\
    > $OUT/minimizer_summary.tsv.tmp
{body}
mv $OUT/filter_summary.tsv.tmp $OUT/filter_summary.tsv
mv $OUT/minimizer_summary.tsv.tmp $OUT/minimizer_summary.tsv
column -t -s$'\\t' $OUT/minimizer_summary.tsv
"""


STAGES = {"index": index_script, "map": map_script, "density": density_script}

RESULT_COLUMNS = [
    "name", "label", "density_group", "reads_set",
    "total_reads", "aligned", "eligible", "correct", "wrong", "accuracy",
    "reads_per_sec_per_thread", "giraffe_mem_gb", "giraffe_wall_s",
    "index_min_gb", "index_zipcodes_gb", "index_build_min", "index_peak_gib",
    "occurrences", "density_graph", "density_reads", "density_reads_source",
    "kept_minimizers_per_read", "mean_hits_per_kept_minimizer",
]


def grep_first(path, pattern, cast=float):
    """Return the first capture group of pattern in a text file, or None."""
    if not path.is_file():
        return None
    match = re.search(pattern, path.read_text(errors="replace"))
    return cast(match.group(1)) if match else None


def wall_seconds(text):
    """Convert /usr/bin/time's h:mm:ss or m:ss.ss to seconds."""
    seconds = 0.0
    for part in text.split(":"):
        seconds = seconds * 60 + float(part)
    return seconds


def latest_index_log(name):
    """Newest index log for this condition that finished writing the index."""
    logs = sorted((DATA_DIR / "logs").glob(f"index_{name}.*.log"), key=lambda p: p.stat().st_mtime)
    done = [p for p in logs if "min size (bytes)" in p.read_text(errors="replace")]
    return done[-1] if done else None


def read_density_summary(path):
    if not path.is_file():
        return None
    with open(path, newline="") as f:
        return next(csv.DictReader(f, delimiter="\t"), None)


def collect_row(c, reads_set, graph_length):
    o = DATA_DIR / "results" / reads_set / c["name"]
    row = {"name": c["name"], "label": c["label"], "density_group": c["density_group"], "reads_set": reads_set}

    row["correct"] = grep_first(o / "gamcompare.txt", r"(\d+) reads correct", int)
    row["eligible"] = grep_first(o / "gamcompare.txt", r"(\d+) reads eligible", int)
    if row["correct"] is not None and row["eligible"]:
        row["wrong"] = row["eligible"] - row["correct"]
        row["accuracy"] = row["correct"] / row["eligible"]
    row["total_reads"] = grep_first(o / "stats.txt", r"Total alignments: (\d+)", int)
    row["aligned"] = grep_first(o / "stats.txt", r"Total aligned: (\d+)", int)
    row["reads_per_sec_per_thread"] = grep_first(o / "giraffe.log", r"Mapping speed: ([\d.]+)")
    row["giraffe_mem_gb"] = grep_first(o / "giraffe.log", r"Memory footprint: ([\d.]+) GB")
    row["giraffe_wall_s"] = grep_first(o / "giraffe.time", r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)",
                                       wall_seconds)

    log = latest_index_log(c["name"])
    if log is not None:
        min_bytes = grep_first(log, r"min size \(bytes\): (\d+)", int)
        zip_bytes = grep_first(log, r"zipcodes size \(bytes\): (\d+)", int)
        row["index_min_gb"] = min_bytes / 1e9 if min_bytes is not None else None
        row["index_zipcodes_gb"] = zip_bytes / 1e9 if zip_bytes is not None else None
        build_s = grep_first(log, r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", wall_seconds)
        row["index_build_min"] = build_s / 60 if build_s is not None else None
        peak_kb = grep_first(log, r"Maximum resident set size \(kbytes\): (\d+)", int)
        row["index_peak_gib"] = peak_kb / 2**20 if peak_kb is not None else None
        row["occurrences"] = grep_first(log, r"Minimizer occurrences: (\d+)", int)
        if row["occurrences"] is not None and graph_length:
            row["density_graph"] = row["occurrences"] / graph_length

    # Seed density is a property of the scheme, so fall back to the default read set's measurement.
    for source in dict.fromkeys([reads_set, DEFAULT_READS.stem]):
        summary = read_density_summary(DATA_DIR / "results" / source / c["name"] / "density" / "density_summary.tsv")
        if summary is not None:
            row["density_reads"] = float(summary["density_pooled"])
            row["density_reads_source"] = source
            break
    for source in dict.fromkeys([reads_set, DEFAULT_READS.stem]):
        summary = read_density_summary(DATA_DIR / "results" / source / c["name"] / "density" / "minimizers.tsv")
        if summary is not None:
            row["kept_minimizers_per_read"] = float(summary["kept_minimizers_per_read"])
            row["mean_hits_per_kept_minimizer"] = float(summary["mean_hits_per_kept_minimizer"])
            break
    return row


def format_value(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def wilson_interval(successes, n, z=1.96):
    """95% Wilson score interval for a binomial proportion."""
    p = successes / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half


def plot_comparison(rows, path, reads_set):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    test_color, reference_color = "#2a78d6", "#8a8984"
    ink, muted = "#0b0b0b", "#52514e"
    rows = [r for r in rows if r.get("accuracy") is not None]
    names = [r["name"] for r in rows]
    colors = [reference_color if r["density_group"] == "reference" else test_color for r in rows]
    x = range(len(rows))

    fig, (ax_err, ax_speed) = plt.subplots(1, 2, figsize=(10, 4.2))
    error = [100 * r["wrong"] / r["eligible"] for r in rows]
    bounds = [wilson_interval(r["wrong"], r["eligible"]) for r in rows]
    yerr = [[e - 100 * lo for e, (lo, hi) in zip(error, bounds)],
            [100 * hi - e for e, (lo, hi) in zip(error, bounds)]]
    ax_err.bar(x, error, color=colors, width=0.6)
    ax_err.errorbar(x, error, yerr=yerr, fmt="none", ecolor=muted, elinewidth=1.5, capsize=4)
    for i, r in zip(x, rows):
        ax_err.annotate(f"{r['wrong']}/{r['eligible']}", (i, 100 * bounds[i][1]), xytext=(0, 4),
                        textcoords="offset points", ha="center", fontsize=9, color=ink)
    ax_err.set_title("Mismapped reads (% of eligible, 95% CI)", fontsize=11, color=ink, loc="left")
    # Leave room above the tallest CI whisker for its count label.
    ax_err.set_ylim(0, 1.2 * max(100 * hi for lo, hi in bounds))

    speed = [r["reads_per_sec_per_thread"] or 0 for r in rows]
    ax_speed.bar(x, speed, color=colors, width=0.6)
    for i, s in zip(x, speed):
        ax_speed.annotate(f"{s:.1f}", (i, s), xytext=(0, 4), textcoords="offset points",
                          ha="center", fontsize=9, color=ink)
    ax_speed.set_title("Mapping speed (reads / s / thread)", fontsize=11, color=ink, loc="left")
    ax_speed.set_ylim(0, 1.15 * max(speed))

    for ax in (ax_err, ax_speed):
        ax.set_xticks(list(x))
        ax.set_xticklabels(names, rotation=20, ha="right", fontsize=9, color=ink)
        ax.set_facecolor("#fcfcfb")
        ax.grid(axis="y", color="#e6e5e0", linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color("#c3c2b7")
        ax.tick_params(colors=muted)
    fig.suptitle(f"{reads_set}  (blue: tested conditions, gray: lab default reference)",
                 fontsize=10, color=muted, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor="#fcfcfb")
    plt.close(fig)


def collect(conditions, reads):
    reads_set = Path(reads).stem
    graph_length = grep_first(DATA_DIR / "graph.length.txt", r"length\s+(\d+)", int)
    rows, skipped = [], []
    for c in conditions:
        row = collect_row(c, reads_set, graph_length)
        (rows if row.get("correct") is not None else skipped).append(row)
    out_dir = DATA_DIR / "results" / reads_set
    out_dir.mkdir(parents=True, exist_ok=True)
    table = out_dir / "results.tsv"
    with open(table, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(RESULT_COLUMNS)
        for row in rows:
            writer.writerow([format_value(row.get(col)) for col in RESULT_COLUMNS])
    print(f"wrote {table} ({len(rows)} conditions)")
    if skipped:
        print(f"skipped (no map results for {reads_set}): {[r['name'] for r in skipped]}")
    if rows:
        plot = out_dir / "comparison.png"
        plot_comparison(rows, plot, reads_set)
        print(f"wrote {plot}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("stage", choices=sorted(STAGES) + ["filters", "collect"])
    parser.add_argument("--conditions", default="conditions.tsv")
    parser.add_argument("--reads", default=str(DEFAULT_READS),
                        help="truth GAM for map/density (default: 1k HG002 HiFi)")
    parser.add_argument("--only", nargs="*", help="run only these condition names")
    parser.add_argument("--submit", action="store_true", help="sbatch the generated scripts")
    args = parser.parse_args()

    conditions = load_conditions(args.conditions)
    if args.only:
        unknown = set(args.only) - {c["name"] for c in conditions}
        if unknown:
            sys.exit(f"unknown condition names: {sorted(unknown)}")
        conditions = [c for c in conditions if c["name"] in args.only]

    if args.stage == "collect":
        collect(conditions, args.reads)
        return
    if args.stage == "filters":
        (DATA_DIR / "logs").mkdir(parents=True, exist_ok=True)
        script = DATA_DIR / "jobs" / f"filters_{Path(args.reads).stem}.sbatch"
        script.write_text(filters_script(conditions, args.reads))
        print(f"[{len(conditions)} conditions] filters -> {script}")
        if args.submit:
            result = subprocess.run(["sbatch", str(script)], check=True, capture_output=True, text=True)
            print("  " + result.stdout.strip())
        return

    jobs_dir = DATA_DIR / "jobs"
    jobs_dir.mkdir(parents=True, exist_ok=True)
    (DATA_DIR / "logs").mkdir(parents=True, exist_ok=True)

    tag = "" if args.stage == "index" else f"_{Path(args.reads).stem}"
    for c in conditions:
        script = jobs_dir / f"{args.stage}{tag}_{c['name']}.sbatch"
        script.write_text(STAGES[args.stage](c, args.reads))
        print(f"[{c['name']}] {args.stage} -> {script}")
        if args.submit:
            result = subprocess.run(["sbatch", str(script)], check=True, capture_output=True, text=True)
            print("  " + result.stdout.strip())


if __name__ == "__main__":
    main()
