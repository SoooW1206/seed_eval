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
  map_mem        sbatch --mem for this condition's map/density jobs (empty: MAP_MEM)
  density_group  only used for result tables and plots
  label          only used for result tables and plots

Stages:
  index    vg minimizer -> indexes/<name>.min + .zipcodes
  map      giraffe -> annotate -> gamcompare -> stats, on --reads; --reps N runs it
           N times into <name>/rep1..N (for timing), keeping the GAMs of rep1 only
  density  giraffe --track-provenance on --reads; per-read seed counts
  filters  one small job: per-filter pass/fail totals and kept minimizers from the density run
  collect  (runs locally, no sbatch) parse logs of all stages into
           results/<reads>/results.tsv and results/<reads>/comparison.png;
           --snapshot also copies the summaries to results_summary/<reads>/ in this repo

Index provenance: each index job writes indexes/<name>.params (vg binary and
version, index_args, graph). index/map/density refuse to run, both when
generating scripts and again inside the job, if that file disagrees with the
condition's row, so a row whose index_args changed without a new name cannot
silently reuse an old index.
"""
import argparse
import csv
import functools
import math
import re
import shutil
import statistics
import subprocess
import sys
from datetime import datetime
from pathlib import Path

REPO_DIR = Path(__file__).resolve().parent
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
OPTIONAL_COLUMNS = ["giraffe_args", "map_mem", "density_group", "label"]


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


@functools.lru_cache(maxsize=None)
def vg_version(vg_binary):
    """First line of `vg version`, e.g. 'vg version v1.77.0-26-gc91cc1d91 "Ruby"'."""
    result = subprocess.run([vg_binary, "version"], check=True, capture_output=True, text=True)
    return result.stdout.splitlines()[0].strip()


@functools.lru_cache(maxsize=None)
def git_commit():
    """Short HEAD commit of this repo, with -dirty if the driver or conditions table has uncommitted edits."""
    def git(*args):
        return subprocess.run(["git", "-C", str(REPO_DIR), *args], check=True, capture_output=True, text=True).stdout
    dirty = git("status", "--porcelain", "--", "seed_eval.py", "conditions.tsv").strip()
    return git("rev-parse", "--short", "HEAD").strip() + ("-dirty" if dirty else "")


def index_params(c):
    """Expected contents of indexes/<name>.params for this condition."""
    return (f"vg_binary={c['vg_binary']}\n"
            f"vg_version={vg_version(c['vg_binary'])}\n"
            f"index_args={c['index_args']}\n"
            f"gbz={GBZ}\n"
            f"dist={DIST}\n")


def check_index_params(c, stage, force):
    """Compare indexes/<name>.params with the row. Returns False if the index stage can be skipped."""
    out = index_prefix(c)
    params = Path(f"{out}.params")
    if params.is_file():
        if params.read_text() != index_params(c):
            sys.exit(f"{c['name']}: {params} does not match conditions.tsv.\n"
                     f"  on disk:\n{params.read_text()}  expected:\n{index_params(c)}"
                     "Give the changed row a new name, or delete the old index files.")
        if stage == "index" and not force:
            print(f"[{c['name']}] index already built with these parameters; skipping (use --force to rebuild)")
            return False
        return True
    if stage == "index":
        if Path(f"{out}.min").exists() and not force:
            sys.exit(f"{c['name']}: {out}.min exists but has no .params; cannot tell how it was built. "
                     "Delete it or pass --force to rebuild.")
        return True
    sys.exit(f"{c['name']}: no {params}; run the index stage first")


def require_index(c):
    """Job-time check that the index on disk was built with this row's parameters."""
    out = index_prefix(c)
    return f"""EXPECTED_PARAMS=$(cat <<'PARAMS'
{index_params(c)}PARAMS
)
[ "$(cat {out}.params 2>/dev/null)" = "$EXPECTED_PARAMS" ] || {{
    echo "{out}.params missing or does not match this condition; rebuild or rename the index" >&2; exit 1; }}
for f in {out}.min {out}.zipcodes; do
    [ -s "$f" ] || {{ echo "missing index file $f; run the index stage first" >&2; exit 1; }}
done
"""


def run_info(c, o):
    """Shell lines recording which code and arguments produced the results in directory o."""
    return f"""printf 'git_commit\\trun_date\\thost\\tjob_id\\tgiraffe_args\\n%s\\t%s\\t%s\\t%s\\t%s\\n' \\
    '{git_commit()}' "$(date -Iseconds)" "$(hostname)" "$SLURM_JOB_ID" '{c['giraffe_args']}' > {o}/run_info.tsv
"""


def giraffe_cmd(c, reads, extra=""):
    out = index_prefix(c)
    return (f"{c['vg_binary']} giraffe -p -b {GIRAFFE_PRESET} -t {THREADS} {c['giraffe_args']} {extra}"
            f" \\\n    -Z {GBZ} -d {DIST} -m {out}.min -z {out}.zipcodes -G {reads}")


def index_script(c, reads):
    out = index_prefix(c)
    log = DATA_DIR / "logs" / f"index_{c['name']}.%j.log"
    return sbatch_header(f"idx-{c['name']}", c["index_mem"], INDEX_TIME, log) + f"""mkdir -p {out.parent}
# Until the new .params is written, map/density treat this index as not built.
rm -f {out}.params
/usr/bin/time -v {c['vg_binary']} minimizer {c['index_args']} -p -t {THREADS} \\
    -d {DIST} -o {out}.min.tmp -z {out}.zipcodes.tmp {GBZ}
mv {out}.min.tmp {out}.min
mv {out}.zipcodes.tmp {out}.zipcodes
echo "min size (bytes): $(stat -c %s {out}.min)"
echo "zipcodes size (bytes): $(stat -c %s {out}.zipcodes)"
cat > {out}.params <<'PARAMS'
{index_params(c)}PARAMS
"""


def results_dir(c, reads):
    return DATA_DIR / "results" / Path(reads).stem / c["name"]


def map_mem(c):
    return c["map_mem"] or MAP_MEM


def map_script(c, reads, rep=None):
    """rep None: single run into <name>/. rep i: replicate i into <name>/rep<i>/; GAMs kept for rep 1 only."""
    o = results_dir(c, reads)
    suffix = ""
    if rep is not None:
        o = o / f"rep{rep}"
        suffix = f"_rep{rep}"
    log = DATA_DIR / "logs" / f"map_{Path(reads).stem}_{c['name']}{suffix}.%j.log"
    vg = c["vg_binary"]
    cleanup = ""
    if rep is not None and rep > 1:
        cleanup = "# Replicates beyond the first are for timing; keep only their statistics.\nrm $O/mapped.gam $O/annot.gam\n"
    return sbatch_header(f"map-{c['name']}{suffix}", map_mem(c), MAP_TIME, log) + require_index(c) + f"""O={o}
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
""" + cleanup + run_info(c, o) + """cat $O/gamcompare.txt
"""


def density_script(c, reads):
    o = results_dir(c, reads) / "density"
    log = DATA_DIR / "logs" / f"density_{Path(reads).stem}_{c['name']}.%j.log"
    # stage.minimizer.results: seeds the scheme picks in the read (before hit-cap filtering)
    # stage.seed.results: seeds giraffe actually uses after filtering
    jq_filter = ('[.name, (.sequence|length), (.annotation.stage.minimizer.results // "NA"), '
                 '(.annotation.stage.seed.results // "NA")] | @tsv')
    return sbatch_header(f"dens-{c['name']}", map_mem(c), MAP_TIME, log) + require_index(c) + f"""command -v jq > /dev/null || {{ echo "jq not found on $(hostname)" >&2; exit 1; }}
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
""" + run_info(c, o) + """cat $O/density_summary.tsv
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
    # Seeds per read (stage.seed.results) from per_read.tsv: median, mean, nearest-rank 90th percentile.
    read SEED_MEDIAN SEED_MEAN SEED_P90 < <(awk -F'\\t' 'NR > 1 && $4 != "NA" {{ print $4 }}' $D/per_read.tsv | sort -n \\
        | awk '{{ a[NR] = $1; s += $1 }} END {{ k = int(0.9 * NR); if (k < 0.9 * NR) k++
              print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2), s / NR, a[k] }}')
    awk -F'\\t' -v OFS='\\t' -v name={c['name']} -v med="$SEED_MEDIAN" -v mean="$SEED_MEAN" -v p90="$SEED_P90" \\
        '$1 == "R" && $3 != "NA" {{ reads++; picked += $2; kept += $3; hits += $4 }}
        END {{ print "name", "reads", "picked_minimizers", "kept_minimizers", "kept_minimizer_hits",
                     "kept_minimizers_per_read", "mean_hits_per_kept_minimizer",
                     "seeds_per_read_median", "seeds_per_read_mean", "seeds_per_read_p90"
               print name, reads, picked, kept, hits, kept / reads, hits / kept, med, mean, p90 }}' $D/filter_lines.tmp > $D/minimizers.tsv
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
printf 'name\\treads\\tpicked_minimizers\\tkept_minimizers\\tkept_minimizer_hits\\tkept_minimizers_per_read\\tmean_hits_per_kept_minimizer\\tseeds_per_read_median\\tseeds_per_read_mean\\tseeds_per_read_p90\\n' \\
    > $OUT/minimizer_summary.tsv.tmp
{body}
mv $OUT/filter_summary.tsv.tmp $OUT/filter_summary.tsv
mv $OUT/minimizer_summary.tsv.tmp $OUT/minimizer_summary.tsv
column -t -s$'\\t' $OUT/minimizer_summary.tsv
"""


STAGES = {"index": index_script, "map": map_script, "density": density_script}

RESULT_COLUMNS = [
    "name", "label", "density_group", "reads_set", "index_args", "giraffe_args", "git_commit", "run_date",
    "total_reads", "aligned", "eligible", "correct", "wrong", "accuracy",
    "reads_per_sec_per_thread", "speed_min", "speed_max", "reps", "reps_consistent",
    "giraffe_mem_gb", "giraffe_wall_s",
    "index_min_gb", "index_zipcodes_gb", "index_build_min", "index_peak_gib",
    "occurrences", "density_graph", "density_reads",
    "seeds_per_read_median", "seeds_per_read_mean", "seeds_per_read_p90",
    "kept_minimizers_per_read", "mean_hits_per_kept_minimizer",
]


SNAPSHOT_FILES = ["results.tsv", "comparison.png", "filter_summary.tsv", "minimizer_summary.tsv"]


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
    # Replicated runs (map --reps) live in rep1..N: accuracy, memory and provenance come from
    # rep1, speed is the median over replicates.
    runs = sorted((p for p in o.glob("rep*") if p.is_dir() and re.fullmatch(r"rep\d+", p.name)
                   and (p / "gamcompare.txt").is_file()), key=lambda p: int(p.name[3:]))
    if runs:
        o = runs[0]
    else:
        runs = [o]
    row = {"name": c["name"], "label": c["label"], "density_group": c["density_group"], "reads_set": reads_set}

    row["correct"] = grep_first(o / "gamcompare.txt", r"(\d+) reads correct", int)
    row["eligible"] = grep_first(o / "gamcompare.txt", r"(\d+) reads eligible", int)
    if row["correct"] is not None and row["eligible"]:
        row["wrong"] = row["eligible"] - row["correct"]
        row["accuracy"] = row["correct"] / row["eligible"]
    # Provenance of the map run. Runs from before run_info.tsv existed get commit NA and the
    # gamcompare output's modification time as an approximate run date.
    info = read_density_summary(o / "run_info.tsv")
    params = Path(f"{index_prefix(c)}.params")
    if params.is_file():
        match = re.search(r"^index_args=(.*)$", params.read_text(), re.M)
        row["index_args"] = match.group(1) if match else None
    if info is not None:
        row["git_commit"], row["run_date"], row["giraffe_args"] = info["git_commit"], info["run_date"], info["giraffe_args"]
    elif (o / "gamcompare.txt").is_file():
        row["giraffe_args"] = c["giraffe_args"]
        mtime = (o / "gamcompare.txt").stat().st_mtime
        row["run_date"] = datetime.fromtimestamp(mtime).astimezone().isoformat(timespec="seconds")
    row["total_reads"] = grep_first(o / "stats.txt", r"Total alignments: (\d+)", int)
    row["aligned"] = grep_first(o / "stats.txt", r"Total aligned: (\d+)", int)
    speeds = [s for s in (grep_first(r / "giraffe.log", r"Mapping speed: ([\d.]+)") for r in runs) if s is not None]
    if speeds:
        row["reads_per_sec_per_thread"] = statistics.median(speeds)
        row["speed_min"], row["speed_max"] = min(speeds), max(speeds)
    row["reps"] = len(speeds)
    # Giraffe should be deterministic, so every replicate must get the same number of reads right.
    corrects = {grep_first(r / "gamcompare.txt", r"(\d+) reads correct", int) for r in runs}
    row["reps_consistent"] = "yes" if len(corrects) == 1 else "no"
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

    # Seed statistics come from the density run on this same read set only; without it they stay NA.
    density_dir = DATA_DIR / "results" / reads_set / c["name"] / "density"
    summary = read_density_summary(density_dir / "density_summary.tsv")
    if summary is not None:
        row["density_reads"] = float(summary["density_pooled"])
        row.update(seed_stats(density_dir / "per_read.tsv"))
    summary = read_density_summary(density_dir / "minimizers.tsv")
    if summary is not None:
        row["kept_minimizers_per_read"] = float(summary["kept_minimizers_per_read"])
        row["mean_hits_per_kept_minimizer"] = float(summary["mean_hits_per_kept_minimizer"])
    return row


def seed_stats(per_read):
    """Median, mean and nearest-rank 90th percentile of seeds per read (stage.seed.results)."""
    with open(per_read, newline="") as f:
        seeds = sorted(int(r["seeds"]) for r in csv.DictReader(f, delimiter="\t") if r["seeds"] != "NA")
    if not seeds:
        return {}
    return {"seeds_per_read_median": statistics.median(seeds), "seeds_per_read_mean": statistics.mean(seeds),
            "seeds_per_read_p90": seeds[math.ceil(0.9 * len(seeds)) - 1]}


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

    # The seeds panel needs this read set's density runs; leave it out if none exist.
    has_seeds = any(r.get("seeds_per_read_median") is not None for r in rows)
    fig, axes = plt.subplots(1, 3 if has_seeds else 2, figsize=(14 if has_seeds else 10, 4.4))
    ax_err, ax_speed = axes[0], axes[1]
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
    ax_speed.set_title("Mapping speed (reads / s / thread, median of reps)", fontsize=11, color=ink, loc="left")
    ax_speed.set_ylim(0, 1.15 * max(speed))

    if has_seeds:
        # Solid bar: median seeds per read; light bar: mean (pulled up by a few repetitive reads).
        ax_seeds = axes[2]
        width = 0.38
        for i, r in zip(x, rows):
            for offset, stat, alpha in ((-width / 2, "median", 1.0), (width / 2, "mean", 0.4)):
                value = r.get(f"seeds_per_read_{stat}")
                if value is None:
                    continue
                ax_seeds.bar(i + offset, value, width=width, color=colors[i], alpha=alpha,
                             edgecolor="#fcfcfb", linewidth=2)
                ax_seeds.annotate(f"{value:.0f}", (i + offset, value), xytext=(0, 3), textcoords="offset points",
                                  ha="center", fontsize=7, color=ink)
        top = max(r.get("seeds_per_read_mean") or 0 for r in rows)
        ax_seeds.set_ylim(0, 1.15 * top if top > 0 else 1)
        ax_seeds.set_title(f"Seeds per read on {reads_set}\n(solid: median, light: mean)",
                           fontsize=11, color=ink, loc="left")

    for ax in axes:
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


def collect(conditions, reads, snapshot=False):
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
    no_density = [r["name"] for r in rows if r.get("seeds_per_read_median") is None]
    if no_density:
        print(f"no density run on {reads_set} for {no_density}: their seed columns are NA"
              + ("" if len(no_density) < len(rows) else " and the seeds panel is left out"))
    inconsistent = [r["name"] for r in rows if r.get("reps_consistent") == "no"]
    if inconsistent:
        print(f"WARNING: replicates disagree on the number of correct reads: {inconsistent}")
    if rows:
        plot = out_dir / "comparison.png"
        plot_comparison(rows, plot, reads_set)
        print(f"wrote {plot}")
    if snapshot:
        dest = REPO_DIR / "results_summary" / reads_set
        dest.mkdir(parents=True, exist_ok=True)
        for name in SNAPSHOT_FILES:
            if (out_dir / name).is_file():
                shutil.copy2(out_dir / name, dest / name)
                print(f"snapshot {dest / name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("stage", choices=sorted(STAGES) + ["filters", "collect"])
    parser.add_argument("--conditions", default="conditions.tsv")
    parser.add_argument("--reads", default=str(DEFAULT_READS),
                        help="truth GAM for map/density (default: 1k HG002 HiFi)")
    parser.add_argument("--only", nargs="*", help="run only these condition names")
    parser.add_argument("--submit", action="store_true", help="sbatch the generated scripts")
    parser.add_argument("--force", action="store_true", help="index: rebuild even if already built")
    parser.add_argument("--reps", type=int, default=1, help="map: number of replicate runs (default 1)")
    parser.add_argument("--snapshot", action="store_true",
                        help="collect: also copy the summaries to results_summary/<reads>/ in this repo")
    args = parser.parse_args()
    if args.reps < 1 or (args.reps > 1 and args.stage != "map"):
        sys.exit("--reps must be >= 1 and only applies to the map stage")

    conditions = load_conditions(args.conditions)
    if args.only:
        unknown = set(args.only) - {c["name"] for c in conditions}
        if unknown:
            sys.exit(f"unknown condition names: {sorted(unknown)}")
        conditions = [c for c in conditions if c["name"] in args.only]

    if args.stage == "collect":
        collect(conditions, args.reads, args.snapshot)
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
    reps = [None] if args.reps == 1 else list(range(1, args.reps + 1))
    for c in conditions:
        if not check_index_params(c, args.stage, args.force):
            continue
        for rep in reps:
            rep_tag = "" if rep is None else f"_rep{rep}"
            script = jobs_dir / f"{args.stage}{tag}_{c['name']}{rep_tag}.sbatch"
            text = map_script(c, args.reads, rep) if args.stage == "map" else STAGES[args.stage](c, args.reads)
            script.write_text(text)
            print(f"[{c['name']}] {args.stage}{rep_tag} -> {script}")
            if args.submit:
                result = subprocess.run(["sbatch", str(script)], check=True, capture_output=True, text=True)
                print("  " + result.stdout.strip())


if __name__ == "__main__":
    main()
