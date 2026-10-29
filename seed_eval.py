#!/usr/bin/env python3
"""Seed-method evaluator for vg giraffe.

Reads a conditions table (one seed-index configuration per row) and generates
one sbatch job per condition for the chosen stage. Scripts are only written and
printed by default; pass --submit to actually submit them.

The driver is data-driven: every condition-specific piece of a command comes
from its row in the table, so adding a condition means adding a row.

Columns:
  name           unique condition name; used for output file names
  vg_binary      absolute path to the vg executable used for every step, or
                 `default` for the vg in the config file
  index_args     arguments passed verbatim to `vg minimizer`
  index_mem      sbatch --mem for the index job (e.g. 80G). Measured at k=31 on
                 the HG002 sampled graph: plain ~56 GiB -> 80G, -W ~150 GiB -> 256G
  giraffe_args   extra arguments passed verbatim to `vg giraffe` (may be empty)
  map_mem        sbatch --mem for this condition's map/density jobs (empty: MAP_MEM)
  density_group  only used for result tables and plots
  label          only used for result tables and plots

Paths (data directory, graph, reads directory, default vg) are read from a
tab-separated config file (--config, default config.tsv next to this script;
see config.example.tsv).

Stages:
  prepare  one-time jobs per graph: check reference paths, build the distance
           index (<graph>.dist next to the .gbz), total graph length
  index    vg minimizer -> indexes/<name>.min + .zipcodes
  map      giraffe -> annotate -> gamcompare -> stats, on --reads; --reps N runs it
           N times into <name>/rep1..N (for timing), keeping the GAMs of rep1 only
           --run-tag TAG writes <name>/TAG/rep1..N instead; --keep-gams keeps every replicate's GAMs
  density  giraffe --track-provenance on --reads; per-read seed counts
  fingerprint  one small job per run with a mapped.gam: per-read MAPQ, score, identity, start and a
           hash of the path (no read bases), into fingerprint.tsv
  determinism  (runs locally) compare fingerprints of the --run-tag replicates by read name, and
           separately together with the earlier rep1; writes results/<reads>/determinism/
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
import collections
import csv
import functools
import hashlib
import math
import re
import shlex
import shutil
import statistics
import subprocess
import sys
from datetime import datetime
from pathlib import Path

REPO_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = REPO_DIR / "config.tsv"
CONFIG_KEYS = ["data_dir", "graph", "reads_dir", "vg"]
# Optional config key: exclude_nodes, a Slurm node list that --fixed-nodes passes to sbatch --exclude
# (to keep runs on one node type).
DEFAULT_READS_NAME = "HG002-sim-hifi-1k.gam"  # default truth GAM, inside reads_dir
# Paths come from the config file; load_config() sets them before anything else runs.
DATA_DIR = GBZ = DIST = READS_DIR = DEFAULT_READS = VG = None
EXCLUDE_NODES = ""
PARTITION = "medium"
THREADS = 16  # vg minimizer: more than 16 threads rarely helps
GIRAFFE_PRESET = "hifi"
COMPARE_RANGE = 200  # bp window for gamcompare to call a read correct
MAP_MEM = "64G"  # measured on 1k reads: giraffe footprint ~41 GB, job peak ~36 GiB
MAP_TIME = "02:00:00"
INDEX_TIME = "12:00:00"
TIMING_ROUNDS = 3  # timing: each job maps every condition this many times, interleaved
TIMING_JOBS = 2  # timing: identical exclusive jobs submitted, ideally landing on different nodes
TIMING_TIME = "04:00:00"  # 12 runs of ~8 min on 100k reads
# timing --real (full real-read FASTQ, lab practice for real-read speed): one round per job, order reversed
# in the second job. Threads and memory as the lab's MAPPER_THREADS and auto_mapping_memory for hifi.
REAL_THREADS = 64
REAL_MEM = "300000M"  # 240000 MB / 64 * threads + 60000 MB
REAL_PARTITION = "long"
REAL_TIME = "30:00:00"  # lab budgets 20 h for one full real-read run; ~6-10 h expected for 4

REQUIRED_COLUMNS = ["name", "vg_binary", "index_args", "index_mem"]
OPTIONAL_COLUMNS = ["giraffe_args", "map_mem", "density_group", "label"]


def load_config(path):
    """Read key/value paths (tab-separated) and set DATA_DIR, GBZ, DIST, READS_DIR, DEFAULT_READS and VG.
    The distance index is expected next to the graph, with .dist in place of .gbz."""
    global DATA_DIR, GBZ, DIST, READS_DIR, DEFAULT_READS, VG, EXCLUDE_NODES
    if not Path(path).is_file():
        sys.exit(f"{path}: config file not found; copy config.example.tsv to config.tsv and fill in your paths")
    config = {}
    with open(path, newline="") as f:
        for row in csv.reader(f, delimiter="\t"):
            if row and row[0].strip() and not row[0].startswith("#") and row[0].strip() != "key":
                config[row[0].strip()] = row[1].strip() if len(row) > 1 else ""
    missing = [key for key in CONFIG_KEYS if not config.get(key)]
    if missing:
        sys.exit(f"{path}: missing or empty keys {missing}")
    GBZ = Path(config["graph"])
    if GBZ.suffix != ".gbz":
        sys.exit(f"{path}: graph must be a .gbz file, got {GBZ}")
    DATA_DIR = Path(config["data_dir"])
    DIST = GBZ.with_suffix(".dist")
    READS_DIR = Path(config["reads_dir"])
    DEFAULT_READS = READS_DIR / DEFAULT_READS_NAME
    VG = config["vg"]
    EXCLUDE_NODES = config.get("exclude_nodes", "")


def load_conditions(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = set(REQUIRED_COLUMNS + OPTIONAL_COLUMNS) - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"{path}: missing columns {sorted(missing)}")
        conditions = [{k: (v or "").strip() for k, v in row.items()} for row in reader]
    validate(conditions, path)
    for c in conditions:
        if c["vg_binary"] == "default":
            c["vg_binary"] = VG
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


def sbatch_header(job_name, mem, time, log, partition=PARTITION, cpus=THREADS, extra=""):
    return f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH -c {cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --output={log}
{extra}set -euo pipefail
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
        missing = [str(record) for record, _, _ in checksum_files(c) if not record.is_file()]
        if stage not in ("index", "checksums") and missing:
            sys.exit(f"{c['name']}: no {', '.join(missing)}; record it with: "
                     f"python3 seed_eval.py checksums --only {c['name']}")
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
""" + require_checksums(c)


COUNTS_OPTION = "--minimaxmer-counts"


def counts_path(c):
    """The minimaxmer counts file named in index_args, or None."""
    match = re.search(rf"{COUNTS_OPTION}\s+(\S+)", c["index_args"])
    return Path(match.group(1)) if match else None


def checksum_files(c):
    """(record, label, file): indexes/<name>.vg_sha256 for the vg binary, and for a minimaxmer condition
    indexes/<name>.counts_sha256 for its counts file. A record holds <label>=<path>, bytes= and sha256=.
    Kept apart from .params, so that adding them changes no .params that queued jobs compare against."""
    out = index_prefix(c)
    files = [(Path(f"{out}.vg_sha256"), "vg_binary", Path(c["vg_binary"]))]
    if counts_path(c):
        files.append((Path(f"{out}.counts_sha256"), "counts", counts_path(c)))
    return files


def write_checksums(c):
    """Index-job lines writing the checksum records of the files the index was just built with."""
    return "".join(f"""printf '{label}=%s\\nbytes=%s\\nsha256=%s\\n' '{path}' "$(stat -c %s {path})" \\
    "$(sha256sum < {path} | cut -d' ' -f1)" > {record}
""" for record, label, path in checksum_files(c))


def require_checksums(c):
    """Job-time check that the vg binary (and counts file) are the ones recorded when the index was built."""
    return "".join(f"""[ -s {record} ] || {{
    echo "missing {record}; record it with: python3 seed_eval.py checksums --only {c['name']}" >&2; exit 1; }}
[ "$(sed -n 's/^{label}=//p' {record})" = "{path}" ] && \\
    [ "$(sha256sum < {path} | cut -d' ' -f1)" = "$(sed -n 's/^sha256=//p' {record})" ] || {{
    echo "{path} is not the {label} recorded in {record}" >&2; exit 1; }}
""" for record, label, path in checksum_files(c))


@functools.lru_cache(maxsize=None)
def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def record_checksums(conditions, force):
    """checksums stage: write missing checksum records for indexes that are already built. A file is only
    recorded if it was last modified before the index's .params was written, i.e. the index was built with
    this very file; otherwise the index has to be rebuilt."""
    for c in conditions:
        params = Path(f"{index_prefix(c)}.params")
        if not params.is_file() or params.read_text() != index_params(c):
            print(f"[{c['name']}] no matching {params}; skipped (the index stage writes the records)")
            continue
        for record, label, path in checksum_files(c):
            if record.is_file() and not force:
                print(f"[{c['name']}] {record} exists; kept (use --force to rewrite)")
                continue
            if path.stat().st_mtime >= params.stat().st_mtime:
                sys.exit(f"{c['name']}: {path} changed after {params} was written; rebuild the index instead")
            record.write_text(f"{label}={path}\nbytes={path.stat().st_size}\nsha256={sha256_file(path)}\n")
            print(f"[{c['name']}] wrote {record}")


def run_info(c, o):
    """Shell lines recording which code and arguments produced the results in directory o."""
    return f"""printf 'git_commit\\trun_date\\thost\\tjob_id\\tgiraffe_args\\n%s\\t%s\\t%s\\t%s\\t%s\\n' \\
    '{git_commit()}' "$(date -Iseconds)" "$(hostname)" "$SLURM_JOB_ID" '{c['giraffe_args']}' > {o}/run_info.tsv
"""


def giraffe_cmd(c, reads, extra="", threads=THREADS, reads_flag="-G"):
    """reads_flag -G: truth GAM; -f: FASTQ (real reads)."""
    out = index_prefix(c)
    return (f"{c['vg_binary']} giraffe -p -b {GIRAFFE_PRESET} -t {threads} {c['giraffe_args']} {extra}"
            f" \\\n    -Z {GBZ} -d {DIST} -m {out}.min -z {out}.zipcodes {reads_flag} {reads}")


def index_script(c, reads):
    out = index_prefix(c)
    log = DATA_DIR / "logs" / f"index_{c['name']}.%j.log"
    return sbatch_header(f"idx-{c['name']}", c["index_mem"], INDEX_TIME, log) + f"""mkdir -p {out.parent}
# Until the new .params is written, map/density treat this index as not built.
rm -f {out}.params {' '.join(str(record) for record, _, _ in checksum_files(c))}
/usr/bin/time -v {c['vg_binary']} minimizer {c['index_args']} -p -t {THREADS} \\
    -d {DIST} -o {out}.min.tmp -z {out}.zipcodes.tmp {GBZ}
mv {out}.min.tmp {out}.min
mv {out}.zipcodes.tmp {out}.zipcodes
echo "min size (bytes): $(stat -c %s {out}.min)"
echo "zipcodes size (bytes): $(stat -c %s {out}.zipcodes)"
{write_checksums(c)}cat > {out}.params <<'PARAMS'
{index_params(c)}PARAMS
"""


def reads_set_name(reads):
    """Read set name: the file name without .gam, or without .fastq.gz/.fq.gz/.fastq/.fq for real reads."""
    name = Path(reads).name
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(ext):
            return name[:-len(ext)]
    return Path(reads).stem


def results_dir(c, reads):
    return DATA_DIR / "results" / Path(reads).stem / c["name"]


def map_mem(c):
    return c["map_mem"] or MAP_MEM


def map_script(c, reads, rep=None, run_tag=None, keep_gams=False):
    """rep None: single run into <name>/. rep i: replicate i into <name>/rep<i>/; GAMs kept for rep 1 only
    unless keep_gams. run_tag puts the runs under <name>/<run_tag>/ so they never mix with the default runs."""
    o = results_dir(c, reads)
    suffix = ""
    if run_tag:
        o = o / run_tag
        suffix = f"_{run_tag}"
    if rep is not None:
        o = o / f"rep{rep}"
        suffix += f"_rep{rep}"
    log = DATA_DIR / "logs" / f"map_{Path(reads).stem}_{c['name']}{suffix}.%j.log"
    vg = c["vg_binary"]
    cleanup = ""
    if rep is not None and rep > 1 and not keep_gams:
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


TIMING_COLUMNS = ["job_id", "node", "round", "order", "condition", "load_seconds", "map_seconds",
                  "reads_per_sec_per_thread", "reads_per_cpu_second", "cpu_seconds", "memory_gb"]
# timing --real only: peak RSS from /usr/bin/time, the rest from real_stats.tsv (MAPQ and identity over aligned reads)
REAL_COLUMNS = ["peak_rss_gib", "reads", "aligned_fraction", "mapq0_fraction", "mapq60_fraction",
                "identity_mean", "identity_median"]
# Per-read fingerprint of a real-read run from `vg filter -T` (jq on the JSON would take about a day per run):
# start and end are the first and last node with strand, n_nodes the number of nodes on the path. No bases.
REAL_TSV_FIELDS = "name;is_aligned;mapping_quality;score;identity;nodes"
REAL_FINGERPRINT_AWK = r'''BEGIN { FS = OFS = "\t"; print "read", "aligned", "mapq", "score", "identity", "start", "end", "n_nodes" }
$0 == "" || /^#/ { next }
{ n = split($6, nodes, ","); if (n && nodes[n] == "") n--
  print $1, ($2 == "True"), $3, $4, $5, (n ? nodes[1] : "NA"), (n ? nodes[n] : "NA"), n }'''
# Arguments: fingerprint.tsv, output directory, number of reads giraffe reports. Writes real_stats.tsv,
# mapq_hist.tsv and identity_hist.tsv (bins of 0.001); numbers only.
REAL_STATS_PY = r'''import collections, statistics, sys
path, out, expected = sys.argv[1], sys.argv[2], int(sys.argv[3])
names, mapq, identity, reads = set(), collections.Counter(), [], 0
with open(path) as f:
    next(f)
    for line in f:
        read, aligned, q, score, ident = line.split("\t")[:5]
        names.add(read)
        reads += 1
        if aligned == "1":
            mapq[int(q)] += 1
            identity.append(float(ident))
if len(names) != reads:
    sys.exit(f"{reads - len(names)} duplicate read names")
if reads != expected:
    sys.exit(f"fingerprint has {reads} reads, giraffe reports {expected}")
aligned = len(identity)
with open(f"{out}/real_stats.tsv", "w") as f:
    f.write("reads\taligned_fraction\tmapq0_fraction\tmapq60_fraction\tidentity_mean\tidentity_median\n")
    f.write(f"{reads}\t{aligned / reads:.6g}\t{mapq[0] / aligned:.6g}\t{mapq[60] / aligned:.6g}\t"
            f"{statistics.fmean(identity):.6g}\t{statistics.median(identity):.6g}\n")
with open(f"{out}/mapq_hist.tsv", "w") as f:
    f.write("mapq\tcount\n" + "".join(f"{q}\t{n}\n" for q, n in sorted(mapq.items())))
bins = collections.Counter(min(int(i * 1000), 999) for i in identity)
with open(f"{out}/identity_hist.tsv", "w") as f:
    f.write("identity_from\tcount\n" + "".join(f"{b / 1000:.3f}\t{n}\n" for b, n in sorted(bins.items())))'''
PLACEMENT_COLUMNS = ["job_id", "condition", "against", "reads", "aligned_both", "aligned_only_condition",
                     "aligned_only_against", "different_start", "higher_score_condition", "higher_score_against",
                     "equal_score"]
# Arguments: timing job directory, job id, then condition:against pairs. Compares the fingerprints of one job
# read by read; different_start counts reads aligned by both with a different first node or strand, and the
# score columns split those reads by which condition scored higher. Writes placement.tsv in the job directory.
PLACEMENT_PY = r'''import sys
t, job, pairs = sys.argv[1], sys.argv[2], [p.split(":") for p in sys.argv[3:]]
cache = {}
def load(name):
    if name not in cache:
        cache[name] = {}
        with open(f"{t}/round1_{name}/fingerprint.tsv") as f:
            next(f)
            for line in f:
                read, aligned, mapq, score, identity, start = line.split("\t")[:6]
                cache[name][read] = (aligned == "1", start, int(score))
    return cache[name]
with open(f"{t}/placement.tsv", "w") as out:
    out.write("COLUMNS\n")
    for cond, against in pairs:
        a, b = load(cond), load(against)
        if a.keys() != b.keys():
            sys.exit(f"{cond} and {against} have different read names")
        n = dict.fromkeys(["both", "only_a", "only_b", "diff", "higher_a", "higher_b", "equal"], 0)
        for read, (aligned_a, start_a, score_a) in a.items():
            aligned_b, start_b, score_b = b[read]
            if aligned_a and aligned_b:
                n["both"] += 1
                if start_a != start_b:
                    n["diff"] += 1
                    n["higher_a" if score_a > score_b else "higher_b" if score_b > score_a else "equal"] += 1
            elif aligned_a or aligned_b:
                n["only_a" if aligned_a else "only_b"] += 1
        out.write("\t".join(map(str, [job, cond, against, len(a), *n.values()])) + "\n")'''.replace(
    "COLUMNS", "\\t".join(PLACEMENT_COLUMNS))


def placement_pairs(conditions):
    """Within each non-reference density group, every other method against the BASELINE_LABEL condition."""
    base = {c["density_group"]: c["name"] for c in conditions if c["label"] == BASELINE_LABEL}
    return [(c["name"], base[c["density_group"]]) for c in conditions
            if c["density_group"] != REFERENCE_GROUP and c["density_group"] in base and c["label"] != BASELINE_LABEL]


def timing_script(conditions, reads, real=False, variant="", exclusive=True, rounds=None, offset=0, sampled=False):
    """One exclusive job that maps the conditions TIMING_ROUNDS times over, reloading the indexes every run.
    Round r starts offset + r - 1 places further along the given order, so no condition is always first or last.
    Only Giraffe runs; each GAM is deleted once written. Speeds are compared within a job.
    real: full real-read FASTQ, REAL_THREADS threads, one round in the given order; before each GAM is deleted,
    a per-read fingerprint and MAPQ/identity statistics are written, and the job ends with placement.tsv.
    sampled (with real): a sampled real-read FASTQ, mapped with the settings of simulated reads (THREADS
    threads, the rounds, map memory) and with real's fingerprints, statistics and placement.tsv.
    exclusive=False: no --exclusive and no --exclude, for many repeated jobs on shared nodes; their runs go to
    results/<reads>/timing_repeats_runs/ instead of timing/, which collect reads.
    rounds: rounds per job instead of TIMING_ROUNDS (full real reads: 1)."""
    reads_set = reads_set_name(reads)
    # Shared-node runs go apart from the exclusive ones, so collect's timing tables never mix the two
    t = DATA_DIR / "results" / reads_set / ("timing" if exclusive else "timing_repeats_runs") / "job$SLURM_JOB_ID"
    log = DATA_DIR / "logs" / f"timing_{reads_set}{variant}.%j.log"
    pinned = f"#SBATCH --exclusive\n#SBATCH --exclude={EXCLUDE_NODES}\n" if exclusive else ""
    columns, rounds, threads = TIMING_COLUMNS + (REAL_COLUMNS if real else []), rounds or TIMING_ROUNDS, THREADS
    if real and not sampled:
        rounds, threads = 1, REAL_THREADS
        text = sbatch_header(f"timing-{reads_set}{variant}", REAL_MEM, REAL_TIME, log, REAL_PARTITION, threads, pinned)
    else:
        mem = f"{max(int(map_mem(c).rstrip('G')) for c in conditions)}G"
        text = sbatch_header(f"timing-{reads_set}{variant}", mem, TIMING_TIME, log, extra=pinned)
    header = "\\t".join(columns)
    formats = "\\t".join(["%s"] * (len(TIMING_COLUMNS) + (2 if real else 0)))
    real_fields = ("" if not real else
                   " \\\n        \"$(awk '/Maximum resident/ {printf \"%.2f\", $NF / 1048576}' $D/giraffe.time)\""
                   " \"$(tail -n 1 $D/real_stats.tsv)\"")
    # A full real-read GAM is ~166 GB: delete it on any exit, including errors and Slurm's SIGTERM; logs stay.
    cleanup = "" if not real else """D=
trap '[ -z "$D" ] || rm -f "$D/mapped.gam"' EXIT
trap 'exit 143' TERM
"""
    text += "".join(require_index(c) for c in conditions) + f"""T={t}
mkdir -p $T
{cleanup}NODE=$(hostname -s)
printf 'git_commit\\trun_date\\thost\\tjob_id\\n%s\\t%s\\t%s\\t%s\\n' \\
    '{git_commit()}' "$(date -Iseconds)" "$NODE" "$SLURM_JOB_ID" > $T/run_info.tsv
printf '{header}\\n' > $T/runs.tsv
# Append one row to runs.tsv from the giraffe.log in $D. Arguments: round, order, condition.
record() {{
    L=$D/giraffe.log
    printf '{formats}\\n' "$SLURM_JOB_ID" "$NODE" "$1" "$2" "$3" \\
        "$(sed -n 's/.*Loading and initialization: \\([0-9.]*\\) seconds.*/\\1/p' $L)" \\
        "$(sed -n 's/.*Mapped [0-9]* reads across [0-9]* threads in \\([0-9.]*\\) seconds.*/\\1/p' $L)" \\
        "$(sed -n 's/.*Mapping speed: \\([0-9.]*\\) reads per second per thread.*/\\1/p' $L)" \\
        "$(sed -n 's/.*Achieved \\([0-9.]*\\) reads per CPU-second.*/\\1/p' $L)" \\
        "$(sed -n 's/.*Used \\([0-9.]*\\) CPU-seconds.*/\\1/p' $L)" \\
        "$(sed -n 's/.*Memory footprint: \\([0-9.]*\\) GB.*/\\1/p' $L)"{real_fields} >> $T/runs.tsv
}}
"""
    order = 0
    for rnd in range(1, rounds + 1):
        shift = (offset + rnd - 1) % len(conditions)
        for c in conditions[shift:] + conditions[:shift]:
            order += 1
            stats = "" if not real else f"""{c['vg_binary']} filter -t {threads} -T "{REAL_TSV_FIELDS}" $D/mapped.gam \\
    | awk '{REAL_FINGERPRINT_AWK}' > $D/fingerprint.tsv
python3 -I -c '{REAL_STATS_PY}' \\
    $D/fingerprint.tsv $D "$(sed -n 's/.*Mapped \\([0-9]*\\) reads across.*/\\1/p' $D/giraffe.log)"
"""
            text += f"""
# round {rnd}, run {order}: {c['name']}
D=$T/round{rnd}_{c['name']}
mkdir -p $D
/usr/bin/time -v -o $D/giraffe.time {giraffe_cmd(c, reads, threads=threads, reads_flag="-f" if real else "-G")} \\
    > $D/mapped.gam 2> $D/giraffe.log
{stats}rm $D/mapped.gam
record {rnd} {order} {c['name']}
"""
    if real:
        pairs = " ".join(f"{a}:{b}" for a, b in placement_pairs(conditions))
        text += f"""
# Read-by-read placement comparison of this job's runs
python3 -I -c '{PLACEMENT_PY}' $T "$SLURM_JOB_ID" {pairs}
cat $T/placement.tsv
"""
    return text + "cat $T/runs.tsv\n"


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


def prepare_scripts():
    """One-time jobs per graph, independent of the reads: (name, sbatch text)."""
    log = DATA_DIR / "logs"
    paths_meta = DATA_DIR / "graph.paths_meta.tsv"
    length = DATA_DIR / "graph.length.txt"
    return [
        # List path metadata of the graph and summarize which sample is REFERENCE.
        ("check_paths", sbatch_header("check-paths", "32G", "00:30:00", log / "check_paths.%j.log",
                                      partition="short", cpus=2) + f"""/usr/bin/time -v {VG} paths -M -x {GBZ} > {paths_meta}

echo "=== path counts by SENSE, SAMPLE"
awk 'NR>1 {{print $2"\\t"$3}}' {paths_meta} | sort | uniq -c
echo "=== REFERENCE paths (first 40)"
awk 'NR==1 || $2=="REFERENCE"' {paths_meta} | head -41
"""),
        # Build the distance index shared by all seed conditions; record peak memory/time.
        ("build_dist", sbatch_header("build-dist", "128G", INDEX_TIME, log / "build_dist.%j.log")
         + f"""/usr/bin/time -v {VG} index -t "$SLURM_CPUS_PER_TASK" -j {DIST}.tmp {GBZ}
mv {DIST}.tmp {DIST}
echo "dist index size (bytes): $(stat -c %s {DIST})"
"""),
        # Total sequence length of the graph: denominator for index-based seed density.
        ("graph_length", sbatch_header("graph-length", "32G", "00:30:00", log / "graph_length.%j.log",
                                       partition="short", cpus=2) + f"""/usr/bin/time -v {VG} stats -l {GBZ} > {length}
cat {length}
"""),
    ]


# One line per alignment: name, MAPQ, score, identity, number of mappings, start (node:offset:strand) and
# a path signature (nodes, offsets, strands, edit lengths). Read bases never leave jq.
FINGERPRINT_JQ = (
    '(.path.mapping // []) as $m | [.name, (.mapping_quality // 0), (.score // 0), (.identity // 0), ($m | length), '
    '($m[0].position | "\\(.node_id // "NA"):\\(.offset // 0):\\(if .is_reverse then "-" else "+" end)"), '
    '($m | map("\\(.position.node_id):\\(.position.offset // 0):\\(if .position.is_reverse then "-" else "+" end)|" '
    '+ ((.edit // []) | map("\\(.from_length // 0),\\(.to_length // 0),\\((.sequence // "") | length)") | join(";"))) '
    '| join(" "))] | @tsv')
# Replace the path signature (last column) by a short hash.
FINGERPRINT_PY = r'''import sys, hashlib
for line in sys.stdin:
    *fields, signature = line.rstrip("\n").split("\t")
    print("\t".join(fields + [hashlib.sha1(signature.encode()).hexdigest()[:16]]))'''
FINGERPRINT_FIELDS = ["mapq", "score", "identity", "n_mappings", "start", "path_hash"]


def fingerprint_script(c, run_dir, label):
    """Job writing run_dir/fingerprint.tsv from run_dir/mapped.gam (Giraffe's own output)."""
    gam, out = run_dir / "mapped.gam", run_dir / "fingerprint.tsv"
    log = DATA_DIR / "logs" / f"fingerprint_{label}.%j.log"
    header = "\\t".join(["read"] + FINGERPRINT_FIELDS)
    return sbatch_header(f"fp-{label}", "8G", MAP_TIME, log, cpus=2) + f"""for tool in jq python3; do
    command -v $tool > /dev/null || {{ echo "$tool not found on $(hostname)" >&2; exit 1; }}
done
[ -s {gam} ] || {{ echo "missing {gam}" >&2; exit 1; }}
{{ printf '{header}\\n'
  {c['vg_binary']} view -aj {gam} | jq -r '{FINGERPRINT_JQ}' | python3 -I -c '{FINGERPRINT_PY}'; }} > {out}.tmp
rows=$(($(wc -l < {out}.tmp) - 1))
unique=$(tail -n +2 {out}.tmp | cut -f1 | sort -u | wc -l)
alignments=$(grep -oP 'Total alignments: \\K\\d+' {run_dir}/stats.txt)
echo "rows=$rows unique_names=$unique total_alignments=$alignments"
[ "$rows" = "$unique" ] || {{ echo "duplicate read names in the fingerprint" >&2; exit 1; }}
[ "$rows" = "$alignments" ] || {{ echo "fingerprint rows differ from vg stats -a" >&2; exit 1; }}
mv {out}.tmp {out}
"""


def fingerprint_runs(c, reads, run_tag):
    """Run directories (with a mapped.gam) of one condition: <name>/<run_tag>/rep* if run_tag, else <name>/rep*."""
    base = results_dir(c, reads) / run_tag if run_tag else results_dir(c, reads)
    runs = sorted((p for p in base.glob("rep*") if re.fullmatch(r"rep\d+", p.name)), key=lambda p: int(p.name[3:]))
    return [p for p in runs if (p / "mapped.gam").is_file()]


def read_compare(path):
    """read name -> (eligible, correct, mapq) from a gamcompare -T table."""
    with open(path, newline="") as f:
        return {r["read"]: (r["eligible"] == "1", r["correct"] == "1", int(r["mq"]))
                for r in csv.DictReader(f, delimiter="\t")}


def write_table(path, columns, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([format_value(row.get(col)) for col in columns])
    print(f"wrote {path} ({len(rows)} rows)")


def read_fingerprints(path):
    """read name -> fingerprint tuple; keyed by name because threaded Giraffe writes reads in any order."""
    with open(path, newline="") as f:
        return {r["read"]: tuple(r[k] for k in FINGERPRINT_FIELDS) for r in csv.DictReader(f, delimiter="\t")}


def classify(prints):
    """How a read's fingerprints differ across runs, most severe first."""
    if len(set(prints)) == 1:
        return "identical"
    if len({p[3] == "0" for p in prints}) > 1:
        return "mapped_flip"
    if len({p[4] for p in prints}) > 1:
        return "start_diff"
    if len({p[5] for p in prints}) > 1:
        return "path_diff_same_start"
    return "score_only"


DETERMINISM_CATEGORIES = ["identical", "score_only", "path_diff_same_start", "start_diff", "mapped_flip"]


def determinism(conditions, reads, run_tag):
    """Compare per-read fingerprints across the run_tag replicates, and separately with the earlier rep1."""
    out_dir = DATA_DIR / "results" / Path(reads).stem / "determinism"
    summary, details = [], []
    for c in conditions:
        new = [p for p in fingerprint_runs(c, reads, run_tag) if (p / "fingerprint.tsv").is_file()]
        old = [p for p in fingerprint_runs(c, reads, None)[:1] if (p / "fingerprint.tsv").is_file()]
        if len(new) < 2:
            print(f"[{c['name']}] skipped: fewer than two {run_tag} runs with a fingerprint")
            continue
        comparisons = [(f"{run_tag} runs", new)]
        if old:
            comparisons.append((f"{run_tag} runs + earlier rep1", new + old))
        for comparison, runs in comparisons:
            labels = [str(p.relative_to(results_dir(c, reads))) for p in runs]
            prints = [read_fingerprints(p / "fingerprint.tsv") for p in runs]
            if len({frozenset(p) for p in prints}) != 1:
                sys.exit(f"[{c['name']}] {comparison}: runs do not contain the same reads")
            correct = [{read: (e, k) for read, (e, k, mq) in read_compare(p / "compare.tsv").items()} for p in runs]
            nodes = [(read_density_summary(p / "run_info.tsv") or {}).get("host", "NA") for p in runs]
            counts = dict.fromkeys(DETERMINISM_CATEGORIES, 0)
            correct_flips = 0
            for read in sorted(prints[0]):
                category = classify([p[read] for p in prints])
                counts[category] += 1
                status = [k.get(read) for k in correct]
                flip = len({s for s in status if s is not None}) > 1
                correct_flips += flip
                if category != "identical" or flip:
                    for label, node, p, s in zip(labels, nodes, prints, status):
                        details.append({"comparison": comparison, "condition": c["name"], "read": read,
                                        "category": category, "correct_flip": "yes" if flip else "no",
                                        "run": label, "node": node, **dict(zip(FINGERPRINT_FIELDS, p[read])),
                                        "correct": "NA" if s is None else ("yes" if s[1] else "no")})
            summary.append({"condition": c["name"], "comparison": comparison, "runs": ",".join(labels),
                            "nodes": ",".join(nodes), "reads": len(prints[0]), **counts,
                            "correct_flip": correct_flips})
    if not summary:
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    write_table(out_dir / "summary.tsv", ["condition", "comparison", "runs", "nodes", "reads",
                                          *DETERMINISM_CATEGORIES, "correct_flip"], summary)
    write_table(out_dir / "reads.tsv", ["comparison", "condition", "read", "category", "correct_flip", "run", "node",
                                        *FINGERPRINT_FIELDS, "correct"], details)


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
    "timed_jobs", "timed_speed_median", "timed_ratio_min", "timed_ratio_max",
    "seeds_per_read_p99", "seeds_per_read_max",
]


SNAPSHOT_FILES = ["results.tsv", "mcnemar.tsv", "comparison.png", "density_sweep.png",
                  "filter_summary.tsv", "minimizer_summary.tsv", "timing.tsv", "timing_runs.tsv", "placement.tsv"]


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


TIMING_SUMMARY_COLUMNS = ["condition", "job_id", "node", "runs", "speed_median", "speed_min", "speed_max",
                          "ratio_vs_ref", "cpu_speed_median", "load_s_median", "map_s_median", "memory_gb_max"]


def timing_summary(conditions, reads_set):
    """One row per condition and timing job: medians over that job's runs, and ratio_vs_ref = this speed median
    over the reference condition's speed median in the same job (NA if the reference did not run there)."""
    names = [c["name"] for c in conditions]
    reference = next((c["name"] for c in conditions if c["density_group"] == REFERENCE_GROUP), None)
    runs = collections.defaultdict(list)
    for path in sorted((DATA_DIR / "results" / reads_set / "timing").glob("job*/runs.tsv")):
        with open(path, newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["condition"] in names and r["reads_per_sec_per_thread"]:  # empty: run did not finish
                    runs[r["job_id"], r["condition"]].append(r)
    medians = {key: statistics.median(float(r["reads_per_sec_per_thread"]) for r in rs) for key, rs in runs.items()}
    rows = []
    for job in sorted({job for job, _ in runs}, key=int):
        for name in names:
            rs = runs.get((job, name))
            if not rs:
                continue
            speeds = [float(r["reads_per_sec_per_thread"]) for r in rs]
            ref = medians.get((job, reference))
            rows.append({"condition": name, "job_id": job, "node": rs[0]["node"], "runs": len(rs),
                         "speed_median": medians[job, name], "speed_min": min(speeds), "speed_max": max(speeds),
                         "ratio_vs_ref": medians[job, name] / ref if ref else None,
                         "cpu_speed_median": statistics.median(float(r["reads_per_cpu_second"]) for r in rs),
                         "load_s_median": statistics.median(float(r["load_seconds"]) for r in rs),
                         "map_s_median": statistics.median(float(r["map_seconds"]) for r in rs),
                         "memory_gb_max": max(float(r["memory_gb"]) for r in rs)})
    return rows


TIMING_TABLES = {"runs.tsv": "timing_runs.tsv", "placement.tsv": "placement.tsv"}  # per job -> per read set


def concat_tables(paths, out):
    """Concatenate TSV tables with the same header into out; remove a stale out when there are none."""
    if not paths:
        out.unlink(missing_ok=True)
        return
    header, rows = None, []
    for path in paths:
        lines = path.read_text().splitlines()
        if header not in (None, lines[0]):
            sys.exit(f"{path}: header differs from {paths[0]}")
        header = lines[0]
        rows += lines[1:]
    out.write_text("\n".join([header] + rows) + "\n")
    print(f"wrote {out} ({len(rows)} rows)")


def seed_stats(per_read):
    """Median, mean, nearest-rank 90th and 99th percentiles and maximum of seeds per read (stage.seed.results)."""
    with open(per_read, newline="") as f:
        seeds = sorted(int(r["seeds"]) for r in csv.DictReader(f, delimiter="\t") if r["seeds"] != "NA")
    if not seeds:
        return {}
    return {"seeds_per_read_median": statistics.median(seeds), "seeds_per_read_mean": statistics.mean(seeds),
            "seeds_per_read_p90": seeds[math.ceil(0.9 * len(seeds)) - 1],
            "seeds_per_read_p99": seeds[math.ceil(0.99 * len(seeds)) - 1], "seeds_per_read_max": seeds[-1]}


def format_value(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


REFERENCE_GROUP = "reference"  # density_group of the condition every other one is tested against
BASELINE_LABEL = "minimizer"  # within a density group, the other methods are tested against this label
# Methods placed in a density group by their measured density rather than the w/s formula; they are also
# tested against every other method of their group.
DENSITY_MATCHED_LABELS = ("minimaxmer",)
MCNEMAR_COLUMNS = ["family", "condition", "against", "eligible_both", "wrong_condition", "wrong_against",
                   "only_condition_correct", "only_against_correct", "net_reads", "net_percent",
                   "p_exact", "m_tests", "p_bonferroni"]


def mcnemar_exact(b, c):
    """Two-sided exact McNemar p-value: if both conditions are equally accurate, b ~ Binomial(b + c, 1/2)."""
    n = b + c
    if n == 0:
        return 1.0
    tail = sum(math.comb(n, i) for i in range(min(b, c) + 1)) / 2 ** n
    return min(1.0, 2 * tail)


def primary_compare(c, reads_set):
    """compare.tsv of the run accuracy is read from: rep1 when the condition has replicates, else the single run."""
    o = DATA_DIR / "results" / reads_set / c["name"]
    return o / "rep1" / "compare.tsv" if (o / "rep1" / "compare.tsv").is_file() else o / "compare.tsv"


def mcnemar_table(rows, reads_set):
    """Paired exact McNemar tests on per-read correctness: every condition against the reference group,
    every other method against BASELINE_LABEL within each density group, and a DENSITY_MATCHED_LABELS method
    also against each other method of its group. Bonferroni over all tests of this read set.
    net_reads = only_condition_correct - only_against_correct (positive: condition is better)."""
    pairs = []
    reference = [r for r in rows if r["density_group"] == REFERENCE_GROUP]
    if reference:
        pairs += [("vs_reference", r, reference[0]) for r in rows if r["density_group"] != REFERENCE_GROUP]
    for group in dict.fromkeys(r["density_group"] for r in rows if r["density_group"] != REFERENCE_GROUP):
        members = [r for r in rows if r["density_group"] == group]
        for base in [r for r in members if r["label"] == BASELINE_LABEL][:1]:
            pairs += [(f"within_{group}", r, base) for r in members if r is not base]
        for matched in [r for r in members if r["label"] in DENSITY_MATCHED_LABELS]:
            pairs += [(f"within_{group}", matched, r) for r in members
                      if r is not matched and r["label"] not in (BASELINE_LABEL, *DENSITY_MATCHED_LABELS)]
    per_read = {}
    for r in rows:
        per_read[r["name"]] = read_compare(primary_compare(r, reads_set))
        wrong = sum(1 for eligible, correct, mq in per_read[r["name"]].values() if eligible and not correct)
        if wrong != r["wrong"]:
            sys.exit(f"{r['name']}: compare.tsv has {wrong} wrong reads but results.tsv has {r['wrong']}")
    tests = []
    for family, a, b in pairs:
        A, B = per_read[a["name"]], per_read[b["name"]]
        shared = [read for read, (eligible, correct, mq) in A.items() if eligible and read in B and B[read][0]]
        only_a = sum(1 for read in shared if A[read][1] and not B[read][1])
        only_b = sum(1 for read in shared if B[read][1] and not A[read][1])
        tests.append({"family": family, "condition": a["name"], "against": b["name"], "eligible_both": len(shared),
                      "wrong_condition": sum(1 for read in shared if not A[read][1]),
                      "wrong_against": sum(1 for read in shared if not B[read][1]),
                      "only_condition_correct": only_a, "only_against_correct": only_b,
                      "net_reads": only_a - only_b, "net_percent": 100 * (only_a - only_b) / len(shared),
                      "p_exact": mcnemar_exact(only_a, only_b)})
    for t in tests:
        t["m_tests"] = len(tests)
        t["p_bonferroni"] = min(1.0, t["p_exact"] * len(tests))
    return tests


def signed(n):
    """+3, −5 (true minus sign) or ±0 for plot labels."""
    return "±0" if n == 0 else (f"+{n}" if n > 0 else f"−{-n}")


METHOD_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # categorical slots 1-4 of the reference palette
REFERENCE_COLOR = "#8a8984"
MAX_LABELLED_CONDITIONS = 6  # above this, bars carry no value labels


def method_colors(rows):
    """One color per method (label) in order of first appearance; the reference group is gray."""
    labels = list(dict.fromkeys(r["label"] for r in rows if r["density_group"] != REFERENCE_GROUP))
    if len(labels) > len(METHOD_COLORS):
        sys.exit(f"plots support at most {len(METHOD_COLORS)} methods, got {labels}")
    return dict(zip(labels, METHOD_COLORS))


def row_color(r, colors):
    return REFERENCE_COLOR if r["density_group"] == REFERENCE_GROUP else colors[r["label"]]


def nominal_density(r):
    """Measured density on reads if known, else the density group read as a fraction (e.g. 2/21)."""
    if r.get("density_reads") is not None:
        return r["density_reads"]
    match = re.fullmatch(r"(\d+)/(\d+)", r["density_group"])
    return int(match.group(1)) / int(match.group(2)) if match else 0.0


def style_axes(ax):
    ax.set_facecolor("#fcfcfb")
    ax.grid(color="#e6e5e0", linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#c3c2b7")
    ax.tick_params(colors="#52514e")


def plot_order(rows, colors):
    """Grouped by method, densest first within a method, reference conditions last."""
    method_rank = {label: i for i, label in enumerate(colors)}
    return sorted(rows, key=lambda r: (r["density_group"] == REFERENCE_GROUP,
                                       method_rank.get(r["label"], len(method_rank)), -nominal_density(r)))


def plot_comparison(rows, path, reads_set, tests=()):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    ink, muted = "#0b0b0b", "#52514e"
    rows = [r for r in rows if r.get("accuracy") is not None]
    method_color = method_colors(rows)
    rows = plot_order(rows, method_color)
    names = [r["name"] for r in rows]
    colors = [row_color(r, method_color) for r in rows]
    x = range(len(rows))
    labelled = len(rows) <= MAX_LABELLED_CONDITIONS

    # The seeds panels need this read set's density runs; leave them out if none exist. No speed panel:
    # speeds are compared with the repeated timing runs (prototype/timing_repeats.py), not single map runs.
    has_seeds = any(r.get("seeds_per_read_median") is not None for r in rows)
    panels = 3 if has_seeds else 1
    # At least 12 inches wide, so the legend (top right) clears the read set name (top left) with one panel.
    fig, axes = plt.subplots(1, panels, figsize=(max(12, panels * max(5, 0.55 * len(rows) + 1.5)), 5.2), squeeze=False)
    axes = axes[0]
    ax_err = axes[0]
    error = [100 * r["wrong"] / r["eligible"] for r in rows]
    ax_err.bar(x, error, color=colors, width=0.6)
    # No error bars: conditions map the same reads, so differences are judged by the paired test.
    vs_reference = {t["condition"]: t for t in tests if t["family"] == "vs_reference"}
    for i, r in zip(x, rows):
        if labelled:
            label = f"{r['wrong']}/{r['eligible']}"
            t = vs_reference.get(r["name"])
            if t:
                label += f"\n{signed(t['wrong_condition'] - t['wrong_against'])} wrong vs ref\n(p={t['p_bonferroni']:.2g})"
            ax_err.annotate(label, (i, error[i]), xytext=(0, 4), textcoords="offset points",
                            ha="center", va="bottom", fontsize=8, color=ink)
    title = "Mismapped reads\n(% of reads with a truth position on CHM13)"
    if vs_reference and labelled:
        reference = next(iter(vs_reference.values()))
        title += (f"\nref = {reference['against']}; p: paired exact McNemar,\n"
                  f"Bonferroni over {reference['m_tests']} tests")
    elif tests:
        title += "\npaired tests: see mcnemar.tsv"
    ax_err.set_title(title, fontsize=10, color=ink, loc="left")
    ax_err.set_ylim(0, max(error) * (1.6 if labelled else 1.15))

    if has_seeds:
        # Solid bar: median seeds per read; light bar: mean (pulled up by a few repetitive reads).
        ax_seeds = axes[1]
        width = 0.38
        for i, r in zip(x, rows):
            for offset, stat, alpha in ((-width / 2, "median", 1.0), (width / 2, "mean", 0.4)):
                value = r.get(f"seeds_per_read_{stat}")
                if value is None:
                    continue
                ax_seeds.bar(i + offset, value, width=width, color=colors[i], alpha=alpha,
                             edgecolor="#fcfcfb", linewidth=2)
                if labelled:
                    ax_seeds.annotate(f"{value:.0f}", (i + offset, value), xytext=(0, 3), textcoords="offset points",
                                      ha="center", fontsize=7, color=ink)
        top = max(r.get("seeds_per_read_mean") or 0 for r in rows)
        ax_seeds.set_ylim(0, 1.15 * top if top > 0 else 1)
        ax_seeds.set_title(f"Seed hits per read (graph positions)\non {reads_set} (solid: median, light: mean)",
                           fontsize=11, color=ink, loc="left")

        # The high end gets its own panel: p99 and the maximum are one to two orders of magnitude above the
        # median and would flatten the panel above. Bar: p99; number above it: the maximum (one read, noisy).
        ax_high = axes[2]
        for i, r in zip(x, rows):
            p99, top_read = r.get("seeds_per_read_p99"), r.get("seeds_per_read_max")
            if p99 is None:
                continue
            ax_high.bar(i, p99, color=colors[i], width=0.6)
            ax_high.annotate(f"{top_read}", (i, p99), xytext=(0, 3), textcoords="offset points", ha="center",
                             va="bottom", fontsize=7 if labelled else 6, rotation=0 if labelled else 90, color=ink)
        top = max((r.get("seeds_per_read_p99") or 0) for r in rows)
        ax_high.set_ylim(0, (1.2 if labelled else 1.6) * top if top > 0 else 1)
        ax_high.set_title(f"Seed hits per read, high end\non {reads_set} (bar: p99, number: maximum)",
                          fontsize=11, color=ink, loc="left")

    for ax in axes:
        ax.set_xticks(list(x))
        ax.set_xticklabels(names, rotation=45, ha="right", rotation_mode="anchor", fontsize=9, color=ink)
        style_axes(ax)
        ax.grid(axis="x", visible=False)
    legend = [Patch(color=color, label=label) for label, color in method_color.items()]
    legend += [Patch(color=REFERENCE_COLOR, label=r["label"]) for r in rows if r["density_group"] == REFERENCE_GROUP]
    fig.legend(handles=legend, loc="upper right", ncol=len(legend), frameon=False, fontsize=9, labelcolor=ink)
    fig.suptitle(reads_set, fontsize=10, color=muted, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor="#fcfcfb")
    plt.close(fig)


def plot_density_sweep(rows, path, reads_set):
    """Measured seed density on reads vs error rate, seeds per read and giraffe memory; one line per
    method, reference conditions as gray points. Returns False (no plot) with fewer than two density groups."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, muted = "#0b0b0b", "#52514e"
    rows = [r for r in rows if r.get("accuracy") is not None and r.get("density_reads") is not None]
    groups = {r["density_group"] for r in rows if r["density_group"] != REFERENCE_GROUP}
    if len(groups) < 2:
        return False
    method_color = method_colors(rows)
    series = [(label, color, sorted((r for r in rows if r["label"] == label and r["density_group"] != REFERENCE_GROUP),
                                    key=lambda r: r["density_reads"]))
              for label, color in method_color.items()]
    series += [(r["label"], REFERENCE_COLOR, [r]) for r in rows if r["density_group"] == REFERENCE_GROUP]
    point = dict(marker="o", ms=7, markeredgecolor="#fcfcfb", markeredgewidth=1.5)

    fig, (ax_err, ax_seeds, ax_high, ax_mem) = plt.subplots(1, 4, figsize=(21, 5.2))
    top_error = 0.0
    for label, color, members in series:
        error = [100 * r["wrong"] / r["eligible"] for r in members]
        top_error = max([top_error] + error)
        ax_err.plot([r["density_reads"] for r in members], error, color=color, lw=2, label=label, **point)
        seeds = [m for m in members if m.get("seeds_per_read_mean") is not None]
        ax_seeds.plot([m["density_reads"] for m in seeds], [m["seeds_per_read_mean"] for m in seeds],
                      color=color, lw=2, **point)
        ax_seeds.plot([m["density_reads"] for m in seeds], [m["seeds_per_read_median"] for m in seeds],
                      color=color, lw=1.5, ls="--", **point)
        high = [m for m in members if m.get("seeds_per_read_p99") is not None]
        ax_high.plot([m["density_reads"] for m in high], [m["seeds_per_read_p99"] for m in high], color=color, lw=2, **point)
        mem = [m for m in members if m.get("giraffe_mem_gb") is not None]
        ax_mem.plot([m["density_reads"] for m in mem], [m["giraffe_mem_gb"] for m in mem], color=color, lw=2, **point)

    titles = ["Mismapped reads\n(% of reads with a truth position on CHM13)",
              f"Seed hits per read (graph positions)\non {reads_set} (solid: mean, dashed: median)",
              f"Seed hits per read, p99\non {reads_set} (own panel: far above the median)",
              "Giraffe memory (GB)"]
    for ax, title in zip((ax_err, ax_seeds, ax_high, ax_mem), titles):
        ax.set_title(title, fontsize=11, color=ink, loc="left")
        ax.set_xlabel("Measured seed density on reads (seeds / bp)", color=muted)
        ax.set_ylim(bottom=0)
        style_axes(ax)
    ax_err.set_ylim(0, 1.15 * top_error)  # keep the highest point (often the reference) off the frame
    handles, labels = ax_err.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=len(labels), frameon=False, fontsize=9, labelcolor=ink)
    fig.suptitle(reads_set, fontsize=10, color=muted, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor="#fcfcfb")
    plt.close(fig)
    return True


def collect(conditions, reads, snapshot=False):
    reads_set = reads_set_name(reads)
    graph_length = grep_first(DATA_DIR / "graph.length.txt", r"length\s+(\d+)", int)
    rows, skipped = [], []
    for c in conditions:
        row = collect_row(c, reads_set, graph_length)
        (rows if row.get("correct") is not None else skipped).append(row)
    out_dir = DATA_DIR / "results" / reads_set
    out_dir.mkdir(parents=True, exist_ok=True)
    # Fixed-hardware speed (timing stage): within each job, then summarised over jobs.
    timing = timing_summary(conditions, reads_set)
    for row in rows:
        jobs = [t for t in timing if t["condition"] == row["name"]]
        if jobs:
            ratios = [t["ratio_vs_ref"] for t in jobs if t["ratio_vs_ref"] is not None]
            row["timed_jobs"] = len(jobs)
            row["timed_speed_median"] = statistics.median(t["speed_median"] for t in jobs)
            if ratios:
                row["timed_ratio_min"], row["timed_ratio_max"] = min(ratios), max(ratios)
    if timing:
        write_table(out_dir / "timing.tsv", TIMING_SUMMARY_COLUMNS, timing)
    else:
        (out_dir / "timing.tsv").unlink(missing_ok=True)  # never snapshot a stale timing table
    # Every timing run, and (real reads) the read-by-read placement comparisons, from all timing jobs.
    for name in ("runs.tsv", "placement.tsv"):
        concat_tables(sorted((out_dir / "timing").glob(f"job*/{name}")), out_dir / TIMING_TABLES[name])
    table = out_dir / "results.tsv"
    if not rows:   # e.g. a real-read set with timing runs only: no empty tables to snapshot
        for name in ("results.tsv", "mcnemar.tsv"):
            (out_dir / name).unlink(missing_ok=True)
        print(f"no map results for {reads_set}: results.tsv and mcnemar.tsv not written (timing tables only)")
    else:
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
    tests = mcnemar_table(rows, reads_set) if rows else []
    if rows:
        write_table(out_dir / "mcnemar.tsv", MCNEMAR_COLUMNS, tests)
        plot = out_dir / "comparison.png"
        plot_comparison(rows, plot, reads_set, tests)
        print(f"wrote {plot}")
        sweep = out_dir / "density_sweep.png"
        if plot_density_sweep(rows, sweep, reads_set):
            print(f"wrote {sweep}")
        else:
            sweep.unlink(missing_ok=True)  # never snapshot a stale sweep from an earlier run
            print(f"density_sweep.png not drawn for {reads_set}: needs results at two or more density groups")
    if snapshot:
        dest = REPO_DIR / "results_summary" / reads_set
        dest.mkdir(parents=True, exist_ok=True)
        for name in SNAPSHOT_FILES:
            if (out_dir / name).is_file():
                shutil.copy2(out_dir / name, dest / name)
                print(f"snapshot {dest / name}")


def submit(script, args):
    """sbatch a generated script; --fixed-nodes and --sbatch-args only change the sbatch command line."""
    cmd = ["sbatch"]
    if args.fixed_nodes:
        if not EXCLUDE_NODES:
            sys.exit("--fixed-nodes needs exclude_nodes in the config")
        cmd.append(f"--exclude={EXCLUDE_NODES}")
    cmd += shlex.split(args.sbatch_args or "")
    result = subprocess.run(cmd + [str(script)], check=True, capture_output=True, text=True)
    print("  " + result.stdout.strip())


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("stage", choices=sorted(STAGES) + ["prepare", "filters", "collect", "fingerprint", "determinism",
                                                           "timing", "checksums"])
    parser.add_argument("--config", default=str(DEFAULT_CONFIG),
                        help="tab-separated paths: data_dir, graph, reads_dir, vg (default: config.tsv in this repo)")
    parser.add_argument("--conditions", default="conditions.tsv")
    parser.add_argument("--reads", help=f"truth GAM for map/density (default: reads_dir/{DEFAULT_READS_NAME})")
    parser.add_argument("--only", nargs="*", help="run only these condition names")
    parser.add_argument("--submit", action="store_true", help="sbatch the generated scripts")
    parser.add_argument("--force", action="store_true",
                        help="index: rebuild even if already built; checksums: rewrite existing records")
    parser.add_argument("--reps", type=int, default=1, help="map: number of replicate runs (default 1)")
    parser.add_argument("--run-tag", help="map/fingerprint/determinism: runs live in <name>/<run-tag>/ "
                                          "(apart from the default runs)")
    parser.add_argument("--keep-gams", action="store_true", help="map: keep the GAMs of every replicate")
    parser.add_argument("--fixed-nodes", action="store_true",
                        help="submit with sbatch --exclude=<exclude_nodes from the config>")
    parser.add_argument("--sbatch-args", help="extra arguments added to the sbatch command line on --submit")
    parser.add_argument("--real", action="store_true",
                        help=f"timing: --reads is a real-read FASTQ; {REAL_THREADS} threads, one round per job, "
                             "fingerprints and MAPQ/identity statistics, two jobs in opposite orders")
    parser.add_argument("--sampled", action="store_true",
                        help=f"timing --real: the FASTQ is a sample (e.g. 100k reads); {THREADS} threads and "
                             f"{TIMING_ROUNDS} rotated rounds as for simulated reads, keeping --real's statistics")
    parser.add_argument("--no-exclusive", action="store_true",
                        help="timing: no --exclusive and no --exclude (shared nodes; for --jobs)")
    parser.add_argument("--jobs", type=int,
                        help="timing (not full --real): write N scripts, job k starting its rounds k-1 places "
                             "further along the order; submit them by hand (not with --submit)")
    parser.add_argument("--rounds", type=int, help=f"timing (not full --real): rounds per job (default {TIMING_ROUNDS})")
    parser.add_argument("--snapshot", action="store_true",
                        help="collect: also copy the summaries to results_summary/<reads>/ in this repo")
    args = parser.parse_args()
    if args.reps < 1 or (args.reps > 1 and args.stage != "map"):
        sys.exit("--reps must be >= 1 and only applies to the map stage")
    if args.keep_gams and args.stage != "map":
        sys.exit("--keep-gams only applies to the map stage")
    if args.real and args.stage != "timing":
        sys.exit("--real only applies to the timing stage")
    if (args.sampled or args.no_exclusive or args.jobs or args.rounds) and args.stage != "timing":
        sys.exit("--sampled, --no-exclusive, --jobs and --rounds only apply to the timing stage")
    if args.sampled and not args.real:
        sys.exit("--sampled needs --real")
    if (args.jobs or args.rounds) and args.real and not args.sampled:
        sys.exit("--jobs and --rounds do not apply to full --real runs (one round, two jobs in opposite orders)")
    if (args.jobs is not None and args.jobs < 1) or (args.rounds is not None and args.rounds < 1):
        sys.exit("--jobs and --rounds must be >= 1")
    if args.jobs and args.submit:
        sys.exit("--jobs: submit the scripts by hand, e.g. in --dependency chains, to limit how many run at once")
    if args.run_tag and (args.stage not in ("map", "fingerprint", "determinism")
                         or not re.fullmatch(r"[A-Za-z0-9_-]+", args.run_tag) or re.fullmatch(r"rep\d+", args.run_tag)):
        sys.exit("--run-tag applies to map/fingerprint/determinism and must be a simple name other than rep<N>")
    if args.stage == "determinism" and not args.run_tag:
        sys.exit("determinism needs --run-tag (the replicate runs to compare)")
    load_config(args.config)
    args.reads = args.reads or str(DEFAULT_READS)

    if args.stage == "prepare":
        for directory in (DATA_DIR / "jobs", DATA_DIR / "logs"):
            directory.mkdir(parents=True, exist_ok=True)
        for name, text in prepare_scripts():
            script = DATA_DIR / "jobs" / f"prepare_{name}.sbatch"
            script.write_text(text)
            print(f"[{name}] prepare -> {script}")
            if args.submit:
                submit(script, args)
        return

    conditions = load_conditions(args.conditions)
    if args.only:
        unknown = set(args.only) - {c["name"] for c in conditions}
        if unknown:
            sys.exit(f"unknown condition names: {sorted(unknown)}")
        conditions = [c for c in conditions if c["name"] in args.only]

    if args.stage == "collect":
        collect(conditions, args.reads, args.snapshot)
        return
    if args.stage == "checksums":
        record_checksums(conditions, args.force)
        return
    if args.stage == "filters":
        for directory in (DATA_DIR / "jobs", DATA_DIR / "logs"):
            directory.mkdir(parents=True, exist_ok=True)
        script = DATA_DIR / "jobs" / f"filters_{Path(args.reads).stem}.sbatch"
        script.write_text(filters_script(conditions, args.reads))
        print(f"[{len(conditions)} conditions] filters -> {script}")
        if args.submit:
            submit(script, args)
        return
    if args.stage == "determinism":
        determinism(conditions, args.reads, args.run_tag)
        return
    if args.stage == "timing":
        if not args.only or not (EXCLUDE_NODES or args.no_exclusive):
            sys.exit("timing needs --only (the conditions, in running order) and, unless --no-exclusive, "
                     "exclude_nodes in the config")
        by_name = {c["name"]: c for c in conditions}
        conditions = [by_name[name] for name in args.only]
        for c in conditions:
            check_index_params(c, "map", False)  # exits if the index is missing or stale
        for directory in (DATA_DIR / "jobs", DATA_DIR / "logs"):
            directory.mkdir(parents=True, exist_ok=True)
        if not args.real or args.sampled:
            rounds = args.rounds or TIMING_ROUNDS
            for k in range(args.jobs or 1):
                variant = f"_job{k + 1}" if args.jobs else ""
                script = DATA_DIR / "jobs" / f"timing_{reads_set_name(args.reads)}{variant}.sbatch"
                script.write_text(timing_script(conditions, args.reads, real=args.real, variant=variant,
                                                exclusive=not args.no_exclusive, rounds=rounds, offset=k,
                                                sampled=args.sampled))
                first = args.only[k % len(args.only):] + args.only[:k % len(args.only)]
                print(f"[{' -> '.join(first)}, rotated each round] x {rounds} rounds -> {script}")
            if args.submit:
                for _ in range(TIMING_JOBS):
                    submit(script, args)
            return
        # Real reads: one round per job, the second job in reverse order.
        for variant, order in (("_order1", conditions), ("_order2", conditions[::-1])):
            script = DATA_DIR / "jobs" / f"timing_{reads_set_name(args.reads)}{variant}.sbatch"
            script.write_text(timing_script(order, args.reads, real=True, variant=variant))
            print(f"[{' -> '.join(c['name'] for c in order)}] x 1 round, {REAL_THREADS} threads -> {script}")
            if args.submit:
                submit(script, args)
        return
    if args.stage == "fingerprint":
        for directory in (DATA_DIR / "jobs", DATA_DIR / "logs"):
            directory.mkdir(parents=True, exist_ok=True)
        for c in conditions:
            runs = fingerprint_runs(c, args.reads, args.run_tag)
            if not runs:
                print(f"[{c['name']}] no runs with a mapped.gam" + (f" under {args.run_tag}/" if args.run_tag else ""))
            for run in runs:
                label = "_".join([Path(args.reads).stem, c["name"]] + ([args.run_tag] if args.run_tag else []) + [run.name])
                script = DATA_DIR / "jobs" / f"fingerprint_{label}.sbatch"
                script.write_text(fingerprint_script(c, run, label))
                print(f"[{c['name']}] fingerprint {run.relative_to(results_dir(c, args.reads))} -> {script}")
                if args.submit:
                    submit(script, args)
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
            rep_tag = ("" if not args.run_tag else f"_{args.run_tag}") + ("" if rep is None else f"_rep{rep}")
            script = jobs_dir / f"{args.stage}{tag}_{c['name']}{rep_tag}.sbatch"
            text = (map_script(c, args.reads, rep, args.run_tag, args.keep_gams) if args.stage == "map"
                    else STAGES[args.stage](c, args.reads))
            script.write_text(text)
            print(f"[{c['name']}] {args.stage}{rep_tag} -> {script}")
            if args.submit:
                submit(script, args)


if __name__ == "__main__":
    main()
