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
"""
import argparse
import csv
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
MAP_MEM = "128G"  # GBZ ~11.5 GiB + dist 5.3 GB + min ~24 GB + zipcodes ~1 GB, plus working memory
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


def sbatch_header(job_name, mem, time, log):
    return f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={PARTITION}
#SBATCH -c {THREADS}
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
        print "name", "reads", "total_bp", "minimizers", "seeds", "density_pooled", "density_mean_per_read", "seeds_per_bp"
        print name, n, len, m, s, m / len, d / n, s / len }}' $O/per_read.tsv > $O/density_summary.tsv
cat $O/density_summary.tsv
"""


STAGES = {"index": index_script, "map": map_script, "density": density_script}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("stage", choices=sorted(STAGES))
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
