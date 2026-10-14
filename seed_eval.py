#!/usr/bin/env python3
"""Seed-method evaluator for vg giraffe.

Reads a conditions table (one seed-index configuration per row) and generates
one sbatch job per condition. Scripts are only written and printed by default;
pass --submit to actually submit them.

The driver is data-driven: every condition-specific piece of a command comes
from its row in the table, so adding a condition means adding a row.

Columns:
  name           unique condition name; used for output file names
  vg_binary      absolute path to the vg executable used for every step
  index_args     arguments passed verbatim to `vg minimizer`
  index_mem      sbatch --mem for the index job (e.g. 80G)
  density_group  only used for result tables and plots
  label          only used for result tables and plots

Currently implemented stage: index (vg minimizer -> .min + .zipcodes).
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

DATA_DIR = Path("/private/groups/patenlab/twang285/seed-eval")
GRAPH_NAME = "hprc-v2.0-mc-chm13-eval-sampled16o_fragmentlinked-for-real-hifi-HG002-full"
PARTITION = "medium"
TIME = "12:00:00"
THREADS = 16  # vg minimizer: more than 16 threads rarely helps

REQUIRED_COLUMNS = ["name", "vg_binary", "index_args", "index_mem"]
OPTIONAL_COLUMNS = ["density_group", "label"]


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


def index_script(c):
    gbz = DATA_DIR / f"{GRAPH_NAME}.gbz"
    dist = DATA_DIR / f"{GRAPH_NAME}.dist"
    out = DATA_DIR / "indexes" / c["name"]
    log = DATA_DIR / "logs" / f"index_{c['name']}.%j.log"
    return f"""#!/bin/bash
#SBATCH --job-name=idx-{c['name']}
#SBATCH --partition={PARTITION}
#SBATCH -c {THREADS}
#SBATCH --mem={c['index_mem']}
#SBATCH --time={TIME}
#SBATCH --output={log}
set -euo pipefail
mkdir -p {out.parent}
/usr/bin/time -v {c['vg_binary']} minimizer {c['index_args']} -p -t {THREADS} \\
    -d {dist} -o {out}.min.tmp -z {out}.zipcodes.tmp {gbz}
mv {out}.min.tmp {out}.min
mv {out}.zipcodes.tmp {out}.zipcodes
echo "min size (bytes): $(stat -c %s {out}.min)"
echo "zipcodes size (bytes): $(stat -c %s {out}.zipcodes)"
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("stage", choices=["index"])
    parser.add_argument("--conditions", default="conditions.tsv")
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

    for c in conditions:
        script = jobs_dir / f"{args.stage}_{c['name']}.sbatch"
        script.write_text(index_script(c))
        cmd = ["sbatch", str(script)]
        print(f"[{c['name']}] {c['index_args']}  mem={c['index_mem']}  -> {script}")
        if args.submit:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print("  " + result.stdout.strip())


if __name__ == "__main__":
    main()
