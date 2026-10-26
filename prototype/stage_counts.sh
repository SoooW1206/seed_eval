#!/bin/bash
#SBATCH --job-name=stage-counts
#SBATCH --partition=medium
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
# Work done per read in each giraffe stage, from the provenance GAM of density runs: every
# annotation.stage.<path>.results count (items coming out of the stage) and the aligner's per-read totals
# (bases aligned, invocations, fallbacks; not its times), one JSON line per read with the read name and
# numbers only (no sequences, no timings). Input for prototype/stage_work.py.
# Usage: sbatch -o LOG prototype/stage_counts.sh RESULTS_DIR VG CONDITION ...
#   writes RESULTS_DIR/<condition>/density/stage_counts.jsonl
set -euo pipefail
R=$1
VG=$2
shift 2
command -v jq > /dev/null || { echo "jq not found" >&2; exit 1; }
for c in "$@"; do
    "$VG" view -aj $R/$c/density/prov.gam | jq -c \
        '{read: .name} + ([.annotation.stage // {} | paths(type == "number") as $p | select($p[-1] == "results")
                           | {key: ($p[:-1] | join(".")), value: getpath($p)}] | from_entries)
        + ([.annotation.aligner_stats.per_read // {} | paths(type == "number") as $p
            | select($p[-1] == "total" and ($p | index("time")) == null)
            | {key: ("aligner." + ($p[:-1] | join("."))), value: getpath($p)}] | from_entries)' \
        > $R/$c/density/stage_counts.jsonl.tmp
    mv $R/$c/density/stage_counts.jsonl.tmp $R/$c/density/stage_counts.jsonl
    echo "$c: $(wc -l < $R/$c/density/stage_counts.jsonl) reads"
done
