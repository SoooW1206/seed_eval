#!/bin/bash
#SBATCH --job-name=minimaxmer-proto
#SBATCH --partition=short
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
# Run the minimaxmer prototype (k=31, w=20) on the reads of a GAM, twice:
#   constant: all counts constant, so the selection must equal vg's (w, k)-minimizers;
#   counted:  reads_count from these reads, genome_count 1.
# Read sequences only pass through pipes; the outputs hold read names and numbers.
# Usage, from the seed-eval repo directory:
#   sbatch -o LOG prototype/run_on_reads.sh TRUTH.gam OUT_DIR VG
set -euo pipefail
GAM=$1
OUT=$2
VG=$3
PROTO="${SLURM_SUBMIT_DIR:-$PWD}/prototype/minimaxmer.py"
for tool in jq python3; do
    command -v $tool > /dev/null || { echo "$tool not found on $(hostname)" >&2; exit 1; }
done
[ -f "$PROTO" ] || { echo "$PROTO not found; submit from the seed-eval repo directory" >&2; exit 1; }
mkdir -p "$OUT"
for mode in constant counted; do
    flag=""
    [ "$mode" = constant ] && flag="--constant-counts"
    "$VG" view -aj "$GAM" | jq -r '[.name, .sequence] | @tsv' \
        | /usr/bin/time -v -o "$OUT/$mode.time" python3 -I "$PROTO" --k 31 --w 20 $flag \
            --per-read "$OUT/per_read_$mode.tsv.tmp" > "$OUT/summary_$mode.txt.tmp"
    mv "$OUT/per_read_$mode.tsv.tmp" "$OUT/per_read_$mode.tsv"
    mv "$OUT/summary_$mode.txt.tmp" "$OUT/summary_$mode.txt"
    echo "== $mode"
    cat "$OUT/summary_$mode.txt"
    grep -E "Elapsed|Maximum resident" "$OUT/$mode.time"
done
