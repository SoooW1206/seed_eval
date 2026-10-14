#!/bin/bash
#SBATCH --job-name=check-paths
#SBATCH --partition=short
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/private/groups/patenlab/twang285/seed-eval/logs/check_paths.%j.log
# List path metadata of the eval graph and summarize which sample is REFERENCE.
set -euo pipefail

VG=~/bin/vg_baseline
DATA=/private/groups/patenlab/twang285/seed-eval
GBZ=$DATA/hprc-v2.0-mc-chm13-eval-sampled16o_fragmentlinked-for-real-hifi-HG002-full.gbz
OUT=$DATA/graph.paths_meta.tsv

/usr/bin/time -v $VG paths -M -x "$GBZ" > "$OUT"

echo "=== path counts by SENSE, SAMPLE"
awk 'NR>1 {print $2"\t"$3}' "$OUT" | sort | uniq -c
echo "=== REFERENCE paths (first 40)"
awk 'NR==1 || $2=="REFERENCE"' "$OUT" | head -41
