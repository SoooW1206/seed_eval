#!/bin/bash
#SBATCH --job-name=build-dist
#SBATCH --partition=medium
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/patenlab/twang285/seed-eval/logs/build_dist.%j.log
# Build the distance index shared by all seed conditions; record peak memory/time.
set -euo pipefail

VG=~/bin/vg_baseline
DATA=/private/groups/patenlab/twang285/seed-eval
NAME=hprc-v2.0-mc-chm13-eval-sampled16o_fragmentlinked-for-real-hifi-HG002-full
GBZ=$DATA/$NAME.gbz
DIST=$DATA/$NAME.dist

/usr/bin/time -v $VG index -t "$SLURM_CPUS_PER_TASK" -j "$DIST.tmp" "$GBZ"
mv "$DIST.tmp" "$DIST"
echo "dist index size (bytes): $(stat -c %s "$DIST")"
