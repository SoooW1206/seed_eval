#!/bin/bash
#SBATCH --job-name=graph-length
#SBATCH --partition=short
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/private/groups/patenlab/twang285/seed-eval/logs/graph_length.%j.log
# Total sequence length of the graph: denominator for index-based seed density.
set -euo pipefail

VG=~/bin/vg_baseline
DATA=/private/groups/patenlab/twang285/seed-eval
GBZ=$DATA/hprc-v2.0-mc-chm13-eval-sampled16o_fragmentlinked-for-real-hifi-HG002-full.gbz

/usr/bin/time -v $VG stats -l "$GBZ" > "$DATA/graph.length.txt"
cat "$DATA/graph.length.txt"
