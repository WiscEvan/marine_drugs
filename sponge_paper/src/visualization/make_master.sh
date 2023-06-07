#!/usr/bin/env bash

SCRIPT="$HOME/sponge_paper/sponge_paper/src/visualization/merge-sponges-binning-and-annotations.py"

contigs="master.contigs.tsv"
clusters="master.clusters.tsv"
binning="$HOME/sponge_paper/sponge_paper/src/visualization/final_binning_filepaths.tsv"
assemblies="$HOME/sponge_paper/sponge_paper/data/interim/assemblies"
annotations="/home/sam/Sponge_Proj/Final/Summary_files"

python $SCRIPT \
    --contigs $contigs \
    --clusters $clusters \
    --binning-filepaths $binning \
    --assemblies $assemblies \
    --annotations $annotations
