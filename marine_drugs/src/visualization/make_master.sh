#!/usr/bin/env bash

SCRIPT="$HOME/marine_drugs/marine_drugs/src/visualization/merge-sponges-binning-and-annotations.py"

contigs="master.contigs.tsv"
clusters="master.clusters.tsv"
binning="$HOME/marine_drugs/marine_drugs/src/visualization/final_binning_filepaths.tsv"
assemblies="$HOME/marine_drugs/marine_drugs/data/interim/assemblies"
annotations="/home/sam/Sponge_Proj/Final/Summary_files"

python $SCRIPT \
    --contigs $contigs \
    --clusters $clusters \
    --binning-filepaths $binning \
    --assemblies $assemblies \
    --annotations $annotations
