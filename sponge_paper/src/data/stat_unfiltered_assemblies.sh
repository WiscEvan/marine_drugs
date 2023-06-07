#!/usr/bin/env bash

export PATH="$PATH:$HOME/sponge_paper/sponge_paper/src/features/"

outdir="$HOME/sponge_paper/sponge_paper/data/raw/assemblies"
metagenomes=`ls ${outdir}/*.fasta`

stat_metagenome.py --write-base-stats=${outdir} ${metagenomes[@]} "${outdir}/unfiltered_assemblies_stats.tsv"