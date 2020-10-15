#!/usr/bin/env bash

export PATH="$PATH:$HOME/marine_drugs/marine_drugs/src/features/"

outdir="$HOME/marine_drugs/marine_drugs/data/raw/assemblies"
metagenomes=`ls ${outdir}/*.fasta`

stat_metagenome.py ${metagenomes[@]} "${outdir}/unfiltered_assemblies_stats.tsv"