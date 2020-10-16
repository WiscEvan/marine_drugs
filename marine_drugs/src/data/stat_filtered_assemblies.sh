#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.filtered_assemblies_stats.err
#SBATCH --output=logs/%J.filtered_assemblies_stats.out


export PATH="$PATH:$HOME/marine_drugs/marine_drugs/src/features/"

outdir="$HOME/marine_drugs/marine_drugs/data/interim/assemblies"
metagenomes=`ls ${outdir}/*.fna`

stat_metagenome.py ${metagenomes[@]} "${outdir}/filtered_assemblies_stats.tsv"