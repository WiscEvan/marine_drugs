#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.filtered_assemblies_stats.err
#SBATCH --output=logs/%J.filtered_assemblies_stats.out


export PATH="$PATH:$HOME/sponge_paper/sponge_paper/src/features/"

outdir="$HOME/sponge_paper/sponge_paper/data/interim/assemblies"
metagenomes=`ls ${outdir}/*.fna`

stat_metagenome.py ${metagenomes[@]} "${outdir}/filtered_assemblies_stats.tsv"