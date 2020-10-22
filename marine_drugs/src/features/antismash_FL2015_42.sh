#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.antismash_FL2015_42.err
#SBATCH --output=logs/%J.antismash_FL2015_42.out



outdir="$HOME/marine_drugs/marine_drugs/data/interim/bgcs/FL2015_42"
fasta="$HOME/marine_drugs/marine_drugs/data/raw/assemblies/FL2015_42.fasta"
cpus=20
# https://phdops.kblin.org/2015-running-antismash-standalone-from-docker.html
run_antismash $fasta $outdir --genefinding-tool prodigal-m --cf-create-clusters --clusterhmmer  --cpus $cpus