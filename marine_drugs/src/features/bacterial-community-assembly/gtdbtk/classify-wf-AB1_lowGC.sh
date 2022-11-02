#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.gtdbtk_classify_wf_AB1_lowGC.err
#SBATCH --output=logs/%J.gtdbtk_classify_wf_AB1_lowGC.out

# source activate sponges


genomeDir="${HOME}/marine_drugs/marine_drugs/data/external"
extension="fasta"
outdir="${HOME}/marine_drugs/marine_drugs/data/processed/gtdbtk"
gtdbtk classify_wf --genome_dir $genomeDir --extension $extension --out_dir $outdir