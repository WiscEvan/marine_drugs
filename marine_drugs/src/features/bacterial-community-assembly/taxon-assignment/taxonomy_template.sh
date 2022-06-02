#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.get_taxonomy_sample.err
#SBATCH --output=logs/%J.get_taxonomy_sample.out

autometa-taxonomy --ncbi /mnt/autometa_databases/ --assembly metagenome.filtered.fna --nucl-orfs orfs.fna --prot-orfs orfs.faa --blast blastp.tsv --lca lca.tsv --method majority_vote taxonomy.tsv