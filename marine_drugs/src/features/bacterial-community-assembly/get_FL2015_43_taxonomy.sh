#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.get_taxonomy_FL2015_43.err
#SBATCH --output=logs/%J.get_taxonomy_FL2015_43.out

autometa-taxonomy --ncbi /mnt/autometa_databases/ --assembly /home/evan/marine_drugs/marine_drugs/data/interim/assemblies/FL2015_43.filtered.fna --nucl-orfs /home/evan/marine_drugs/marine_drugs/data/interim/binning/FL2015_43.orfs.fna --prot-orfs /home/evan/marine_drugs/marine_drugs/data/interim/binning/FL2015_43.orfs.faa --blast /home/evan/marine_drugs/marine_drugs/data/interim/binning/FL2015_43.blastp.tsv --lca /home/evan/marine_drugs/marine_drugs/data/interim/binning/FL2015_43.lca.tsv --method majority_vote /home/evan/marine_drugs/marine_drugs/data/interim/binning/FL2015_43.taxonomy.tsv