#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.get_taxonomy_FL2015_44.err
#SBATCH --output=logs/%J.get_taxonomy_FL2015_44.out

autometa-taxonomy --ncbi /mnt/autometa_databases/ --assembly /home/evan/sponge_paper/sponge_paper/data/interim/assemblies/FL2015_44.filtered.fna --nucl-orfs /home/evan/sponge_paper/sponge_paper/data/interim/binning/FL2015_44.orfs.fna --prot-orfs /home/evan/sponge_paper/sponge_paper/data/interim/binning/FL2015_44.orfs.faa --blast /home/evan/sponge_paper/sponge_paper/data/interim/binning/FL2015_44.blastp.tsv --lca /home/evan/sponge_paper/sponge_paper/data/interim/binning/FL2015_44.lca.tsv --method majority_vote /home/evan/sponge_paper/sponge_paper/data/interim/binning/FL2015_44.taxonomy.tsv