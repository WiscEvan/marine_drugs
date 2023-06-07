#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.rnaseq_genefamilies_associations.err
#SBATCH --output=logs/%J.rnaseq_genefamilies_associations.out

REPO="/media/bigdrive2/evan"
script="${REPO}/sponge_paper/sponge_paper/src/features/pathway-profiling/sponge_phenotype_association.R"
profiles="${REPO}/sponge_paper/sponge_paper/data/processed/pathway-profiling/master_community_genefamilies_profiles.munged.tsv"
metadata="${REPO}/sponge_paper/sponge_paper/src/features/pathway-profiling/sponge_metadata.rnaseqsponges.munged.tsv"
output="${REPO}/sponge_paper/sponge_paper/data/processed/pathway-profiling/genefamilies-associations"
cores=24

Rscript $script --profiles=$profiles --metadata=$metadata --output=$output --cores=$cores