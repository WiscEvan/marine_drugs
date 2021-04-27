#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.rnaseq_genefamilies_associations.err
#SBATCH --output=logs/%J.rnaseq_genefamilies_associations.out

REPO="/media/bigdrive2/evan"
script="${REPO}/marine_drugs/marine_drugs/src/features/pathway-profiling/sponge_phenotype_association.R"
profiles="${REPO}/marine_drugs/marine_drugs/data/processed/pathway-profiling/master_community_genefamilies_profiles.munged.tsv"
metadata="${REPO}/marine_drugs/marine_drugs/src/features/pathway-profiling/sponge_metadata.rnaseqsponges.munged.tsv"
output="${REPO}/marine_drugs/marine_drugs/data/processed/pathway-profiling/genefamilies-associations"
cores=24

Rscript $script --profiles=$profiles --metadata=$metadata --output=$output --cores=$cores