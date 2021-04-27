#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.rnaseq_pathabundance_associations.err
#SBATCH --output=logs/%J.rnaseq_pathabundance_associations.out




REPO="/media/bigdrive2/evan"
script="${REPO}/marine_drugs/marine_drugs/src/features/pathway-profiling/sponge_phenotype_association.R"
profiles="${REPO}/marine_drugs/marine_drugs/data/processed/pathway-profiling/master_community_pathabundance_profiles.munged.tsv"
metadata="${REPO}/marine_drugs/marine_drugs/src/features/pathway-profiling/sponge_metadata.rnaseqsponges.munged.tsv"
output="${REPO}/marine_drugs/marine_drugs/data/processed/pathway-profiling/pathabundance-associations"
cores=24

if [[ ! -f $profiles ]]
then
    # Required prior to running:
    python munge_profiles.py \
        --input "${REPO}/marine_drugs/marine_drugs/data/processed/pathway-profiling/master_community_pathabundance_profiles.tsv" \
        --output $profiles
fi

Rscript $script --profiles=$profiles --metadata=$metadata --output=$output --cores=$cores