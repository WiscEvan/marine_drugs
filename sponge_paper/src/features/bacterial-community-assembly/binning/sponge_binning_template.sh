#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.sample_bacteria_binning.err
#SBATCH --output=logs/%J.sample_bacteria_binning.out

# Shared Data
binning="$HOME/sponge_paper/sponge_paper/data/interim/binning"
coverages="${binning}/sample.coverages.tsv"
taxonomy="${binning}/sample.taxonomy.tsv"
norm_kmers="${binning}/sample.kmers.bacteria.am_clr.tsv"

embed_methods=("bhsne" "umap")
for embed_method in ${embed_methods[@]};do
    embedded_kmers="${binning}/sample.kmers.am_clr.${embed_method}.tsv"
    # Bacteria specific
    out="${binning}/sample.bacteria.binning.${embed_method}.tsv"
    markers="${binning}/sample.bacteria.markers.tsv"
    /home/evan/miniconda3/envs/autometa/bin/autometa-binning $norm_kmers $coverages $markers $out \
        --embedded-kmers $embedded_kmers \
        --taxonomy $taxonomy \
        --clustering-method hdbscan \
        --domain bacteria
done