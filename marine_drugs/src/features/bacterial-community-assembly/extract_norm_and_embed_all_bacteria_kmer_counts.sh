#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.extract_norm_and_embed_all_bacteria_kmer_counts.err
#SBATCH --output=logs/%J.extract_norm_and_embed_all_bacteria_kmer_counts.out

# source activate sponges


embed_methods=("bhsne" "umap")

assemblies="$HOME/marine_drugs/marine_drugs/data/interim/assemblies"
binning="$HOME/marine_drugs/marine_drugs/data/interim/binning"
src="$HOME/marine_drugs/marine_drugs/src/features/subset_kmer_counts.py"

for fasta in `ls ${assemblies}/*.filtered.fna`;do
    # Get filepaths corresponding to sponge.
    filename=$(basename ${fasta})
    sponge=${filename/.filtered.fna/}
    counts="${binning}/${sponge}.kmers.tsv"
    taxonomy="${binning}/${sponge}.taxonomy.tsv"
    kmers="${binning}/${sponge}.kmers.bacteria.tsv"
    if [ ! -f $kmers ];
    then python $src --counts $counts --taxonomy $taxonomy --out $kmers
    else echo "$kmers already exists. Skipping..."
    fi
    # Now perform normalization and embedding on bacteria counts.
    for embed_method in ${embed_methods[@]};do
        normalized="${kmers/.tsv/.am_clr.tsv}"
        embedded="${normalized/.tsv/.$embed_method}.tsv"
        if [ ! -f $embedded ];
        then autometa-kmers --fasta $fasta --kmers $kmers --size 4 --norm-method am_clr --normalized $normalized --embedded $embedded --embed-method $embed_method --pca-dimensions 50 --do-pca
        else echo "$embedded already exists. skipping..."
        fi
    done
done