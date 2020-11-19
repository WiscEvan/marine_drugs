#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.embed_FL2015_4_neighborhood_and_bacteria.err
#SBATCH --output=logs/%J.embed_FL2015_4_neighborhood_and_bacteria.out

# source activate sponges


norm_methods=("ilr" "am_clr")
embed_methods=("bhsne" "umap")

fasta="$HOME/marine_drugs/marine_drugs/data/interim/assemblies/FL2015_4.filtered.fna"
binning="$HOME/marine_drugs/marine_drugs/data/interim/binning"
nbr_counts="${binning}/FL2015_4.kmers.neighborhood.tsv"
bact_counts="${binning}/FL2015_4.kmers.bacteria.tsv"

data=($nbr_counts $bact_counts)
for kmers in ${data[@]};do
    filename=$(basename ${kmers})
    for norm_method in ${norm_methods[@]};do
        for embed_method in ${embed_methods[@]};do
            normalized="${binning}/${filename/.tsv/.$norm_method}.tsv"
            embedded="${normalized/.tsv/.$embed_method}.tsv"
            autometa-kmers --fasta $fasta --kmers $kmers --size 4 --norm-method $norm_method --normalized $normalized --embedded $embedded --embed-method $embed_method --pca-dimensions 50 --do-pca
        done
    done
done