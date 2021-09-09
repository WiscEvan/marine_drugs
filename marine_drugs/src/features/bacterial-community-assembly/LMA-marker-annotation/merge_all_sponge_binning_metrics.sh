#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.merge_all_sponge_binning_metrics.err
#SBATCH --output=logs/%J.merge_all_sponge_binning_metrics.out



HMMSCAN_FILES="${HOME}/marine_drugs/marine_drugs/data/interim/LMA-marker-annotation/*.hmmscan.tsv"
BIN_FILES="${HOME}/marine_drugs/marine_drugs/data/interim/binning/final_binning_filepaths.tsv"
SRC="${HOME}/marine_drugs/marine_drugs/src/features/bacterial-community-assembly/LMA-marker-annotation/merge_all_sponge_binning_metrics.py"

INDIR="/home/evan/marine_drugs/marine_drugs/data/interim/LMA-marker-annotation"
OUTDIR="/home/evan/marine_drugs/marine_drugs/data/processed/LMA-marker-annotation"
if [ ! -d $OUTDIR ]
then mkdir -p $OUTDIR
fi

for rank in phylum class;do
    out="${OUTDIR}/LMA_sponges_${rank}_markers_cluster_metrics.tsv"
    python $SRC --indir $INDIR --rank $rank --out $out
done