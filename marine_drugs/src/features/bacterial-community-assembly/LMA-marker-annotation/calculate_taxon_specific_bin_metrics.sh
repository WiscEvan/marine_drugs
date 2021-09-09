#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.calculate_taxon_specific_bin_metrics.err
#SBATCH --output=logs/%J.calculate_taxon_specific_bin_metrics.out



SRC="${HOME}/marine_drugs/marine_drugs/src/features/bacterial-community-assembly/LMA-marker-annotation/calculate_taxon_specific_bin_metrics.py"
DATA="${HOME}/marine_drugs/marine_drugs/data/interim/binning"
OUTDIR="${HOME}/marine_drugs/marine_drugs/data/interim/LMA-marker-annotation"
CPUS=8


if [ ! -d $OUTDIR ]
then mkdir -p $OUTDIR
fi

EXTERNAL="${HOME}/marine_drugs/marine_drugs/data/external"
for hmmscan in `~/marine_drugs/marine_drugs/data/interim/LMA-marker-annotation/*.hmmscan.tsv`;do
    # --hmmscan --> path to hmmscan.tsv output on taxon-specific marker set
    sponge=$(basename ${hmmscan} | cut -f1,2 -d"_")
    
    # --markers --> path to write markers.tsv output of taxon-specific marker set
    markers="${hmmscan/.hmmscan/.markers}"
    # --metrics --> path to write binning_metrics.tsv
    metrics="${hmmscan/.hmmscan/.metrics}"
    
    # --binning --> path to binning.tsv table
    binning="${DATA}/${sponge}/"
    
    
    python $SRC --hmmscan $hmmscan --markers $markers --binning $binning --metrics $metrics
done