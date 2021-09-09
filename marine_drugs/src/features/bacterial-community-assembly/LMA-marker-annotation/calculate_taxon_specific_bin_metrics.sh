#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.calculate_taxon_specific_bin_metrics.err
#SBATCH --output=logs/%J.calculate_taxon_specific_bin_metrics.out



HMMSCAN_FILES="${HOME}/marine_drugs/marine_drugs/data/interim/LMA-marker-annotation/*.hmmscan.tsv"
BIN_FILES="${HOME}/marine_drugs/marine_drugs/data/interim/binning/final_binning_filepaths.tsv"
SRC="${HOME}/marine_drugs/marine_drugs/src/features/bacterial-community-assembly/LMA-marker-annotation/calculate_taxon_specific_bin_metrics.py"

for hmmscan in `ls $HMMSCAN_FILES`;do
    # --hmmscan --> path to hmmscan.tsv output on taxon-specific marker set
    sponge=$(basename ${hmmscan} | cut -f1,2 -d"_")
    # --binning --> path to binning.tsv table
    # e.g. FL2015-4        /path/to/binning.tsv
    binning=$(grep "${sponge/_/-}\s" $BIN_FILES | awk '{print $2}')
    
    # --output-markers --> path to write markers.tsv output of taxon-specific marker set
    markers="${hmmscan/.hmmscan/.markers}"
    # --output-metrics --> path to write binning_metrics.tsv
    metrics="${hmmscan/.hmmscan/.metrics}"

    python $SRC --hmmscan $hmmscan --binning $binning --output-markers $markers --output-metrics $metrics

done