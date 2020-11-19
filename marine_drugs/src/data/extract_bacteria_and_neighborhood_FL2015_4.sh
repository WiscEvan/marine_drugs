#!/usr/bin/env bash

# 
#   --counts COUNTS       File with which to extract sequences
#   --embedded EMBEDDED   Path to embedded kmers file
#   --taxonomy TAXONOMY   Path to taxonomy file
#   --bacteria-out BACTERIA_OUT
#                         File to write out subsetted counts
#   --neighborhood-out NEIGHBORHOOD_OUT
#                         File to write out subsetted counts

src="/home/evan/marine_drugs/marine_drugs/src/data/extract_bacteria_and_neighborhood_FL2015_4.py"
binning="$HOME/marine_drugs/marine_drugs/data/interim/binning"
counts="${binning}/FL2015_4.kmers.tsv"
embedded="${binning}/FL2015_4.kmers.ilr.umap.tsv"
taxonomy="${binning}/FL2015_4.taxonomy.tsv"
bacteria_out="${binning}/FL2015_4.kmers.bacteria.tsv"
neighborhood_out="${binning}/FL2015_4.kmers.neighborhood.tsv"
python $src \
    --counts $counts \
    --embedded $embedded \
    --taxonomy $taxonomy \
    --bacteria-out $bacteria_out \
    --neighborhood-out $neighborhood_out