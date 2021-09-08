#!/usr/bin/env bash

# conda install -c bioconda checkm-genome -y
# conda activate sponges

# To look-up rank and taxon use checkm taxon_list > taxon_list.txt and then open taxon_list.txt

SRC="/home/evan/marine_drugs/marine_drugs/src/features/bacterial-community-assembly/LMA-marker-annotation/retrieve_checkm_taxon_set_markers.py"
CHECKM_HMM="/home/evan/miniconda3/envs/sponges/checkm_data/hmms/checkm.hmm"


rank="phylum"
taxon="Proteobacteria"
outfile="/home/evan/marine_drugs/marine_drugs/data/external/${rank}_${taxon}_markers.txt"
checkm taxon_set $rank $taxon $outfile
python $SRC --taxon-set $outfile --outdir ${outfile/.txt/} --checkm-hmm $CHECKM_HMM

rank="class"
taxon="Alphaproteobacteria"
outfile="/home/evan/marine_drugs/marine_drugs/data/external/${rank}_${taxon}_markers.txt"
checkm taxon_set $rank $taxon $outfile
python $SRC --taxon-set $outfile --outdir ${outfile/.txt/} --checkm-hmm $CHECKM_HMM

rank="class"
taxon="Gammaproteobacteria"
outfile="/home/evan/marine_drugs/marine_drugs/data/external/${rank}_${taxon}_markers.txt"
checkm taxon_set $rank $taxon $outfile
python $SRC --taxon-set $outfile --outdir ${outfile/.txt/} --checkm-hmm $CHECKM_HMM