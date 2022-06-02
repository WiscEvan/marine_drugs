#!/usr/bin/env bash

# conda install -c bioconda checkm-genome -y
# conda activate sponges

# To look-up rank and taxon use checkm taxon_list > taxon_list.txt and then open taxon_list.txt

SRC="/home/evan/marine_drugs/marine_drugs/src/features/bacterial-community-assembly/LMA-marker-annotation/retrieve_checkm_taxon_set_markers.py"
CHECKM_HMM="/home/evan/miniconda3/envs/sponges/checkm_data/hmms/checkm.hmm"
EXTERNAL="/home/evan/marine_drugs/marine_drugs/data/external"

rank="phylum"
taxon="Proteobacteria"
outfile="${EXTERNAL}/${rank}_${taxon}_markers.txt"
if [ ! -f $outfile ]
then
    checkm taxon_set $rank $taxon $outfile
else
    echo "$outfile found, skipping checkm taxon_set ..."
fi
python $SRC --taxon-set $outfile --outdir ${outfile/.txt/} --checkm-hmm $CHECKM_HMM

rank="class"
taxon="Alphaproteobacteria"
outfile="${EXTERNAL}/${rank}_${taxon}_markers.txt"
if [ ! -f $outfile ]
then
    checkm taxon_set $rank $taxon $outfile
else
    echo "$outfile found, skipping checkm taxon_set ..."
fi
python $SRC --taxon-set $outfile --outdir ${outfile/.txt/} --checkm-hmm $CHECKM_HMM

rank="class"
taxon="Gammaproteobacteria"
outfile="${EXTERNAL}/${rank}_${taxon}_markers.txt"
if [ ! -f $outfile ]
then
    checkm taxon_set $rank $taxon $outfile
else
    echo "$outfile found, skipping checkm taxon_set ..."
fi
python $SRC --taxon-set $outfile --outdir ${outfile/.txt/} --checkm-hmm $CHECKM_HMM

### now press all generated/retrieved HMM models
for hmm in `ls ${EXTERNAL}/*/*.hmm`;do
    hmmpress -f $hmm
done
