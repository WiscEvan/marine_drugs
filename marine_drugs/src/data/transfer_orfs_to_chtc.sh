#!/usr/bin/env bash

declare -a samples
samples=(FL2014_9 FL2015_43 FL2015_5 FL2015_9)
binning="$HOME/marine_drugs/marine_drugs/data/interim/binning"
chtc="/home/erees/marine_drugs/marine_drugs/data/interim/binning"
for sample in ${samples[@]};do
    orfs="${binning}/${sample}.orfs.faa"
    rsync -azvP -e "ssh -i ~/.ssh/submit1" $orfs erees@submit-1.chtc.wisc.edu:"${chtc}/"
done