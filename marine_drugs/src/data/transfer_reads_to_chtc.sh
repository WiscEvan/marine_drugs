#!/usr/bin/env bash

assemblies="$HOME/marine_drugs/marine_drugs/data/raw/assemblies"
chtc="/home/erees/marine_drugs/marine_drugs/data/raw/assemblies"
for reads in FL2015_8_1P.00.0_0.cor.fastq.gz FL2015_8_2P.00.0_0.cor.fastq.gz;do
    rsync -azvP -e "ssh -i ~/.ssh/submit1" "${assemblies}/${reads}" erees@submit-1.chtc.wisc.edu:"${chtc}/"
done