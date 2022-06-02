#!/usr/bin/env bash

binning="$HOME/marine_drugs/marine_drugs/data/interim/binning"
chtc="/home/erees/marine_drugs/marine_drugs/data/interim/binning"
rsync -azP -e "ssh -i ~/.ssh/submit1" erees@submit-1.chtc.wisc.edu:"${chtc}/*.blastp.tsv" "${binning}/."
