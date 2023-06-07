#!/usr/bin/env bash

binning="$HOME/sponge_paper/sponge_paper/data/interim/binning"
chtc="/home/erees/sponge_paper/sponge_paper/data/interim/binning"
rsync -azP -e "ssh -i ~/.ssh/submit1" erees@submit-1.chtc.wisc.edu:"${chtc}/*.blastp.tsv" "${binning}/."
