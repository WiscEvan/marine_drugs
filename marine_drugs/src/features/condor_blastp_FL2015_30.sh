#!/bin/bash

# First, copy the nr database from /staging into the working directory,
cp /staging/groups/kwan_group/nr.dmnd ./

./diamond blastp \
    --query FL2015_30.orfs.faa \
    --out FL2015_30.blastp.tsv \
    --db nr.dmnd \
    --evalue 1e-5 \
    --max-target-seqs 200 \
    --outfmt 6 \
    --threads 16 \
    --memory-limit 100

#
# Before the script exits, make sure to remove the file(s) from the working directory
rm nr.dmnd