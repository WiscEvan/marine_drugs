#!/usr/bin/env python

import argparse
import os
import glob

from typing import Dict
from Bio import SeqIO



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", default="sponge_paper/data/interim/mitogenomes/concat")
    parser.add_argument("--outdir", default="sponge_paper/data/interim/mitogenomes/cleaned")
    args = parser.parse_args()
    search_str = os.path.join(args.indir, "*.faa")
    fps = [fp for fp in glob.glob(search_str)]
    print(f"Found {len(fps)} genes containing all 12 sponges")
    samples : Dict[str,set] = {}
    for fp in fps:
        gene = os.path.basename(fp).replace(".faa", "")
        organism_set = set()
        for rec in SeqIO.parse(fp, 'fasta'):
            organism = rec.description.replace(" ", "_")
            organism_set.add(organism)
        samples[gene] = organism_set

    popped_gene, min_orgs = samples.popitem()
    for gene,organism_set in samples.items():
        min_orgs = min_orgs.intersection(organism_set)

    filtered_records = {}
    for fp in fps:
        gene = os.path.basename(fp).replace(".faa", "")
        records = []
        for rec in SeqIO.parse(fp, 'fasta'):
            organism = rec.description.replace(" ", "_")
            if organism not in min_orgs:
                continue
            rec.id = organism
            rec.name = organism
            rec.description = ''
            records.append(rec)
        filtered_records[gene] = records
    
    print(f"{len(min_orgs)} organisms across {len(filtered_records)} genes (all containing all sponge mtdna)")

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)

    for gene,records in filtered_records.items():
        outfname = f"{gene}.faa"
        outfpath = os.path.join(args.outdir, outfname)
        n_written = SeqIO.write(records, outfpath, 'fasta')
        print(f"{n_written} records to {outfname}")

if __name__ == '__main__':
    main()