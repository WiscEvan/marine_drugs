#!/usr/bin/env python

import os
import argparse
import pandas as pd
from typing import List
from Bio import SeqIO
from Bio.Seq import Seq


def extract_seqrecords(gbk:str, organism_name: str)->List[SeqIO.SeqRecord]:
    records = []
    for record in SeqIO.parse(gbk, "genbank"):
        for feat in record.features:
            if feat.type != "CDS":
                continue
            protein_id = feat.qualifiers.get('protein_id')[0]
            product = feat.qualifiers.get('product')[0]
            seq = Seq(feat.qualifiers.get("translation")[0])
            description = f"{protein_id} - {product} - {organism_name}"
            seqrecord = SeqIO.SeqRecord(seq, id=protein_id, name=product, description=description)
            records.append(seqrecord)
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gbk", help="Path to gbk file", required=True, nargs='+')
    parser.add_argument("--gbk-orgs", help="Path to gbk_organisms.tsv table", required=True)
    parser.add_argument("--outdir", help="Directory path to write fasta files", required=True)
    args = parser.parse_args()
    # Accession, Organism
    df = pd.read_table(args.gbk_orgs, index_col='Accession')
    for gbk in args.gbk:
        # get organism name
        accession = os.path.basename(gbk).replace(".gbk", "")
        try:
            organism_name = df.loc[accession, "Organism"]
        except KeyError:
            print(f"Failed to find organism name for {accession}")
        # get seqs
        records = extract_seqrecords(gbk, organism_name)
        # filepath handling
        outname = os.path.basename(gbk).replace(".gbk", ".faa")
        outfpath = os.path.join(args.outdir, outname)
        # write fasta
        n_written = SeqIO.write(records, outfpath, 'fasta')
        print(f"Wrote {n_written} CDS to {outfpath}")

if __name__ == '__main__':
    main()