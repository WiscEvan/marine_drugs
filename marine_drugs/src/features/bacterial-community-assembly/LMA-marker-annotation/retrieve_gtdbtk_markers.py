#!/usr/bin/env python

import argparse
import os
import pandas as pd
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtdb-taxonomy", default="/media/bigdrive1/Databases/release202/taxonomy/gtdb_taxonomy.tsv", help='GTDB release taxonomy/gtdb_taxonomy.tsv')
    parser.add_argument("--taxon-filter", default="p__proteobacteria", help='taxon to filter for retrieval of accesssions')
    # c__alphaproteobacteria
    # c__gammaproteobacteria
    args = parser.parse_args()

    df = pd.read_csv(args.gtdb_taxonomy, sep='\t', header=None).rename(columns={0:"accession", 1:"taxonomy"}).set_index("accession")
    dff = df.taxonomy.str.split(";").explode().str.lower()

    dff = dff[dff.str.contains(args.taxon_filter)]

    records = [record for record in SeqIO.parse(args.msa, 'fasta') if record.id in dff.index]




if __name__ == '__main__':
    main()