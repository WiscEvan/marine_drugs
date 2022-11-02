#!/usr/bin/env python

import argparse
import os
import pandas as pd

from Bio import Phylo

def rename_sponge_terminals(tree, sponge_metadata_df):
    for terminal in tree.get_terminals():
        if "FL" in terminal.name:
            taxon = sponge_metadata_df.loc[terminal.name, 'Taxonomic Classification']
            taxon = taxon.replace('(determined w/mtDNA)', '').strip()
            sponge_id = terminal.name.replace('_', '-')
            terminal_name = f"{taxon} - {sponge_id}"
            terminal.name = terminal_name
    return tree

def main():
    parser = argparse.ArgumentParser(description="Convert sponge ID tree terminals to their scientific names.")
    parser.add_argument("--tree", help="Path to tree file", required=True)
    parser.add_argument("--format", help="Format of tree file", required=False, default="newick")
    parser.add_argument("--out", help="Path to write renamed tree file", required=True)
    parser.add_argument("--out-format", help="Format of output tree file", required=False, default="newick")
    parser.add_argument("--metadata", help="Path to sponge metadata table", required=True)
    # args.metadata = marine_drugs/data/raw/sponge_metadata.tsv
    args = parser.parse_args()
    tree = Phylo.read(args.tree, args.format)
    df = pd.read_table(args.metadata, index_col='Sponge specimen')
    df.index = df.index.map(lambda x: x.replace('-', '_'))

    tree = rename_sponge_terminals(tree, df)
    Phylo.write(tree, args.out, args.out_format)

if __name__ == "__main__":
    main()
