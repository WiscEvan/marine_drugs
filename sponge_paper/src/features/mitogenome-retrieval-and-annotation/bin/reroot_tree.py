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
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Convert sponge ID tree terminals to their scientific names."
    )
    parser.add_argument("--tree", help="Path to tree file", required=True)
    parser.add_argument("--format", help="Format of tree file", required=False, default="newick")
    parser.add_argument("--outgroup-targets",
        help="Reroot this tree with the outgroup clade containing outgroup_targets.",
        nargs='+',
        default=['Rhizopus_arrhizus', 'Podila_verticillata', 'Allomyces_macrogynus'],
    )
    parser.add_argument("--out", help="Path to write renamed tree file", required=True)
    parser.add_argument("--out-format", help="Format of output tree file", required=False, default="newick")
    args = parser.parse_args()
    tree = Phylo.read(args.tree, args.format)
    print(f"attempting to root tree using outgroup targets: {', '.join(args.outgroup_targets)}")
    tree.root_with_outgroup(args.outgroup_targets)
    Phylo.write(tree, args.out, args.out_format)
    print(f"Wrote rooted tree to {args.out} in {args.out_format} format")

if __name__ == "__main__":
    main()
