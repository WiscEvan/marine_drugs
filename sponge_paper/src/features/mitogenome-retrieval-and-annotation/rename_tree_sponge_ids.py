#!/usr/bin/env python

import argparse
import os
import pandas as pd

from Bio import SeqIO
from Bio import Phylo


def main():
	parser = argparse.ArgumentParser(description="Convert sponge ID tree terminals to their scientific names.")
	parser.add_argument("--tree", help="Path to tree file", required=True)
	parser.add_argument("--format", help="Format of tree file", required=False, default="newick")
	parser.add_argument("--out", help="Path to write renamed tree file", required=True)
	parser.add_argument(    "--out-format", help="Format of output tree file", required=False, default="newick")
	args = parser.parse_args()
	tree = Phylo.read(args.tree, args.format)
	names = [terminal.name for terminal in tree.get_terminals()]
	genomes = {}
	for name in names:
		match = re.search(r"(GCA_\d+)\.\d?", name)
		if match:
			genome = match.group()
			print(genome)
if __name__ == "__main__":
	main()
