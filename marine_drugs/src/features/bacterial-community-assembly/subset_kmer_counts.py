#!/usr/bin/env python

import argparse
import os

import pandas as pd
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(
        description="Extract bins from `fasta` using provided `binning` file into `outdir`."
    )
    parser.add_argument(
        "--counts", help="File with which to extract sequences", required=True
    )
    parser.add_argument("--taxonomy", help="Path to taxonomy file", required=True)
    parser.add_argument(
        "--out", help="File to write out subsetted counts", required=True
    )

    args = parser.parse_args()

    counts = pd.read_csv(args.counts, sep="\t", index_col="contig")
    taxonomy = pd.read_csv(args.taxonomy, sep="\t", index_col="contig")

    bacteria = taxonomy[taxonomy["superkingdom"] == "bacteria"]
    bacteria_counts = counts.loc[counts.index.isin(bacteria.index)]
    bacteria_counts.to_csv(args.out, sep="\t", index=True, header=True)
    print(f"Wrote {bacteria_counts.shape[0]} bacterial contigs to {args.out}")


if __name__ == "__main__":
    main()
