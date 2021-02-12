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
    embedding = pd.read_csv(args.embedded, sep="\t", index_col="contig")

    df = pd.merge(embedding, taxonomy, left_index=True, right_index=True, how="left")

    bacteria = df[df["superkingdom"] == "bacteria"]
    bacteria_counts = counts.loc[counts.index.isin(bacteria.index)]
    bacteria_counts.to_csv(args.bacteria_out, sep="\t", index=True, header=True)
    # See notebook 02-WiscEvan-kmer-embedding-methods-comparisons.ipynb
    neighborhood = df[df["x"] > 5]
    neighborhood_counts = counts.loc[counts.index.isin(neighborhood.index)]
    neighborhood_counts.to_csv(args.neighborhood_out, sep="\t", index=True, header=True)


if __name__ == "__main__":
    main()
