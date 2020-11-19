#!/usr/bin/env python

import argparse
import os

import pandas as pd
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(
        description="Extract bins from `fasta` using provided `binning` file into `outdir`."
    )
    parser.add_argument("--binning", help="Autometa binning file", required=True)
    parser.add_argument(
        "--fasta", help="File with which to extract sequences", required=True
    )
    parser.add_argument(
        "--outdir",
        help="Directory path to write clusters. (will create `outdir` if it does not exists)",
        required=True,
    )
    args = parser.parse_args()

    # Read in fasta to prepare for writing clusters
    records = [record for record in SeqIO.parse(args.fasta, "fasta")]

    if not os.path.isdir(args.outdir):
        print(f"outdir does not exists. Creating: {args.outdir}")
        os.makedirs(args.outdir)
    # Read in binning
    df = pd.read_csv(args.binning, sep="\t", index_col="contig")
    # Group contigs by their cluster
    for cluster, dff in df.groupby("cluster"):
        bin_contigs = dff.index.tolist()
        bin_records = [record for record in records if record.id in bin_contigs]
        bin_filepath = os.path.join(args.outdir, f"{cluster}.fna")
        SeqIO.write(bin_records, bin_filepath, "fasta")

    print(f"Wrote {df.cluster.nunique()} clusters to {args.outdir}")


if __name__ == "__main__":
    main()