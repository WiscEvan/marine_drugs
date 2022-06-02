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
    parser.add_argument(
        "--use-latest-refinement",
        help="Will attempt to use latest refinement column from table output by Automappa. If refinement columns do not exists, will default to 'cluster'.",
        action="store_true",
        required=False,
    )
    args = parser.parse_args()

    # Read in fasta to prepare for writing clusters
    records = [record for record in SeqIO.parse(args.fasta, "fasta")]

    if not os.path.isdir(args.outdir):
        print(f"outdir does not exists. Creating: {args.outdir}")
        os.makedirs(args.outdir)
    # Read in binning
    sep = ","
    if args.binning.endswith(".tsv") or args.binning.endswith(".tsv"):
        sep = "\t"
    df = pd.read_csv(args.binning, sep=sep, index_col="contig")
    # Group contigs by their cluster
    if args.use_latest_refinement:
        try:
            latest_refinement = [col for col in df.columns if "refinement" in col][-1]
            grouped_clusters = df.groupby(latest_refinement)
        except IndexError:
            print(
                "Could not find refinement columns, proceeding with 'cluster' as binning column."
            )
            grouped_clusters = df.groupby("cluster")
    else:
        grouped_clusters = df.groupby("cluster")
    n_written = 0
    for cluster, dff in grouped_clusters:
        bin_contigs = dff.index.tolist()
        bin_records = [record for record in records if record.id in bin_contigs]
        bin_filepath = os.path.join(args.outdir, f"{cluster}.fna")
        SeqIO.write(bin_records, bin_filepath, "fasta")
        n_written += 1

    print(f"Wrote {n_written} clusters to {args.outdir}")


if __name__ == "__main__":
    main()