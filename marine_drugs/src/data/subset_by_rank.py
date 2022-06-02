#!/usr/bin/env python

import os
import argparse

import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Subset counts by provided `rank`.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--taxonomy", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--rank", required=False, default="phylum")
    args = parser.parse_args()

    counts = pd.read_csv(args.counts, sep='\t', index_col='contig')
    kmer_cols = counts.columns.tolist()
    taxa = pd.read_csv(args.taxonomy, sep='\t', index_col='contig')
    df = pd.merge(counts, taxa, how='left', left_index=True, right_index=True)
    for rank_name, rank_df in df.groupby(args.rank):
        rank_name = rank_name.replace(" ", "_")
        if not os.path.isdir(args.output):
            os.makedirs(args.output)
        out = os.path.join(args.output, f"{rank_name}.tsv")
        rank_df[kmer_cols].to_csv(out, sep='\t', index=True, header=True)
        print(f"wrote {rank_df.shape[0]} contigs to {out}")


if __name__ == "__main__":
    main()