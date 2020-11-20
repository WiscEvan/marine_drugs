#!/usr/bin/env python

import argparse
import os
import glob
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Stat sponge metagenome(s)")
    parser.add_argument(
        "binning",
        help="path to binning directory (will glob *.bacteria.binning.{umap,bhsne}.tsv)",
    )
    parser.add_argument(
        "out",
        help="path to write sponge metagenomes table of statistics",
    )
    args = parser.parse_args()
    search_string = os.path.join(args.binning, "*.bacteria.binning.*.tsv")
    # We only want *.bacteria.binning.{umap,bhsne}.tsv, not bacteria.binning.hdbscan.tsv
    binning_files = [
        filepath for filepath in glob.glob(search_string) if "hdbscan" not in filepath
    ]
    completeness = 1
    purity = 2
    samples = []
    for filepath in binning_files:
        clusters = []
        sample = os.path.basename(filepath).split(".bacteria")[0]
        embed_method = os.path.basename(filepath).split(".")[3]
        df = pd.read_csv(filepath, sep="\t", index_col="contig")
        n_contigs = df.shape[0]
        for cluster, dff in df.groupby("cluster"):
            clusters.append(
                {
                    "cluster": cluster,
                    "completeness": dff.iloc[0, completeness],
                    "purity": dff.iloc[0, purity],
                }
            )
        clusters_df = pd.DataFrame(clusters)
        n_bins = clusters_df.shape[0]
        samples.append(
            {
                "sample": sample,
                "embed method": embed_method,
                "median completeness": clusters_df.completeness.median(),
                "median purity": clusters_df.purity.median(),
                "mean completeness": clusters_df.completeness.mean(),
                "mean purity": clusters_df.purity.mean(),
                "num. clusters": clusters_df.shape[0],
                "num. contigs": n_contigs,
            }
        )
    samples_df = pd.DataFrame(samples).set_index(["sample", "embed method"])
    samples_df.sort_index().to_csv(args.out, sep="\t", index=True, header=True)
    print(f"wrote: {args.out}")


if __name__ == "__main__":
    main()
