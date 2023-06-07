#!/usr/bin/env python

import argparse
import glob
import os

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--indir",
        help="Path to directory containing sponges *metrics.tsv files",
        required=True,
    )
    parser.add_argument("--out", help="Path to write merged metrics.tsv", required=True)
    parser.add_argument(
        "--rank",
        help="Path to write merged metrics.tsv",
        default="phylum",
        choices=["phylum", "class"],
    )
    args = parser.parse_args()

    rank_marker_set = args.rank.lower()
    # all_markers -> encompass current rank as well as higher ranks
    # e.g. if class is Gammaproteobacteria -> then phylum markers (proteobacteria) and kingdom markers (bacteria) are also included
    search_str = os.path.join(args.indir, f"*{rank_marker_set}*all_markers.metrics.tsv")
    dfs = []
    for fp in glob.glob(search_str):
        df = pd.read_csv(fp, sep="\t", index_col="cluster")
        lineage_marker_set = (
            os.path.basename(fp)
            .split("_markers_")[0]
            .split(f"_{rank_marker_set}_")[-1]
            .lower()
        )
        df["lineage_marker_set"] = lineage_marker_set
        df["sponge"] = os.path.basename(fp).rsplit(f"_{rank_marker_set}_", 1)[0]
        dfs.append(df)
    df = pd.concat(dfs)

    cols = [col for col in df.columns]
    sponge_col = cols.pop()
    lineage_marker_set_col = cols.pop()
    df["rank_marker_set"] = rank_marker_set
    # Sort columns so sponge, rank and lineage are next to cluster index
    for col in [lineage_marker_set_col, "rank_marker_set", sponge_col]:
        cols.insert(0, col)
    df = df[cols]
    df.to_csv(args.out, sep="\t", index=True, header=True)
    unclustered_count = df.loc["unclustered"].shape[0]
    marker_set_count = df[lineage_marker_set_col].nunique()
    cluster_count = int((df.shape[0] - unclustered_count) / marker_set_count)
    n_sponges = int(len(dfs) / marker_set_count)
    print(
        f"Wrote {cluster_count:,} clusters of {n_sponges:,} LMA sponges to {args.out}"
    )


if __name__ == "__main__":
    main()