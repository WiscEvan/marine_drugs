#!/usr/bin/env python

import argparse
import pandas as pd
import glob
import os
from tqdm import tqdm


def read_kegg_counts(fp: str) -> pd.DataFrame:
    with open(fp) as fh:
        kos = []
        for line in fh:
            llist = line.strip().split("\t")
            orf, kegg = (llist[0], llist[1]) if len(llist) >= 2 else (llist[0], pd.NA)
            ko_hit = {"orf": orf, "knum": kegg}
            kos.append(ko_hit)
    df = pd.DataFrame(kos).dropna(subset=["knum"])
    return df.knum.value_counts().to_frame()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        help="Path to output directory of kofamscan with tables containing '*.kofamscan.tsv' in their filename",
        required=True,
    )
    parser.add_argument(
        "--output",
        help="Path to write merged counts of kofamscan results",
        required=False,
        default="ko_counts.tsv",
    )
    args = parser.parse_args()
    table_dirpath = os.path.join(args.input, "*kofamscan.tsv")
    filepaths = glob.glob(table_dirpath)
    counts_dfs = []
    for filepath in tqdm(
        filepaths, total=len(filepaths), desc="Parsing kofamscan results"
    ):
        sponge = os.path.basename(filepath).split(".")[0].replace("_", "-")
        df = read_kegg_counts(filepath).rename(columns={"knum": sponge})
        df.index.name = "knum"
        counts_dfs.append(df)
    ko_df = pd.concat(counts_dfs, axis=1).convert_dtypes().fillna(0)
    ko_df.to_csv(args.output, sep="\t", index=True, header=True)
    print(
        f"Wrote {ko_df.shape[0]:,} knums of {ko_df.shape[1]:,} kofamscan results to {args.output}"
    )


if __name__ == "__main__":
    main()
