#!/usr/bin/env python

import argparse
import pandas as pd
import glob
import os
from tqdm import tqdm


def read_kofamscan(fp: str) -> pd.DataFrame:
    with open(fp) as fh:
        kos = []
        for line in fh:
            llist = line.strip().split("\t")
            orf, kegg = (llist[0], llist[1]) if len(llist) >= 2 else (llist[0], pd.NA)
            ko_hit = {"orf": orf, "knum": kegg}
            kos.append(ko_hit)
    return pd.DataFrame(kos).dropna(subset=["knum"])


def read_blastp_amphimedon_hits(fp: str) -> pd.DataFrame:
    df = pd.read_csv(fp, sep="\t")
    df["sponge"] = df.file.map(
        lambda x: os.path.basename(x).split(".")[0].replace("_", "-")
    )
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--indir",
        help="Path to output directory of kofamscan with tables containing '*.kofamscan.tsv' in their filename",
        required=True,
    )
    parser.add_argument(
        "--blastp-hits",
        help="Path to sponge_proteins_amphimedon_blastp_hits.tsv",
        required=True,
    )
    parser.add_argument(
        "--output",
        help="Path to write counts of kofamscan results subset by hits to amphimedon queenslandica",
        required=False,
        default="amphimedon_queenslandica_hits_ko_counts.tsv",
    )
    args = parser.parse_args()
    table_dirpath = os.path.join(args.indir, "*kofamscan.tsv")
    dff = read_blastp_amphimedon_hits(args.blastp_hits)
    filepaths = glob.glob(table_dirpath)
    counts_dfs = []
    for filepath in tqdm(
        filepaths, total=len(filepaths), desc="Parsing kofamscan results"
    ):
        sponge = os.path.basename(filepath).split(".")[0].replace("_", "-")
        df = read_kofamscan(filepath)
        if sponge in dff.sponge.unique():
            # subset by amphimedon queenslandica hits
            hits_df = dff.loc[dff.sponge.eq(sponge)]
            df = df.loc[df.orf.isin(hits_df.qseqid)].copy()
        counts_df = df.knum.value_counts().to_frame().rename(columns={"knum": sponge})
        counts_df.index.name = "knum"
        counts_dfs.append(counts_df)

    ko_df = pd.concat(counts_dfs, axis=1).convert_dtypes().fillna(0)
    ko_df.to_csv(args.output, sep="\t", index=True, header=True)
    print(
        f"Wrote {ko_df.shape[0]:,} knums of {ko_df.shape[1]:,} kofamscan results to {args.output}"
    )


if __name__ == "__main__":
    main()
