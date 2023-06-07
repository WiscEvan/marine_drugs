#!/usr/bin/env python

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Path to humann profiles table", required=True)
    parser.add_argument("--output", help="Path to write munged table", required=True)
    args = parser.parse_args()
    df = pd.read_csv(args.input, sep="\t")
    # First remove Remove RNA-seq filename data and reduce to just sample name to match sponge metadata.
    # R reads '-' in as '.' so we explictly replace here.
    df.columns = df.columns.map(
        lambda c: c.split("_")[0].replace("-", ".").replace("#", "")
    )
    # Now get rid of extra whitespace for '# Pathway' or '# Gene Family'
    # and Add FL2014 to FL20 samples. Need to insert 0 char for single ints as well
    df.columns = df.columns.map(
        lambda c: c.strip()
        if "Pathway" in c or "FL2015" in c or "Gene" in c
        else f"FL2014.{int(c.split('.')[-1])}"
    )


    index_col = [col for col in df.columns if "FL20" not in col][0]
    df.set_index(index_col, inplace=True)
    df.index = df.index.map(lambda x: x.replace("##", "--"))

    if "Pathway" in index_col:
        df = df.transpose()
    df.to_csv(args.output, sep="\t", index=True, header=True)
    print(f"Wrote munged table to {args.output}")


if __name__ == "__main__":
    main()