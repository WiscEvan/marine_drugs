#!/usr/bin/env python

# Here are a few notes from the binning that I can't do using Automappa:
# FL2015_37: Split bin3 by class: 3a = alphaproteos, 3b = unclassified (This is true for both bhtsne and umap)
# FL2015_42: Split bin8 into 2 bins at coverage 30X.
# FL2015_43: No additional fiddling necessary
# FL2015_44: Bin10 - remove all contigs with cov >20X
# FL2015_5: No additional fiddling necessary

# These files are now in /home/evan/sponge_paper/sponge_paper/data/interim/refined

import os
import pandas as pd

HOME = os.path.expanduser("~")
interim_binning = os.path.join(
    HOME, "sponge_paper", "sponge_paper", "data", "interim", "binning"
)
refinements_from_sam = os.path.join(
    HOME,
    "sponge_paper",
    "sponge_paper",
    "data",
    "interim",
    "refined",
    "before_sams_notes",
)
outdir = os.path.join(
    HOME,
    "sponge_paper",
    "sponge_paper",
    "data",
    "interim",
    "refined",
    "after_sams_notes",
)
if not os.path.isdir(outdir):
    os.makedirs(outdir)


def refine_FL2015_37():
    # FL2015_37: Split refinement_3 by class: 3a = alphaproteos, 3b = unclassified (This is true for both bhtsne and umap)
    sponge = "FL2015_37"
    filepath = os.path.join(
        refinements_from_sam, f"{sponge}.bacteria.binning.umap.refined.csv"
    )
    master = pd.read_csv(filepath, index_col="contig")
    taxonomy = os.path.join(interim_binning, f"{sponge}.taxonomy.tsv")
    dff = pd.read_csv(taxonomy, sep="\t", index_col="contig")
    mdf = pd.merge(master, dff, how="left", left_index=True, right_index=True)
    # Grab latest refinement column.
    bin_col = [col for col in mdf.columns if "refinement" in col][-1]
    cols = [bin_col, "class"]
    mdf = mdf[cols]
    # Now we retrieve the latest refinement and retrieve refinement_3
    refinement = mdf[mdf[bin_col] == "refinement_3"]
    print(
        f"Splitting refinement_3's {refinement.shape[0]:,} contigs into two bins of alphaproteos and non-alphaproteos."
    )
    # Now we split this bin by contigs that are/are not classified as alphaproteobacteria.
    non_alphaproteos = refinement[refinement["class"] != "alphaproteobacteria"]
    alphaproteobacteria = refinement[refinement["class"] == "alphaproteobacteria"]
    # Apply final results to table and write to output directory.
    refinement_num = int(bin_col.split("_")[-1])
    # Here we are removing contigs above the 20x coverage from the refinement_10 bin.
    for refined_df in [non_alphaproteos, alphaproteobacteria]:
        print(f"Adding {refined_df.shape[0]:,} contigs into their own bin")
        refinement_num += 1
        refinement_col = f"refinement_{refinement_num}"
        master.loc[refined_df.index, refinement_col] = refinement_col
        master.fillna(axis="columns", method="ffill", inplace=True)
    filename = os.path.basename(filepath)
    out = os.path.join(outdir, filename)
    master.to_csv(out, sep=",", index=True, header=True)
    print(f"Wrote refinement to {out}")


def refine_FL2015_42():
    # FL2015_42: Split refinement_8 into 2 bins at coverage 30X.
    sponge = "FL2015_42"
    filepath = os.path.join(
        refinements_from_sam, f"{sponge}.bacteria.binning.umap.refined.csv"
    )
    master = pd.read_csv(filepath, index_col="contig")
    coverages = os.path.join(interim_binning, f"{sponge}.coverages.tsv")
    dff = pd.read_csv(coverages, sep="\t", index_col="contig")
    mdf = pd.merge(master, dff, how="left", left_index=True, right_index=True)
    # From inspection of the file latest_refinement is refinement_17
    # bin_col = "refinement_17"
    # Grab latest refinement column.
    bin_col = [col for col in mdf.columns if "refinement" in col][-1]
    cols = [bin_col, "coverage"]
    mdf = mdf[cols]
    # Now we retrieve the latest refinement and retrieve refinement_8
    refinement = mdf[mdf[bin_col] == "refinement_8"]
    print(
        f"Splitting refinement_8 {refinement.shape[0]:,} contigs into two bins of above and below 30x coverage."
    )
    # Now we retrieve contigs above 30x coverage and below 30x coverage.
    above_cov = refinement[refinement.coverage >= 30.0]
    below_cov = refinement[refinement.coverage < 30.0]
    # Apply final results to table and write to output directory.
    refinement_num = int(bin_col.split("_")[-1])
    for refined_df in [above_cov, below_cov]:
        print(f"Adding {refined_df.shape[0]:,} contigs into their own bin")
        refinement_num += 1
        refinement_col = f"refinement_{refinement_num}"
        master.loc[refined_df.index, refinement_col] = refinement_col
        master.fillna(axis="columns", method="ffill", inplace=True)

    filename = os.path.basename(filepath)
    out = os.path.join(outdir, filename)
    master.to_csv(out, sep=",", index=True, header=True)
    print(f"Wrote refinement to {out}")


def refine_FL2015_44():
    # FL2015_44: refinement_10 - remove all contigs with cov >20X
    sponge = "FL2015_44"
    filepath = os.path.join(
        refinements_from_sam, f"{sponge}.bacteria.binning.umap.refined.csv"
    )
    master = pd.read_csv(filepath, index_col="contig")
    coverages = os.path.join(interim_binning, f"{sponge}.coverages.tsv")
    dff = pd.read_csv(coverages, sep="\t", index_col="contig")
    mdf = pd.merge(master, dff, how="left", left_index=True, right_index=True)
    # Grab latest refinement column.
    bin_col = [col for col in mdf.columns if "refinement" in col][-1]
    cols = [bin_col, "coverage"]
    mdf = mdf[cols]
    # Now we retrieve the latest refinement and retrieve refinement_10
    refinement = mdf[mdf[bin_col] == "refinement_10"]
    print(
        f"Checking refinement_10's {refinement.shape[0]:,} contigs above 20x coverage"
    )
    above_cov = refinement[refinement.coverage > 20.0]
    # Apply final results to table and write to output directory.
    refinement_num = int(bin_col.split("_")[-1])
    refinement_num += 1
    refinement_col = f"refinement_{refinement_num}"
    # Here we are removing contigs above the 20x coverage from the refinement_10 bin.
    master.loc[above_cov.index, refinement_col] = "unclustered"
    print(f"unclustered {above_cov.shape[0]} contigs from refinement_10")
    # The other contigs in refinement_10 will be forward filled to the next latest refinement column.
    master.fillna(axis="columns", method="ffill", inplace=True)
    filename = os.path.basename(filepath)
    out = os.path.join(outdir, filename)
    master.to_csv(out, sep=",", index=True, header=True)
    print(f"Wrote refinement to {out}")


if __name__ == "__main__":

    refine_FL2015_37()
    refine_FL2015_42()
    refine_FL2015_44()