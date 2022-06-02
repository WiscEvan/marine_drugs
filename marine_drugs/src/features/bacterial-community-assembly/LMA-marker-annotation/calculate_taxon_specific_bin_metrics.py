#!/usr/bin/env python

import argparse
import os
import pandas as pd


def read_hmmscan(infpath, outfpath):
    col_indices = [0, 1, 2, 5]
    df = pd.read_csv(infpath, sep="\s+", usecols=col_indices, header=None, comment="#")
    hmmtab_header = ["sname", "sacc", "orf", "score"]
    columns = {i: k for i, k in zip(col_indices, hmmtab_header)}
    df = df.rename(columns=columns).dropna()
    # Get contigs from annotated ORFs to perform bin completeness/purity calculations
    # Convert ORFs to Contigs
    df["contig"] = df.orf.map(lambda x: x.rsplit("_", 1)[0])
    df.to_csv(outfpath, sep="\t", index=False, header=True)
    grouped_df = df.groupby("contig")["sacc"]
    return grouped_df.value_counts().unstack().convert_dtypes()


def read_binning(fpath):
    sep = "\t" if fpath.endswith(".tsv") else ","
    bin_df = pd.read_csv(fpath, sep=sep, index_col="contig").fillna(
        axis=1, method="ffill"
    )
    if "gc_content" in bin_df.columns and "cluster" in bin_df.columns:
        bin_col = "cluster"
    elif "gc_content" in bin_df.columns and "recruited_cluster" in bin_df.columns:
        bin_col = "recruited_cluster"
    else:
        bin_col = bin_df.columns[-1]
    # Retrieve final refinement column and rename to cluster
    return bin_df[[bin_col]].rename(columns={bin_col: "cluster"})


def get_bin_metrics(df, markers_df):
    metrics = []
    # markers_df.shape :=> (rows=contigs, cols=markers)
    expected_number = markers_df.shape[1]
    for cluster, dff in df.groupby("cluster"):

        markers = markers_df[markers_df.index.isin(dff.index)].sum()

        ## Get copy number marker counts
        copy_num_counts = {}
        copy_nums = [1, 2, 3, 4, 5]
        for i, copy_num in enumerate(copy_nums):
            # Check if last in the list of copy_nums
            if i + 1 == len(copy_nums):
                copy_num_filter = markers >= copy_num
            else:
                copy_num_filter = markers == copy_num
            copy_num_marker_count = markers[copy_num_filter].count()
            copy_num_pfams = ",".join(
                sorted(
                    marker
                    for marker in markers[copy_num_filter].index
                    if "PF" in marker
                )
            )
            copy_num_tigrfams = ",".join(
                sorted(
                    marker
                    for marker in markers[copy_num_filter].index
                    if "TIGR" in marker
                )
            )
            # At the last copy number indicate that this operation was greater than or equal
            if i + 1 == len(copy_nums):
                copy_num_counts.update(
                    {
                        f"{copy_num}+": copy_num_marker_count,
                        f"({copy_num}+)-copy PFAMS": copy_num_pfams,
                        f"({copy_num}+)-copy TIGRFAMS": copy_num_tigrfams,
                    }
                )
            # Otherwise put the copy number
            else:
                copy_num_counts.update(
                    {
                        copy_num: copy_num_marker_count,
                        f"{copy_num}-copy PFAMS": copy_num_pfams,
                        f"{copy_num}-copy TIGRFAMS": copy_num_tigrfams,
                    }
                )

        ## Calculate Bin Metrics
        is_present = markers.ge(1)
        present_marker_count = markers[is_present].count()

        # Calculate Completeness
        completeness = present_marker_count / expected_number * 100
        completeness = round(completeness, 2)

        # Calculate purity
        # NOTE: Protect from divide by zero
        single_copy_marker_count = copy_num_counts[1]
        if present_marker_count == 0:
            purity = pd.NA
        else:
            purity = single_copy_marker_count / present_marker_count * 100
            purity = round(purity, 2)

        # Get total marker sum
        marker_sum = markers.sum()

        # get all present markers PFAMS and TIGRFAMS
        pfams = ",".join(
            sorted(marker for marker in markers[is_present].index if "PF" in marker)
        )
        tigrfams = ",".join(
            sorted(marker for marker in markers[is_present].index if "TIGR" in marker)
        )

        # Get multi-copy marker sum
        multi_copy_marker_sum = 0
        for copy_num in copy_nums:
            if copy_num == 1:
                continue
            if copy_num in copy_num_counts:
                multi_copy_marker_count = copy_num_counts[copy_num]
            elif f"{copy_num}+" in copy_num_counts:
                multi_copy_marker_count = copy_num_counts[f"{copy_num}+"]
            else:
                raise ValueError(
                    f"copy num ({copy_num}) not in copy_num_counts keys: {copy_num_counts.keys()}"
                )
            if multi_copy_marker_count:
                # Account for 1 of each multi-copy marker in copy nums except single-copy markers
                multi_copy_marker_count - 1
            multi_copy_marker_sum += multi_copy_marker_count

        metrics.append(
            {
                "cluster": cluster,
                "completeness": completeness,
                "purity": purity,
                "present_marker_count": present_marker_count,
                "single_copy_marker_count": single_copy_marker_count,
                "multi_copy_marker_sum": multi_copy_marker_sum,
                "marker_sum": marker_sum,
                "PFAMS": pfams,
                "TIGRFAMS": tigrfams,
                **copy_num_counts,
            }
        )
    metrics_df = pd.DataFrame(metrics).set_index("cluster")
    fam_cols = [
        col
        for col in metrics_df.columns
        if isinstance(col, str) and ("PFAM" in col or "TIGRFAM" in col)
    ]
    metric_cols = [
        "completeness",
        "purity",
        "present_marker_count",
        "single_copy_marker_count",
        "multi_copy_marker_sum",
        "marker_sum",
    ]
    mult_cols = [
        col for col in metrics_df.columns if col in copy_nums or col.endswith("+")
    ]
    outcols = mult_cols + metric_cols + fam_cols
    return metrics_df[outcols].convert_dtypes()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--hmmscan",
        help="path to hmmscan.tsv output on taxon-specific marker set",
        required=True,
    )
    parser.add_argument("--binning", help="path to binning.tsv table", required=True)
    parser.add_argument(
        "--output-markers",
        help="path to write markers.tsv output of taxon-specific marker set",
        required=True,
    )
    parser.add_argument(
        "--output-metrics",
        help="path to write binning_metrics.tsv",
        required=True,
    )
    # $HOME/marine_drugs/marine_drugs/data/interim/binning/final_binning_filepaths.tsv
    args = parser.parse_args()

    # "/home/evan/marine_drugs/marine_drugs/data/interim/LMA-marker-annotation/FL2015_4_phylum_Proteobacteria_markers_all_markers.hmmscan.tsv"
    markers_df = read_hmmscan(infpath=args.hmmscan, outfpath=args.output_markers)

    # Read binning to assess metrics
    # binning = "/home/evan/marine_drugs/marine_drugs/data/interim/refined/first_round_of_refinements/before_sams_notes/FL2015_4.bacteria.binning.umap.refined.csv"
    bin_df = read_binning(args.binning)

    metrics_df = get_bin_metrics(df=bin_df, markers_df=markers_df)

    metrics_df.to_csv(args.output_metrics, sep="\t", index=True, header=True)
    print(f"Wrote {metrics_df.shape[0]:,} clusters' metrics to {args.output_metrics}")


if __name__ == "__main__":
    main()