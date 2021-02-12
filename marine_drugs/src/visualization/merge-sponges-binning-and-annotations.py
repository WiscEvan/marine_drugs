#!/usr/bin/env python

import argparse
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import SeqUtils


def get_binning(filepath: str, sep: str, sponge: str) -> pd.DataFrame:
    # Now read in binning filepath
    binning = pd.read_csv(filepath, sep=sep, index_col="contig")
    # Select appropriate binning column
    refinement = False
    for col in binning.columns:
        if "refinement_" in col:
            refinement = True
            break
    # Subset the columns down to only one final cluster
    bin_col = binning.columns.tolist()[-1] if refinement else "cluster"
    cols = ["sponge", bin_col]
    binning["sponge"] = sponge
    sponge_master = binning[cols]
    if bin_col != "cluster":
        sponge_master.rename(columns={bin_col: "cluster"}, inplace=True)
    return sponge_master


def get_orfs_info(filepath: str) -> pd.DataFrame:
    orfs = pd.DataFrame(
        [
            {
                "contig": record.id.rsplit("_", 1)[0],
                "orf": record.id,
                "orf_len": len(record.seq),
            }
            for record in SeqIO.parse(filepath, "fasta")
        ]
    )
    # Overwrite orfs df with grouped orfs df
    orfs = (
        orfs.groupby("contig")[["orf", "orf_len"]]
        .agg({"orf": "count", "orf_len": "mean"})
        .rename(columns={"orf": "orf_count", "orf_len": "mean_orf_length"})
    )
    return orfs


def get_contig_stats(assembly: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "contig": record.id,
                "GC": SeqUtils.GC(record.seq),
                "length": len(record.seq),
                "coverage": record.id.rsplit("_", 1)[-1],
            }
            for record in SeqIO.parse(assembly, "fasta")
        ]
    ).set_index("contig")


def aggregate_master_table(
    assemblies_dirpath: str, binning_dirpath: str
) -> pd.DataFrame:
    final_binning_filepaths = os.path.join(
        binning_dirpath, "final_binning_filepaths.tsv"
    )
    df = pd.read_csv(
        final_binning_filepaths, sep="\t", header=None, names=["sponge", "filepath"]
    )
    n_contigs = 0
    masters = []
    for i, row in df.iterrows():
        # Now read in binning filepath
        sep = "," if row.filepath.endswith(".csv") else "\t"
        sponge = row.sponge.replace("-", "_")
        sponge_master = get_binning(filepath=row.filepath, sep=sep, sponge=sponge)
        # Then we will merge other annotations to the final clustering dataframe. e.g. x, y, taxonomy, coverage, GC, num.orfs, etc.
        # Now let's retrieve each contig's number of ORFs and mean orf length
        orf_filepath = os.path.join(binning_dirpath, f"{sponge}.orfs.faa")
        orfs = get_orfs_info(filepath=orf_filepath)
        # Now get contig assembly information
        assembly = os.path.join(assemblies_dirpath, f"{sponge}.filtered.fna")
        contig_stats = get_contig_stats(assembly=assembly)
        # Find kmer coordinates filepath
        kmer_filepath = os.path.join(
            binning_dirpath, f"{sponge}.kmers.bacteria.am_clr.umap.tsv"
        )
        # Find taxonomy annotations
        taxonomy_fp = os.path.join(binning_dirpath, f"{sponge}.taxonomy.tsv")
        annotation_filepaths = [kmer_filepath, taxonomy_fp]
        annotations = [orfs, contig_stats]
        for fp in annotation_filepaths:
            annotation_df = pd.read_csv(fp, sep="\t", index_col="contig")
            annotations.append(annotation_df)
        for annotation in annotations:
            sponge_master = pd.merge(
                sponge_master, annotation, how="left", left_index=True, right_index=True
            )
        n_contigs += sponge_master.shape[0]
        print(
            f"{sponge_master.shape[0]:,} contigs added to master table from {sponge}. (Total: {n_contigs:,})"
        )
        masters.append(sponge_master)
    master = pd.concat(masters)
    assert n_contigs == master.shape[0], "error merging binning files"
    cluster_criterion = master.cluster != "unclustered"
    return master[cluster_criterion]


def fragmentation_metric(df, quality_measure=0.50):
    """Describes the quality of assembled genomes that are fragmented in
    contigs of different length.
    For more information see:
        http://www.metagenomics.wiki/pdf/definition/assembly/n50
    Parameters
    ----------
    quality_measure : 0 < float < 1
        Description of parameter `quality_measure` (the default is .50).
        I.e. default measure is N50, but could use .1 for N10 or .9 for N90
    Returns
    -------
    int
        Minimum contig length to cover `quality_measure` of genome (i.e. percentile contig length)
    """
    target_size = df.length.sum() * quality_measure
    lengths = 0
    for length in df.length.sort_values(ascending=False):
        lengths += length
        if lengths > target_size:
            return length


def load_markers(filepath: str) -> pd.DataFrame:
    # index=contig, cols=[domain sacc,..] (default)
    df = pd.read_csv(filepath, sep="\t", index_col="contig")
    grouped_df = df.groupby("contig")["sacc"]
    return grouped_df.value_counts().unstack()


def get_metabin_stats(bin_df: pd.DataFrame, binning_dirpath: str) -> pd.DataFrame:
    """Retrieve statistics for all clusters recovered from Autometa binning.
    Parameters
    ----------
    bin_df : pd.DataFrame
        Autometa binning table. index=contig, cols=['cluster','length', 'GC', 'coverage', ...]
    markers_fpath : str
        </path/to/{domain}.markers.tsv

    Returns
    -------
    pd.DataFrame
        dataframe consisting of various metabin statistics indexed by cluster.
    """
    stats = []
    markers_dfs = {
        sponge: load_markers(
            os.path.join(binning_dirpath, f"{sponge}.bacteria.markers.tsv")
        )
        for sponge in bin_df.sponge.unique().tolist()
    }
    bin_df = bin_df.convert_dtypes()
    for sponge_cluster, dff in bin_df.groupby(["sponge", "cluster"]):
        sponge, cluster = sponge_cluster
        markers_df = markers_dfs[sponge]
        length_weighted_coverage = np.average(
            a=dff.coverage,
            weights=dff.length / dff.length.sum(),
        )
        length_weighted_gc = np.average(a=dff.GC, weights=dff.length / dff.length.sum())
        # num_expected_markers = markers_df.shape[1]
        # We are hardcoding this, because some markers were not recovered in some of the assemblies.
        num_expected_markers = 139
        pfam_counts = markers_df.loc[markers_df.index.isin(dff.index)].sum()
        if pfam_counts.empty:
            total_markers = 0
            num_single_copy_markers = 0
            num_markers_present = 0
            completeness = pd.NA
            purity = pd.NA
        else:
            total_markers = pfam_counts.sum()
            num_single_copy_markers = pfam_counts[pfam_counts == 1].count()
            num_markers_present = pfam_counts[pfam_counts >= 1].count()
            completeness = num_markers_present / num_expected_markers * 100
            purity = num_single_copy_markers / num_markers_present * 100
        nseqs = dff.shape[0]
        stats.append(
            {
                "sponge": sponge,
                "cluster": cluster,
                "nseqs": nseqs,
                "size (bp)": dff.length.sum(),
                "N90": fragmentation_metric(dff, quality_measure=0.9),
                "N50": fragmentation_metric(dff, quality_measure=0.5),
                "N10": fragmentation_metric(dff, quality_measure=0.1),
                "length_weighted_gc": length_weighted_gc,
                "min_GC": dff.GC.min(),
                "max_GC": dff.GC.max(),
                "std_GC": dff.GC.std(),
                "length_weighted_coverage": length_weighted_coverage,
                "min_coverage": dff.coverage.min(),
                "max_coverage": dff.coverage.max(),
                "std_coverage": dff.coverage.std(),
                "total_orf_count": dff.orf_count.sum(),
                "min_mean_orf_length": dff.mean_orf_length.min(),
                "max_mean_orf_length": dff.mean_orf_length.max(),
                "std_mean_orf_length": dff.mean_orf_length.std(),
                "completeness": completeness,
                "purity": purity,
                "num_total_markers": total_markers,
                f"num_unique_markers (expected {num_expected_markers})": num_markers_present,
                "num_single_copy_markers": num_single_copy_markers,
            }
        )
    return pd.DataFrame(stats).set_index(["sponge", "cluster"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("contig_out", help="path to write output master table.")
    parser.add_argument("cluster_out", help="path to write output master table.")
    parser.add_argument("assemblies", help="directory path filtered assemblies.")
    parser.add_argument(
        "annotations", help="directory path interim binning annotations."
    )
    args = parser.parse_args()
    if os.path.exists(args.contig_out) and os.stat(args.contig_out):
        contig_master = pd.read_csv(args.contig_out, sep="\t", index_col="contig")
    else:
        contig_master = aggregate_master_table(
            assemblies_dirpath=args.assemblies,
            binning_dirpath=args.annotations,
        )
    cluster_master = get_metabin_stats(contig_master, args.annotations)
    cluster_master.to_csv(args.cluster_out, sep="\t", index=True, header=True)
    contig_master.to_csv(args.contig_out, sep="\t", index=True, header=True)


if __name__ == "__main__":
    main()