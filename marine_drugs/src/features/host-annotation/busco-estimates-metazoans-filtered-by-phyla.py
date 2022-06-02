#!/usr/bin/env python

# First we are going to perform a BUSCO annotation on the entire metagenome
#  using the parameter: `busco --lineage metazoa ...`
# Then we filter *contigs* containing metazoan markers that have been assigned a taxonomy by Autometa (taxonomy.tsv)
# Note: Phyla to filter based on busco --list-datasets: arthropoda, mollusca, nematoda, vertebrata 
# Above phyla are all under metazoa lineage and we know our sponge should *not* be assigned these taxa
# After we have filtered out the contigs that match any of the above phyla, we re-calculate our BUSCO values from full_table.tsv
# all paths should match: "**/run_metazoa_odb10/full_table.tsv"

import argparse
import os
import glob
import pandas as pd

HOME = os.path.expanduser("~")
DEFAULT_DIR = os.path.join(HOME, "marine_drugs/marine_drugs/data/interim/host-annotation/busco-marker-assessment")
DEFAULT_TAXA_DIR = os.path.join(HOME, "marine_drugs/marine_drugs/data/interim/binning")
DEFAULT_SPONGE_METADATA = os.path.join(HOME, "marine_drugs/marine_drugs/data/raw/sponge_metadata.tsv")

def filter_metazoans(full_table: str, taxonomy: str, unwanted_phyla: list =["arthropoda", "mollusca", "nematoda", "vertebrata"]) -> pd.DataFrame:
    """
    Filter out contigs that are hitting as:
    
    1. arthropoda
    2. mollusca
    3. nematoda
    4. vertebrata
    
    Note: These phyla are selected because busco --list-datasets 
    has them listed under the metazoan markers dataset AND
    we know the contigs pertaining to our sponges should not be assigned these phyla.
    """

    # Busco id	Status	Sequence	Gene Start	Gene End	Strand	Score	Length
    names = ["Busco id","Status","Sequence","Gene Start","Gene End","Strand","Score","Length"]
    df = pd.read_csv(full_table, sep='\t', index_col='Sequence', comment='#', names=names)
    taxa = pd.read_csv(taxonomy, sep='\t', index_col='contig')
    # Now we grab all contigs that hit to the above phyla
    phyla_contigs = set(taxa[taxa.phylum.isin(unwanted_phyla)].index.tolist())
    print(f"{len(phyla_contigs):,} contigs hitting to unwanted phyla")
    # Use the phyla contigs to filter our BUSCO table.
    # First remove any rows that do not have a contig designation (any BUSCO markers that are identified as missing)
    num_prefilter_markers = df.dropna(axis=0).shape[0]
    df = df[~df.index.isin(phyla_contigs)]
    num_postfilter_markers = df.dropna(axis=0).shape[0]
    num_markers_filtered = num_prefilter_markers - num_postfilter_markers
    print(f"{num_markers_filtered} markers removed b/c they were assigned unwanted phyla")
    return df

def calculate_busco_metrics(df: pd.DataFrame) -> dict:
    # We need to check if any 'duplicates' were removed from the taxonomic filter.
    # If no, we can continue on business as usual.
    # Otherwise, we will need to reassign the status to complete rather than duplicate
    # All busco ids under 'duplicated' status could have more than one respective pair, 
    # Really the annotation should be multi-copy rather than single-copy
    # So we are just going to reassign any duplicated Busco ids that are 'duplicated' to duplicated (as in multi-copy)
    duplicate_df = df[df.Status == "Duplicated"]
    is_duplicated = duplicate_df.duplicated(subset="Busco id", keep=False)
    new_single_copies = duplicate_df[~is_duplicated]
    # Now we need to reassign the original df with Status == "Complete" instead of Duplicated
    df.loc[new_single_copies.index, "Status"] = "Complete"
    df.loc[df.duplicated(subset="Busco id", keep=False), "Status"] = "Duplicated"   
    # Now that we have resolved duplicates, we should be able to proceed for all metrics.
    counts = df.Status.value_counts()
    duplicates = df[df.Status == "Duplicated"]["Busco id"].nunique()
    complete_single_copy = counts.Complete
    complete = complete_single_copy + duplicates
    num_metazoan_markers = 954
    return {
        "Complete BUSCOs (C)": complete,
        "Percentage Complete": complete / num_metazoan_markers * 100,
        "Complete and single-copy BUSCOs (S)": complete_single_copy,
        "Complete and duplicated BUSCOs (D)": duplicates,
        "Fragmented BUSCOs (F)": counts.Fragmented,
        "Missing BUSCOs (M)": counts.Missing,
        "Total BUSCO groups searched": num_metazoan_markers,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--busco", help="busco marker assessment directory", default=DEFAULT_DIR)
    parser.add_argument("--taxa", help="Genome marker assessment directory", default=DEFAULT_TAXA_DIR)
    parser.add_argument("--sponge-metadata", help="Path to sponge_metadata.tsv", default=DEFAULT_SPONGE_METADATA)
    parser.add_argument("--output", help="Output BUSCO metrics table path", default="busco_metrics.tsv")
    args = parser.parse_args()
    
    # fullpath $HOME/marine_drugs/marine_drugs/data/interim/host-annotation/busco-marker-assessment/FL2015_8_metagenome_busco_genome_metazoa_odb10/run_metazoa_odb10/full_table.tsv
    # we will need to find sample from two directories up (base output directory of busco):
    # e.g. FL2015_8_metagenome_busco_genome_metazoa_odb10/run_metazoa_odb10/full_table.tsv
    search_string = os.path.join(args.busco, "FL*_metagenome*", "run_metazoa_odb10", "full_table.tsv")
    full_tables = [fpath for fpath in glob.glob(search_string, recursive=True)]
    metrics = []
    lineage = "metazoa_odb10"
    for full_table in full_tables:
        # We know we will always have run_metazoa_obd10/full_table.tsv
        sample = os.path.basename(os.path.dirname(os.path.dirname(os.path.abspath(full_table)))).split("_metagenome")[0]
        taxonomy = os.path.join(args.taxa, f"{sample}.taxonomy.tsv")
        filtered_df = filter_metazoans(full_table=full_table, taxonomy=taxonomy)
        metric = calculate_busco_metrics(filtered_df)
        metric.update({"sample": sample, "lineage": lineage})
        metrics.append(metric)
    
    master = pd.DataFrame(metrics).set_index("sample")
    metadata = pd.read_csv(args.sponge_metadata, sep='\t', index_col="Sponge specimen")
    metadata.index = metadata.index.map(lambda sample: sample.replace("-", "_"))
    master = pd.merge(master, metadata, how='left', left_index=True, right_index=True)
    master.sort_values("Percentage Complete", ascending=False, inplace=True)
    master.to_csv(args.output, sep='\t', header=True, index=True)
    print(f"Wrote {master.shape[0]:,} BUSCO summaries to {args.output}")

if __name__ == "__main__":
    main()