#!/usr/bin/env python

import argparse
import os
import glob
import pandas as pd

HOME = os.path.expanduser("~")
DEFAULT_DIR = os.path.join(HOME, "sponge_paper/sponge_paper/data/interim/host-annotation/busco-marker-assessment/")
DEFAULT_SPONGE_METADATA = os.path.join(HOME, "sponge_paper/sponge_paper/data/raw/sponge_metadata.tsv")

def retrieve_summaries(dirpath):
    metrics = []
    # No recursive here b/c I think there are other summaries in subdirectories that we do *NOT* want to grab.
    # short_summary.generic.eukaryota_odb10.GCA_013339895.1_Emu_genome_v1_genomic_busco_genome_auto_lineage_euk.txt
    search_string = os.path.join(dirpath, "**", "short_summary.*.txt")
    summaries = glob.glob(search_string, recursive=True)
    for fpath in summaries:
        # short_summary.specific.metazoa_odb10.FL2015_9_busco_genome_metazoa_odb10.txt
        # short_summary.specific.metazoa_odb10.GCF_000090795.1_v1.0_genomic_busco_genome_metazoa_odb10.txt
        if "short_summary.specific." in os.path.basename(fpath):
            lineage, sample = os.path.basename(fpath).split("short_summary.specific.")[1].strip(".txt").split(".", 1)
        elif "short_summary.generic." in os.path.basename(fpath):
            lineage, sample = os.path.basename(fpath).split("short_summary.generic.")[1].strip(".txt").split(".", 1)
        else:
            raise LookupError(f"{os.path.basename(fpath)} does not appear to be a file that should be parsed!")
        # Now we want to pull the busco mode out from the file naming, e.g. 'genome' or 'proteins'
        # ./FL2015_34_busco_genome_auto_lineage_euk
        # ./FL2015_44_busco_proteins_auto_lineage_euk
        # ./FL2015_42_metagenome_busco_genome_metazoa_odb10
        dirname = os.path.basename(os.path.dirname(fpath))
        if "_metagenome_busco_genome_" in dirname:
            sample = sample.split("_busco_genome_")[0]
            busco_mode = "metagenome"
        elif "_busco_genome_" in dirname:
            sample = sample.split("_busco_genome_")[0]
            busco_mode = "genome"
        elif "_busco_proteins_" in dirname:
            sample = sample.split("_busco_proteins_")[0]
            busco_mode = "proteins"
        elif "_markers" in dirname:
            sample = sample.split("_markers")[0]
            busco_mode = "proteins"
        else:
            raise LookupError(f"{os.path.basename(fpath)} does not appear to be a file that should be parsed!")

        with open(fpath) as fh:
            for line in fh:
                line = line.strip()
                if "Complete BUSCOs" in line:
                    C = line.split()[0]
                if "Complete and single-copy BUSCOs" in line:
                    S = line.split()[0]
                if "Complete and duplicated BUSCOs" in line:
                    D = line.split()[0]
                if "Fragmented BUSCOs" in line:
                    F = line.split()[0]
                if "Missing BUSCOs" in line:
                    M = line.split()[0]
                if "Total BUSCO groups searched" in line:
                    total_searched = line.split()[0]
        metric = {
            "sample": sample,
            "busco mode": busco_mode,
            "lineage": lineage,
            "Complete BUSCOs (C)": C,
            "Complete and single-copy BUSCOs (S)": S,
            "Complete and duplicated BUSCOs (D)": D,
            "Fragmented BUSCOs (F)": F,
            "Missing BUSCOs (M)": M,
            "Total BUSCO groups searched": total_searched,
        }
        metrics.append(metric)

    df = pd.DataFrame(metrics)
    columns = df.columns.tolist()
    df["Percentage Complete"] =  df["Complete BUSCOs (C)"].astype(int) / df["Total BUSCO groups searched"].astype(int) * 100
    columns.insert(2, "Percentage Complete")
    df = df[columns]
    df.set_index("sample", inplace=True)
    df.drop_duplicates(inplace=True)
    df.sort_values(["lineage", "Percentage Complete"], ascending=[True, False], inplace=True)
    return df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Genome marker assessment directory", default=DEFAULT_DIR)
    parser.add_argument("--output", help="Output BUSCO metrics table path", default="busco_metrics.tsv")
    parser.add_argument("--sponge-metadata", help="Path to sponge_metadata.tsv", default=DEFAULT_SPONGE_METADATA)
    args = parser.parse_args()
    
    df = retrieve_summaries(dirpath=args.input)
    metadata = pd.read_csv(args.sponge_metadata, sep='\t', index_col="Sponge specimen")
    metadata.index = metadata.index.map(lambda sample: sample.replace("-", "_"))
    master = pd.merge(df, metadata, how='left', left_index=True, right_index=True)
    master.index.name = 'sample'
    master.sort_values(["lineage", "Percentage Complete", "busco mode"], ascending=[True, False, False], inplace=True)
    master.index = master.index.map(lambda x: x.replace("_genomic","").replace("_protein","") if "GCF_000090795.1_v1.0" in x else x)

    master.to_csv(args.output, sep='\t', header=True, index=True)
    print(f"Wrote all BUSCO summaries to {args.output}")

    out = os.path.join(os.path.dirname(args.output), "busco_metazoa_odb10_metrics.tsv")
    metazoa = master[master.lineage == "metazoa_odb10"]
    metazoa.to_csv(out, sep='\t', header=True, index=True)
    print(f"Wrote {metazoa.shape[0]:,} metazoa BUSCO summaries to {out}")

if __name__ == "__main__":
    main()