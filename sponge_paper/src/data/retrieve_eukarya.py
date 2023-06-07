#!/usr/bin/env python

# Retrieve numbers of contigs broken down by superkingdom for each provided taxonomy.tsv

import argparse
import os
import glob
import pandas as pd

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(
        description="Retrieve contig counts of superkingdoms for each provided taxonomy.tsv file"
    )
    parser.add_argument(
        "--input",
        help="directory containing files with *.taxonomy.tsv extension",
        required=True,
    )
    parser.add_argument(
        "--fastas",
        help="directory containing fasta files to write superkingdoms."
        " Will write each superkingdom into its own sponge directory",
        required=True,
    )
    parser.add_argument(
        "--output",
        help="Directory to write output table and kingdom fastas.",
        required=True,
    )
    args = parser.parse_args()
    search_string = os.path.join(args.input, "*.taxonomy.tsv")
    filepaths = glob.glob(search_string)
    print(f"aggregating {len(filepaths)} taxonomy files.")

    if not os.path.isdir(args.output):
        os.makedirs(args.output)
        print(f"{args.output} does not exists. Created...")

    taxonomies = []
    for filepath in filepaths:
        filename = os.path.basename(filepath)
        sponge = filename.split(".taxonomy.tsv")[0]
        df = pd.read_csv(filepath, sep="\t", index_col="contig")
        taxonomy = {"sponge": sponge}
        taxonomy.update(df.groupby("superkingdom")["superkingdom"].count().to_dict())
        taxonomies.append(taxonomy)
        if args.fastas:
            try:
                fasta = glob.glob(os.path.join(args.fastas, f"{sponge}.fasta"))[0]
                records = [record for record in SeqIO.parse(fasta, "fasta")]
                print(f"{len(records):,} contigs in {sponge}.fasta")
                for kingdom, dff in df.groupby("superkingdom"):
                    outdir = os.path.join(args.output, sponge)
                    if not os.path.isdir(outdir):
                        os.makedirs(outdir)
                    output = os.path.join(outdir, f"{kingdom}.fna")
                    kingdom_contigs = set(dff.index.tolist())
                    kingdom_records = [
                        record for record in records if record.id in kingdom_contigs
                    ]
                    num_written_records = SeqIO.write(kingdom_records, output, "fasta")
                    print(f"Wrote {num_written_records:,} to {output}")
            except IndexError:
                print(f"failed to locate {sponge} fasta file. skipping...")
                continue

    df = pd.DataFrame(taxonomies)
    output_table = os.path.join(args.output, "sponges_superkingdom_counts.tsv")
    df.to_csv(output_table, sep="\t", index=False, header=True)
    print(f"Wrote sponge superkingdoms contig counts to {output_table}")


if __name__ == "__main__":
    main()