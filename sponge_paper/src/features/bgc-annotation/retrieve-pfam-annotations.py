#!/usr/bin/env python

import argparse
import pandas as pd

from Bio import SeqIO


def retrieve_pfam_annotations(filepath:str)->pd.DataFrame:
    print(f"Parsing {filepath} for PFAM_domain feature")
    pfams = []
    for record in SeqIO.parse(filepath, "genbank"):
        for feature in record.features:
            if feature.type == "PFAM_domain":
                desc = feature.qualifiers.get("description", ["???"])[0]
                pfam_db = feature.qualifiers.get("database", ["???"])[0]
                db_xref = feature.qualifiers.get("db_xref", ["???"])[0]
                domain_id = feature.qualifiers.get("domain_id", ["???"])[0]
                label = feature.qualifiers.get("label", ["???"])[0]
                locus_tag = feature.qualifiers.get("locus_tag", ["???"])[0]
                prot_start = feature.qualifiers.get("protein_start", ["???"])[0]
                prot_end = feature.qualifiers.get("protein_end", ["???"])[0]
                translation = feature.qualifiers.get("translation", ["???"])[0]

                pfams.append({
                    "contig":record.description,
                    "PFAM": db_xref,
                    "description":desc,
                    "PFAM database": pfam_db,
                    "antismash contig":record.id,
                    "antismash domain ID":domain_id,
                    "antismash label":label,
                    "antismash locus_tag": locus_tag,
                    "protein start":prot_start,
                    "protein end": prot_end,
                    "translation":translation,
                    "feature start": str(feature.location.start),
                    "feature end": str(feature.location.end),
                    "feature strand": str(feature.location.strand),
                })
    return pd.DataFrame(pfams)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Path to main antismash genbank output file",required=True)
    parser.add_argument("--output", help="Path to write pfam annotations table",required=True)
    args = parser.parse_args()

    df = retrieve_pfam_annotations(args.input)

    df.to_csv(args.output, sep='\t', header=True, index=False)
    print(f"Wrote {df.shape[0]:,} PFAM domain annotations to {args.output}")

if __name__ == "__main__":
    main()