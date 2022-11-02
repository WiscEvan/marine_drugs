#!/usr/bin/env python3
"""
Takes arwen output and converts to genbank format
"""

import argparse
import os
import re

# See https://caretdashcaret.com/2016/06/15/how-to-create-genbank-files-with-biopython/
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation


from datetime import datetime

TODAY = datetime.today().strftime("%d-%b-%Y").upper()


def mtdna_length(arwen_fpath):
    seq_len_desc = "nucleotides in sequence"
    with open(arwen_fpath) as fh:
        for line in fh:
            if seq_len_desc in line:
                return int(line.split(seq_len_desc)[0])


def get_trna(arwen_fpath):
    """Retrieves trna as Bio.SeqRecord...
    ReturnType: list
    Returns: [trna_feature,...]

    [SeqFeature(start,end,type,qualifiers)]
    feature = SeqFeature(start=x,end=y),type='tRNA'
         tRNA            2492..2555
                         /product="tRNA-Asn"
    """
    features = []
    trna_count = 0
    feature = None
    seq = None
    with open(arwen_fpath) as fh:
        getseq = False
        for line in fh:
            cond1 = "mtRNA" in line
            cond2 = "Overlap" not in line
            cond3 = "Primary" not in line
            if cond1 and cond2 and cond3:
                trna_count += 1
                mtrna = line.replace(" ", "-").lstrip("-")
                line = next(fh)
                line = next(fh)
                m = re.search("\[(\d+),(\d+)\]", line)
                start, end = map(int, [m.group(1), m.group(2)])
                start -= 1
                strand = line.strip().split()[1]
                strand = -1 if "c" in strand else 1
                feature = SeqFeature(
                    FeatureLocation(start=start, end=end, strand=strand),
                    type="tRNA",
                    qualifiers={"product": mtrna, "note": "arwen"},
                )
                features.append(feature)
    return features


def main():
    parser = argparse.ArgumentParser("Converts arwen output to genbank")
    parser.add_argument(
        "--out", help="</path/to/output/genbank>", default="arwen_output.gbk"
    )
    parser.add_argument("arwen", help="output from arwen")
    parser.add_argument("mtdna", help="mtDNA record")
    args = parser.parse_args()
    features = get_trna(args.arwen)
    record = [rec for rec in SeqIO.parse(args.mtdna, "fasta")][0]
    for feature in features:
        record.features.append(feature)
    record.annotations["molecule_type"] = "DNA"
    SeqIO.write(record, args.out, "genbank")
    print(f"Written: {args.out}")


if __name__ == "__main__":
    main()
