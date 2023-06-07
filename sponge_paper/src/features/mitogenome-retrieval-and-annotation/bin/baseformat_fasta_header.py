#!/usr/bin/env python

import argparse
from typing import List
from Bio import SeqIO


def remove_descriptions(fpath: str)-> List[SeqIO.SeqRecord]:
    records = []
    for rec in SeqIO.parse(fpath, 'fasta'):
        rec.description = ''
        rec.name = rec.id
        records.append(rec)
    return records

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    records = remove_descriptions(args.fasta)
    SeqIO.write(records, args.out, 'fasta')

if __name__ == '__main__':
    main()
