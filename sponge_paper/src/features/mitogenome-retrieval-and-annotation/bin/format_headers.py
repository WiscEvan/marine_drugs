#!/usr/bin/env python3
"""
Format headers to so taxa from each partition will match... This will require
the Fehlauer-Ale table with columns, Taxa, COI5P and 16S to allow tranlsation of accessions to Taxa
"""

import argparse
import os

from glob import glob
import pandas as pd

from Bio import SeqIO


def get_translations(fpath):
    df = pd.read_csv(fpath, sep='\t', index_col='Taxa')
    coi = {v:k.rstrip('_') for k,v in df['COI5P'].dropna().to_dict().items()}
    rrnL = {v:k.rstrip('_') for k,v in df['16S'].dropna().to_dict().items()}
    return {'COXI':coi, '16S':rrnL}

def translate_headers(infpath, translation):
    records = []
    for record in SeqIO.parse(infpath, 'fasta'):
        record_baseid = record.id.split('.')[0]
        new_header = translation.get(record_baseid)
        if new_header:
            record.id = new_header
        else:
            print(
                f'\n{8*"-"}Current Record:{8*"-"}\n'
                f'record id:{record.id}\n'
                f'record name:{record.name}\n'
                f'record description:{record.description}\n')
            record.id = input('record id:')
            # print('\n'.join(['-- Record i.d. set --> ',record.id]))
        records.append(record)
    return records

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seqs',
        help='</path/to/[COXI|rrnL]/sequences.fna>', required=True)
    parser.add_argument('--gene-name',
        help='name of gene in sequences.fna',
        choices=['COXI','16S'],required=True)
    parser.add_argument('--transltab',
        help='</path/to/translations/table.tsv>', required=True)
    parser.add_argument('--outdir', help='</path/to/output/directory>',
        default=os.curdir)
    args = parser.parse_args()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    fpaths = [args.seqs]
    translations = get_translations(args.transltab)
    for fpath in fpaths:
        # Translate records
        records = translate_headers(
            fpath,
            translations.get(args.gene_name))
        # Write out translated records
        outfname, ext = os.path.splitext(os.path.basename(fpath))
        outfname = f'{outfname}.reformatted{ext}'
        outfpath = os.path.join(args.outdir, outfname)
        n_written = SeqIO.write(records, outfpath, 'fasta')


if __name__ == '__main__':
    main()
