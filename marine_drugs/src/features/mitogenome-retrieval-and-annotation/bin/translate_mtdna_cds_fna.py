#!/usr/bin/env python


import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# from Bio.Data.CodonTable import CodonTable

def translate_dna(fasta: str, out_faa: str = None, translation_table: int = 4) -> str:
    # codon_table =CodonTable.generic_by_id[translation_table]
    aa_records = []
    for dna_record in SeqIO.parse(fasta, 'fasta'):
        # use both fwd and rev sequences
        dna_seqs = [dna_record.seq, dna_record.seq.reverse_complement()]

        # generate all translation frames
        aa_seqs = (s[i:].translate(table=translation_table, to_stop=True) for i in range(3) for s in dna_seqs)

        # select the longest one
        max_aa = max(aa_seqs, key=len)

        # write new record
        record_id = dna_record.id.rsplit("_", 1)[0].replace(";", "")
        description = dna_record.description.replace(record_id, "").replace(";", "").replace("_", "")
        description = f"{record_id} - {description}"
        aa_record = SeqRecord(max_aa, id=record_id, description=description)
        aa_records.append(aa_record)
    
    if out_faa:
        n_written = SeqIO.write(aa_records, out_faa, 'fasta')
        print(f"Wrote {n_written} seqs to {out_faa}")
    return aa_records


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--fasta", help="Path to mtdna.cds.fna file (used for downloaded mtdna from plese et al. 2021", required=True)
    parser.add_argument("--out", help="Path to write output fasta file", required=True)
    parser.add_argument("--code", help="translation table code", choices=range(1, 34), type=int, default=4)
    args = parser.parse_args()
    translate_dna(args.fasta, out_faa=args.out, translation_table=args.code)


if __name__ == '__main__':
    main()