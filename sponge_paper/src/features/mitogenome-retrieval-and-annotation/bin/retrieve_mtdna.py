#!/usr/bin/env python
import os
from Bio import SeqIO

from os.path import expanduser

def main():

    home = expanduser("~")
    assemblies = os.path.join(home, "sponge_paper", "sponge_paper", "data", "interim", "assemblies")

    mitogenomes = os.path.join(home, "sponge_paper", "sponge_paper", "data", "interim", "mitogenomes")
    fpath = os.path.join(mitogenomes, "mtdna.txt")
    with open(fpath) as fh:
        mtdna = {}
        for line in fh:
            sample,contig = line.strip().split('\t')
            mtdna.update({sample:contig})

    for sample, contig in mtdna.items():
        assembly = os.path.join(assemblies, f"{sample}.filtered.fna")
        record = [record for record in SeqIO.parse(assembly, "fasta") if record.id == contig]
        out = os.path.join(mitogenomes, f"{sample}.mtdna.fna")
        SeqIO.write(record, out, "fasta")
        print(f"wrote {contig} to {out}")

if __name__ == '__main__':
    main()