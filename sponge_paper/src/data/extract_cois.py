#!/usr/bin/env python
import os
from Bio import SeqIO

from os.path import expanduser
from glob import glob

COI_IDS = [
    "MK105443.1",
    "MK105444.1",
    "MK105445.1",
    "JX306085.1",
    "JX306086.1",
    "MH749402.1",
    "AJ704977.1",
    "AJ704976.1",
    "YP_001648644.1",
]

home = expanduser("~")
external = os.path.join(home, "sponge_paper", "sponge_paper", "data", "external")
search_path = os.path.join(external, "*.fna")
fpaths = [fp for fp in glob(search_path) if "sponges_COI.fna" not in os.path.basename(fp)]

recs = {}
for fpath in fpaths:
    records = [rec for rec in SeqIO.parse(fpath, "fasta") if rec.id in COI_IDS]
    for record in records:
        if record.id in recs:
            continue
        recs.update({record.id: record})

records = [record for __,record in recs.items()]

outfpath = os.path.join(external, "sponges_COI.fna")
n_seqs = SeqIO.write(records, outfpath, "fasta")
print(f"wrote {n_seqs} sequences to {outfpath}")