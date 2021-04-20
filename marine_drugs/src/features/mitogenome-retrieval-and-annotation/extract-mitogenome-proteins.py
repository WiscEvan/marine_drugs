#!/usr/bin/env python

# Extract mtDNA genes from references and annotated sponge mitogenomes to their
# respective fasta file in sorted order across fastas.
# 1. Retrieve sequences from references (We will need to know what sequences to keep for downstream phylogenetic analysis)
# 2. Retrieve sequences from annotated sponges and keep sorted order across genes
# 3. Write out sample ordered fasta files for each gene
# Following previous mtDNA phylogenetic analysis methods:
# See here: https://mynotebook.labarchives.com/share/Jason%2520Kwan/NDM4LjF8Mzc5NzMxLzMzNy0xMjU1L1RyZWVOb2RlLzE3ODIwMjYyMjF8MTExMi4x
# * Next steps after this extraction:
#   - Alignment of each fasta file to MSA (muscle|mafft)
#   - Trimming of each MSA (BMGE)
#   - Create nexus file of concatenated alignments with information on gene ranges across concatenation
#   - Construct bayesian inference tree with Mr.Bayes using nexus file as input


import os
import glob

from Bio import SeqIO




def main():
    INDIR = "/Users/rees/marine_drugs/marine_drugs/data/interim/mitogenomes"
    # FL2015_37.mitos/eiwIsrSS/result.fas
    search_string = os.path.join(INDIR, "**", "result.faa")
    # Paths to mitogenome annotatoin results
    prots = [
        (
            os.path.basename(os.path.dirname(os.path.dirname(filepath))).replace(
                ".mitos", ""
            ), filepath
        )
        for filepath in glob.glob(search_string, recursive=True)
        if "FL2015_5_code5"
        not in filepath  # Account for instance where mitochondrial genetic code of 5 was used
    ]
    search_string = os.path.join(INDIR, "**", "result.fas")
    nucls = [
        (
            os.path.basename(os.path.dirname(os.path.dirname(filepath))).replace(
                ".mitos", ""
            ), filepath
        )
        for filepath in glob.glob(search_string, recursive=True)
        if "FL2015_5_code5"
        not in filepath  # Account for instance where mitochondrial genetic code of 5 was used
    ]

    samples = {}
    # sponge : protein : {"amino_acid": A.A. seqrecord, "nucleotide": nucl. seqrecord}
    # Now organize all of our sequences together for phylogenetic analysis
    for sponge, filepath in nucls:
        # example: record.id >NODE_4354_length_16774_cov_3815.904588; 1-752; -; rrnL-a
        for record in SeqIO.parse(filepath, 'fasta'):
            protein = record.description.split(";")[-1].strip()
            if protein in samples:
                samples[protein]["nucleotide"].append(record)
            else:
                samples[protein] = {"nucleotide": [record], "amino_acid":[]}

    for sponge, filepath in prots:
        for record in SeqIO.parse(filepath, 'fasta'):
            protein = record.description.split(";")[-1].strip()
            if protein in samples:
                samples[protein]["amino_acid"].append(record)
            else:
                samples[protein] = {"amino_acid": [record], "nucleotide": []}

if __name__ == "__main__":
    main()