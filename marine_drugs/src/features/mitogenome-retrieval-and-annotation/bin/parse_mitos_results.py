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


import argparse
import os
import glob
from typing import Dict

from Bio import SeqIO


def get_nucleotide_seqs(mitos: str, samples: Dict[str,Dict[str,Dict[str,SeqIO.SeqRecord]]] = None) -> Dict[str,Dict[str,Dict[str,SeqIO.SeqRecord]]]:
    # FL2015_37.mitos/eiwIsrSS/result.fas
    nucls_search_string = os.path.join(mitos, "**", "result.fas")
    nucls = [fp for fp in glob.glob(nucls_search_string, recursive=True)]
    samples = {}
    # sponge : protein : {"amino_acid": A.A. seqrecord, "nucleotide": nucl. seqrecord}
    for fp in nucls:
        job_identifier = os.path.basename(os.path.dirname(os.path.dirname(fp)))
        sponge_id = job_identifier.replace("_code4_refseq89_metazoa", "")
        for record in SeqIO.parse(fp, "fasta"):
            protein_name = record.description.split(';')[-1].strip()
            if sponge_id not in samples:
                samples[sponge_id] = {protein_name: {"nucleotide": [record], "amino_acid":[]}}
            elif protein_name not in samples[sponge_id]:
                samples[sponge_id][protein_name] = {"nucleotide": [record], "amino_acid":[]}
            else:
                samples[sponge_id][protein_name]["nucleotide"].append(record)
    return samples


def get_protein_seqs(mitos: str, samples: Dict[str,Dict[str,Dict[str,SeqIO.SeqRecord]]] = None) -> Dict[str,Dict[str,Dict[str,SeqIO.SeqRecord]]]:
    prots_search_string = os.path.join(mitos, "**", "result.faa")
    prots = [fp for fp in glob.glob(prots_search_string, recursive=True)]
    for fp in prots:
        job_identifier = os.path.basename(os.path.dirname(os.path.dirname(fp)))
        sponge_id = job_identifier.replace("_code4_refseq89_metazoa", "")
        for record in SeqIO.parse(fp, "fasta"):
            protein_name = record.description.split(';')[-1].strip()
            # record.id is contig_id
            # record.description: NODE_703_length_19394_cov_1408.958694; 4902-5081; +; atp8
            if sponge_id not in samples:
                samples[sponge_id] = {protein_name: {"nucleotide": [], "amino_acid":[record]}}
            elif protein_name not in samples[sponge_id]:
                samples[sponge_id][protein_name] = {"nucleotide": [], "amino_acid":[record]}
            else:
                samples[sponge_id][protein_name]["amino_acid"].append(record)
    return samples


def sort_samples(samples):
    sorted_proteins = {}
    for sponge,proteins in samples.items():
        print(f"{sponge} contains {len(proteins)} annotations")
        for protein_name,records in proteins.items():
            nucl_record = records['nucleotide']
            aa_record = records['amino_acid']
            if protein_name in sorted_proteins:
                sorted_proteins[protein_name]["nucleotide"].extend(nucl_record)
                sorted_proteins[protein_name]["amino_acid"].extend(aa_record)
            else:
                sorted_proteins[protein_name] = {
                    "amino_acid":aa_record,
                    "nucleotide":nucl_record,
                }
    return sorted_proteins


def write_samples(sorted_proteins, outdir="sorted_mtdna"):
    os.makedirs(outdir, exist_ok=True)
    for protein_name, records in sorted_proteins.items():
        nucl_out = os.path.join(outdir, f"{protein_name}.fna")
        aa_out = os.path.join(outdir, f"{protein_name}.faa")
        n_aa_written = SeqIO.write(records['amino_acid'], aa_out,'fasta')
        n_nucl_written = SeqIO.write(records['nucleotide'], nucl_out,'fasta')
        print(f"Wrote {n_aa_written} a.a. seqs")
        print(f"Wrote {n_nucl_written} nucl. seqs")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mitos',
        help="base directory path holding mitos results (will recursively search for result.fa{a,s} in this directory)",
        default="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/data/interim/mitogenomes/mitos_results")
    parser.add_argument('--outdir',
        help="Path to output directory",
        default="mtdna_annotations")
    args = parser.parse_args()
    samples = get_nucleotide_seqs(args.mitos)
    samples = get_protein_seqs(args.mitos, samples)
    # Now organize all of our sequences together for phylogenetic analysis
    sorted_samples = sort_samples(samples)
    # protein -> (sorted in same order as sponge)
    # Get minimal set
    write_samples(sorted_samples, outdir=args.outdir)

    
    # [fp for fp in glob.glob(args.outdir, "*.faa") if len([rec for rec in SeqIO.parse(fp, 'fasta')]) == 12]

if __name__ == "__main__":
    main()