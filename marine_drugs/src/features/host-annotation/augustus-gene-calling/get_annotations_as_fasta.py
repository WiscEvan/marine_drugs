#!/usr/bin/env python
"""
Example annotation from Augustus GFF file:

# start gene g16742
NODE_61320_length_901_cov_22.192020	AUGUSTUS	gene	2	901	0.11	+	.	g16742
NODE_61320_length_901_cov_22.192020	AUGUSTUS	transcript	2	901	0.11	+	.	g16742.t1
NODE_61320_length_901_cov_22.192020	AUGUSTUS	intron	761	819	0.48	+	.	transcript_id "g16742.t1"; gene_id "g16742";
NODE_61320_length_901_cov_22.192020	AUGUSTUS	CDS	2	760	0.11	+	2	transcript_id "g16742.t1"; gene_id "g16742";
NODE_61320_length_901_cov_22.192020	AUGUSTUS	exon	2	760	.	+	.	transcript_id "g16742.t1"; gene_id "g16742";
NODE_61320_length_901_cov_22.192020	AUGUSTUS	CDS	820	884	0.49	+	2	transcript_id "g16742.t1"; gene_id "g16742";
NODE_61320_length_901_cov_22.192020	AUGUSTUS	exon	820	901	.	+	.	transcript_id "g16742.t1"; gene_id "g16742";
NODE_61320_length_901_cov_22.192020	AUGUSTUS	stop_codon	882	884	.	+	0	transcript_id "g16742.t1"; gene_id "g16742";
NODE_61320_length_901_cov_22.192020	AUGUSTUS	tts	901	901	.	+	.	transcript_id "g16742.t1"; gene_id "g16742";
# coding sequence = [ccgcattcaaatccttgagttcccccggcagtactccagggggtgtgataggaggggtggaaaaggacacagacgagct
# gttttgttcgcacactgtgaaagagattaaacagattgagcagcagacgcggctggacattgataggagaaaggaagacctgagacaaatggttggtg
# agaggtaccgtgatctaattgatgctgcagactcaattgtggacatgaaaaaatgctcagagatgatgctgagcagtgtgtccagtatggaacaggtg
# tgccagcgtctccaacagacctacttcttgaaaggctctggagggatgaggggacctcctgaaaaggacaagagccatgaagtggcagcccagttgaa
# gctactagccgatgcccctgagaagatttggcagtgtctcgaagatcagagctacctccaggctgcacagtattttctgctctccagacacatcttct
# ctcacctggagttgagtggagcaggctcagatgttgtaaagaaactctggcactccgtctctagcttcaaggataccataatggagtgttgcagctgc
# catcttcagtgctcagaccttgaggaaggtccaatggtccagtccctgtgtgccatgatgctcctagagggttgctcacctaggcaagcctttgccaa
# gtttcttgtggctaggaaggcagccattcagtctgtgttccaggcagccattttgggggaggggccccaccgcagtgtcaaggcacaagtgtcaaagg
# atgccagatgggagcttaggtccttctctactgctgcaaaagctgacacagttcactag]
# protein sequence = [AFKSLSSPGSTPGGVIGGVEKDTDELFCSHTVKEIKQIEQQTRLDIDRRKEDLRQMVGERYRDLIDAADSIVDMKKCS
# EMMLSSVSSMEQVCQRLQQTYFLKGSGGMRGPPEKDKSHEVAAQLKLLADAPEKIWQCLEDQSYLQAAQYFLLSRHIFSHLELSGAGSDVVKKLWHSV
# SSFKDTIMECCSCHLQCSDLEEGPMVQSLCAMMLLEGCSPRQAFAKFLVARKAAIQSVFQAAILGEGPHRSVKAQVSKDARWELRSFSTAAKADTVH]
# Evidence for and against this transcript:
# % of transcript supported by hints (any source): 0
# CDS exons: 0/2
# CDS introns: 0/1
# 5'UTR exons and introns: 0/0
# 3'UTR exons and introns: 0/1
# hint groups fully obeyed: 0
# incompatible hint groups: 0
# end gene g16742
"""


import os

from Bio import SeqIO

def parse_original_fasta(fasta):
    return [record for record in SeqIO.parse(fasta, "fasta")]


def retrieve_annotations(gff):
    parsing = False
    with open(gff) as fh:
        # One line will have 'coding sequence = [...'
        # another line will have 'protein sequence = [...', 
        # For both should check for ']' in same and subsequent lines
        # Record recorded prior to 

        for line in fh:
            line = line.strip()
            assert llist == 8
            if parsing:
                llist = line.split("\t")
                seqname = llist[0]
                feature = llist[1]
                start = llist[2]
                end = llist[3]
                strand = llist[4]
                frame = llist[5]
                transcript_id = llist[6]
                seqname = llist[7]
            if "start gene" in line:
                parsing = not parsing
            if "end gene" in line:
                parsing = not parsing
                break
            




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta")


    records = parse_original_fasta


if __name__ == "__main__":
    main()