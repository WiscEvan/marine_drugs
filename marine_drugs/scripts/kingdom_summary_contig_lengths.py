import pandas as pd
import sys
import os
from Bio import SeqIO
import numpy as np

input_directory = '/media/bigdrive1/sam/Sponge_RNA_alignments/Scaffolds/FL2015_*'

master_df = pd.read_csv(input_directory+"/FL2015_*_TPM_master.tab",sep = '\t', header = 0)
contig_length_sum = master_df.pivot_table(index='superkingdom', values='Length', aggfunc=np.sum).reset_index()
contig_length_sum
