import sys
from Bio import SeqIO
import pandas as pd
from Bio import SeqIO

input_directory = sys.argv[1]
origin_file = sys.argv[2]
tax_file = sys.argv[3]

tax_df = pd.read_csv(input_directory+'/'+tax_file, sep="\t", header = 0)
only_bact_df = tax_df[tax_df['superkingdom'].str.contains("bacteria")] #create subset of contigs that are classified as bacteria

wanted = only_bact_df['contig'].to_list() #Convert contig column to list of contigs

records = (r for r in SeqIO.parse(input_directory+'/'+origin_file, "fasta") if r.id in wanted) #Extract bacterial contigs from scaffold file

output_file = input_directory+'/'+origin_file.split(".")[0] + "_bacterial_only.fasta" #Create output file name
SeqIO.write(records, output_file, "fasta") #Write to output file
