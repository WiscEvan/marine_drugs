import sys
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC

input_directory = sys.argv[1]

master_df = pd.DataFrame(columns=['Bin','Size','Cov','GC','#Contigs','N50','Longest_contig'])

bins = [file for file in os.listdir(input_directory) if file.endswith(".fasta")]

for file in bins:
        bin = input_directory+'/'+file
        bin_name = str(file).replace('.fasta','')

        #Reset all variables
        num_of_contigs = 0
        records =[]
        seq_length_list = []
        sequence_total = 0

        #Count no. of contigs
        records = list(SeqIO.parse(bin, "fasta"))
        num_of_contigs = len(records)
        #Add up length of all contigs and then format to bi in Mbp
        seq_length_list = [len(rec) for rec in SeqIO.parse(bin, "fasta")]
        assembly_size = sum(seq_length_list)
        half_size  = assembly_size/2
        size = assembly_size/1000000
        longest_contig = max(seq_length_list)
        #Find N50. This was tricky.
        #You gotta sort your contig lengths and then working from biggest to smallest, start adding them together.
        #Once that summed total reaches half the size of the genome, we've found our bin N50.
        sorted_list = sorted(seq_length_list, reverse = True)
        for num in sorted_list:
                sequence_total += int(num)
                if sequence_total > half_size:
                        n50 = num
                        break

        #Reset GC and coverages
        GC_total = 0
        weighted_cov = 0
        #Calculate GC and coverage corrected for contig length
        for record in SeqIO.parse(bin, "fasta"):
                header_string = str(record.id)
                node,node_no,leng,length,co,cov = header_string.replace('>','').split('_')
                coverage = float(cov)
                length = len(record)
                frac_length = length/assembly_size
                weighted_cov += coverage * frac_length
                GC_total += GC(record.seq)*frac_length
        #Add info to master dataframe
        master_df = master_df.append({'Bin':bin_name,'Size':size,'Cov':weighted_cov,'GC':GC_total,'#Contigs':num_of_contigs,'N50':n50,'Longest_contig':longest_contig},ignore_index=True)

#Print it all out to a file
master_df.to_csv(input_directory+'/bin_summary.txt', sep='\t', index=False)
