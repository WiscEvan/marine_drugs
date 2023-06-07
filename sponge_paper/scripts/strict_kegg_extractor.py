import sys
from Bio import SeqIO
import os
import pandas as pd

input_directory = sys.argv[1]

file = [file for file in os.listdir(input_directory) if file.endswith("kegg_output_detailed")][0]
kegg_file = input_directory+"/"+file

output_file = kegg_file+"_reliable"
output = open(output_file,'a')
for line in open(kegg_file):
   if line.startswith('*'):
      output.write(line)

output.close()

reliable_kegg_file_txt = [file for file in os.listdir(input_directory) if file.endswith("_reliable")][0] #Read in newly created file with only reliable annotations
#Replace whitespaces with a tab
filepath_in = input_directory+"/"+reliable_kegg_file_txt
filepath_out = input_directory+"/"+reliable_kegg_file_txt+"_tabbed"
with open(filepath_in, 'r') as file_in, open(filepath_out, 'w') as file_out:
    for line in file_in:
        data = line.split()  # splits the content removing spaces and final newline
        line_w_tabs = "\t".join(data) + '\n'
        file_out.write(line_w_tabs)

#Read in and parse the tabbed file
tabbed_kegg_file = [file for file in os.listdir(input_directory) if file.endswith("_tabbed")][0] 
kegg_df = pd.read_csv(input_directory+"/"+tabbed_kegg_file,sep = '\t',usecols=[1,2,4],names=["ChrID","KO","Score"]) #Read in tabbed file into pandas dataframe
kegg_df_dups_removed = kegg_df.groupby(["ChrID","KO"], as_index=False, sort=False)["Score"].max() #Keep annotation with highest score
kegg_df_dups_removed[['Geneid']] = kegg_df_dups_removed.ChrID.str.split('_', expand=True)[1]+"_"+kegg_df_dups_removed.ChrID.str.split('_', expand=True)[6] #Generate Geneid for later merging
kegg_subset_df = kegg_df_dups_removed[["Geneid","KO"]] #Convert to mapper output format
final_file = input_directory+"/"+tabbed_kegg_file+"_final" #Generate empty output file
kegg_subset_df.to_csv(final_file, sep='\t', index=False) #Write subsetted df to file
