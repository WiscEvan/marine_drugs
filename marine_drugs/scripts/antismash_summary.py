#Use python antismash_summary.py /path/to/gbk/files
import sys
import os
import pandas as pd
from Bio import SeqIO

input_directory = sys.argv[1]
df = pd.DataFrame(columns=['BGC','Cluster','Spades_Node', 'BGC_start','BGC_end'])
gbk_files = [file for file in os.listdir(input_directory) if file.endswith(".gbk")]
cluster_no = 0
for gbk in gbk_files:
        cluster_no += 1
        filename = str(gbk).replace(".gbk","")
        with open(input_directory+"/"+gbk, "r") as input:
                for i, line in enumerate(input):
                        if i == 11 and 'NODE' in line:
                                node = str((line.strip('\n').split('::'))[1])
                        elif i ==13 and 'Orig' in line:
                                start = int((line.strip('\n').split('::'))[1])
                        elif i == 14 and 'Orig' in line:
                                end = int((line.strip('\n').split('::'))[1])
                        elif i >14:
                                break
        df = df.append({'BGC':filename, 'Cluster': cluster_no, 'Spades_Node': node, 'BGC_start':start,'BGC_end': end}, ignore_index=True)
df_type = pd.read_csv(input_directory+'/BGC_type.txt', sep='\t', header=0)
df_final = pd.merge(df,df_type,on='BGC', how='inner')
df_final.to_csv(input_directory+'/antismash_summary.txt', sep="\t", index=False)
