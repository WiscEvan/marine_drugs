import pandas as pd
import sys
import os

input_directory = sys.argv[1]

#First, let's make a df of ALL genes
gff_files = [file for file in os.listdir(input_directory) if file.endswith(".gff")]
df1 = pd.DataFrame(columns=['Geneid']) #New dataframe to put all our info in
for file in gff_files:
        with open(input_directory+'/'+file,'r') as input:
                for line in input:
                        if not line.startswith('##'):
                                        if 'Prodigal' in line: #gff files double up on the info - this prevents duplication
                                                info_list = line.strip('\n').split(";") #Pull geneid and annotation here
                                                gene_id_index = [idx for idx, s in enumerate(info_list) if 'Parent' in s][0]
                                                geneid_string = str(info_list[gene_id_index])
                                                geneid= geneid_string.split("=")[1]
                                                df1 = df1.append({'Geneid': geneid},ignore_index=True) #Add the data to our dataframe

#Not all genes have been annotated with a COG cat, so we can't just pull the data to match.
#Instead let's make a second dataframe and we can merge them afterwards
tab_files = [file for file in os.listdir(input_directory) if file.endswith(".tab")]
df2 = pd.DataFrame(columns=['Geneid','COG_ID', 'COG_cat'])
for tab in tab_files:
        with open(input_directory+'/'+tab,'r') as input:
                for line in input:
                        qseqid,stitle,pident,evalue,qlen,slen = line.strip('\n').split('\t')
                        gene_id = str(qseqid)+'_gene'
                        if len(line.strip('\n').split('|')) > 3:
                                print(str(line))
                        else:
                                trash,cog_id,cog_cat = line.strip('\n').split('|')
                                detabbed_cogcat = cog_cat.split('\t')[0]
                                cogid = str(cog_id)
                                cogcat = str(detabbed_cogcat)
                                df2 = df2.append({'Geneid':gene_id,'COG_ID':cogid, 'COG_cat':cogcat}, ignore_index=True)

#Merge dataframes. Geneids without an annotation should come up as NA
COG_df = pd.merge(df1,df2,on='Geneid', how='inner')

#write df to file
COG_df.to_csv(input_directory+'/COG_summary.txt', sep='\t', index=False)
