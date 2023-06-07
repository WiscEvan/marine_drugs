import pandas as pd
import sys
import os
from Bio import SeqIO
import numpy as np
from functools import reduce

input_directory = sys.argv[1] 
base_name = sys.argv[2] 

#First let's transform each count table into a dataframe
DNA_file = input_directory+"/"+str([file for file in os.listdir(input_directory) if file.endswith("DNA_count.txt")][0])
DNA_count_df = pd.read_csv(DNA_file, sep="\t", names=['Geneid','Chr','Start', 'End', 'Strand','Length','DNA_count'])

print("DNA counts converted to database")

#Next, let's get the corresponding RNA counts. We only need the gene ID and the count.
RNA_file = input_directory+"/"+str([file for file in os.listdir(input_directory) if file.endswith("RNA_count.txt")][0])
tmpRNA_count_df = pd.read_csv(RNA_file, sep="\t", names=['Geneid','Chr','Start', 'End', 'Strand','Length','RNA_count'])
RNA_count_df = tmpRNA_count_df.drop(['Chr','Start', 'End', 'Strand','Length'], axis=1)

print("RNA counts converted to database")

#Now, let's merge that data
DNA_and_RNA_df = pd.merge(DNA_count_df,RNA_count_df,on='Geneid', how='inner')

#Add kingdom and phylum classification data (should all be eukaryotes -mostly porifera!)
taxonomy_df = pd.read_csv(input_directory+'/'+base_name+'.taxonomy.tsv',sep = '\t', header = 0)
taxonomy_df.rename({'contig': 'Chr'}, axis=1, inplace=True)
tax_df = taxonomy_df[['Chr','superkingdom','phylum']]
#Merge data and keep only porifera genes
merge1_df = pd.merge(DNA_and_RNA_df,tax_df,on='Chr', how='inner') 
euk_df = merge1_df.loc[merge1_df['phylum'] == 'porifera']

print("DNA and RNA databases have been successfully merged")

#Next, we'll start calculating the cov weighted and unweighted TPMs
#First, let's get the coverage/abundance of each gene. This is the DNA counts.
euk_df['RPK'] = euk_df['DNA_count']/euk_df['Length']*1000 #Correct for gene length and generate new column
RPK_sum =euk_df["RPK"].sum() #Calculate sum of all length-corrected gene counts
DNA_scaling_factor = RPK_sum/1000000 #Calculate scaling factor
euk_df["RPM"] = euk_df["RPK"]/DNA_scaling_factor #Calculate RPM (Reads per million) for each gene. This corrects for bacterial sequencing depth.

print("DNA counts converted to RPM")

#Now let's do the same thing for the transcripts
euk_df["TPK"] = euk_df["RNA_count"]/euk_df["Length"]*1000 #Correct for gene length and generate new column
TPK_sum =euk_df["TPK"].sum() #Calculate sum of all length-corrected transcript counts
RNA_scaling_factor = TPK_sum/1000000 #Calculate scaling factor
euk_df["TPM"] = euk_df["TPK"]/RNA_scaling_factor #Calculate RPM (Reads per million) for each gene. This corrects for bacterial sequencing depth.

print("RNA counts converted to TPM")

#Create new column with expression (beware /0 errors)
def divide_by_zero_check():
        try:
                euk_df["Expression"] = euk_df["TPM"]/euk_df["RPM"]
        except ZeroDivisionError:
                euk_df["Expression"] = 0

divide_by_zero_check()

#If RPM is zero, make expression value null
euk_df.replace([np.inf, -np.inf], np.nan, inplace=True)

#Let's print this out to a file for a sanity check 
euk_df.to_csv(input_directory+'/'+base_name+'_TPM_interim.tab', sep='\t', index=False)

print("TPM_interim file generated")

#Attach KO metadata
KO_file = input_directory+"/"+str([file for file in os.listdir(input_directory) if file.endswith("_final")][0])
KO_df = pd.read_csv(KO_file, sep='\t',header=0)
master_df = pd.merge(euk_df,KO_df,on='Geneid', how='inner')
master_df.to_csv(input_directory+'/'+base_name+'_eukary_KO_gene_expression.txt', sep='\t', index=False)

#Getting average and median expression per KO pathway
control_KO_df = pd.read_csv(input_directory+'/KEGG_control.txt',sep = '\t', names=['Major Category','Minor category','Pathway','KO','Abbrev','Description'])
annots_df = pd.merge(master_df,control_KO_df,on='KO', how='outer')
pathway_summ_tmp = annots_df[['Pathway','Expression']]
pathway_summ_df = pathway_summ_tmp.fillna(0)
pathway_sum = pathway_summ_df.pivot_table(index='Pathway', values='Expression', aggfunc=np.sum).reset_index()
pathway_mean = pathway_summ_df.pivot_table(index='Pathway', values='Expression', aggfunc=np.mean).reset_index() 
pathway_median = pathway_summ_df.pivot_table(index='Pathway', values='Expression', aggfunc=np.median).reset_index()

#Merge the tables
exp_tables = [pathway_sum,pathway_mean,pathway_median]
df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['Pathway'],how='outer'), exp_tables)

#Rename headers
df_merged  = df_merged.rename(columns={'Expression_x': 'Sum_of_expression'})
df_merged  = df_merged.rename(columns={'Expression_y': 'Average_expression'})
df_merged  = df_merged.rename(columns={'Expression': 'Median_expression'})

#Remove rows that sum to zero
df_final = df_merged.loc[(df_merged.sum(axis=1) != 0)]

df_final.to_csv(input_directory+'/'+base_name+'_eukary_KO_gene_expression_summary.txt', sep='\t', index=False)

#Adding in COG annotations
COG_file = [file for file in os.listdir(input_directory) if file.endswith("_eggNOG_output.emapper.annotations")][0]
COG_df_init = pd.read_csv(input_directory+'/'+COG_file, skiprows=4, sep="\t", header=0)
COG_df_init = COG_df_init.rename(columns={'#query': 'Query'})
COG_df_init.drop(COG_df_init.tail(3).index,inplace = True)
COG_df_init[['Geneid']] = COG_df_init.Query.str.split('.', expand=True)[0]
COG_df = COG_df_init[['Geneid','COG_category']]

#Merge with master and then extract expression per gene COG annot
tmp_COG = pd.merge(euk_df,COG_df, on='Geneid',how = 'inner')
COG_mini1 = tmp_COG[['COG_category','Expression']]
COG_mini1 = COG_mini1.fillna(0)
exploded_cats_df = COG_mini1['COG_category'].apply(list).explode() #Generate new df where category string has been split into new rows
COG_mini=pd.merge(COG_mini1, exploded_cats_df, left_index=True, right_index=True) #Add back in expression data
COG_mini = COG_mini.drop(['COG_category_x'], axis=1)
COG_mini = COG_mini.rename(columns={'COG_category_y': 'COG_category'})
COG_mini.to_csv(input_directory+'/'+base_name+'_eukary_COG_gene_expression.txt', sep='\t', index=False)

#Get sample expression per COG category
COG_expression_matrix_sum = COG_mini.pivot_table(index='COG_category', values='Expression', aggfunc=np.sum).reset_index()
COG_expression_matrix_mean = COG_mini.pivot_table(index='COG_category', values='Expression', aggfunc=np.mean).reset_index()
COG_expression_matrix_median = COG_mini.pivot_table(index='COG_category', values='Expression', aggfunc=np.median).reset_index()

#Merge the tables
exp_tables = [COG_expression_matrix_sum,COG_expression_matrix_mean,COG_expression_matrix_median]
COG_expression_matrix_df = reduce(lambda  left,right: pd.merge(left,right,on=['COG_category'],how='outer'), exp_tables)

#Rename headers
COG_expression_matrix_df  = COG_expression_matrix_df.rename(columns={'Expression_x': 'Sum_of_expression'})
COG_expression_matrix_df  = COG_expression_matrix_df.rename(columns={'Expression_y': 'Average_expression'})
COG_expression_matrix_df  = COG_expression_matrix_df.rename(columns={'Expression': 'Median_expression'})

COG_expression_matrix_df.to_csv(input_directory+'/'+base_name+'_eukary_COG_gene_expression_summary.txt', sep='\t', index=False)

#Counting genes annotated against KEGG and COG databases
all_gene_counts = euk_df['Geneid'].nunique() #Count genes per bin
kegg_gene_counts =master_df['KO'].count() #Count genes with KEGG annotations
COG_gene_counts = tmp_COG['COG_category'].count() #Count genes with COG annotations


#Calculate % annotated genes per bin
Percentage_KO_annotated = kegg_gene_counts/all_gene_counts*100
Percentage_COG_annotated = COG_gene_counts/all_gene_counts*100

#print it out to file
output_file = input_directory+'/'+base_name+'_percent_annots_per_bin.tab'
with open(output_file,'w') as out:
        out.write('Total genes: '+str(all_gene_counts)+'\n'+'Total KEGG annotated genes: '+str(kegg_gene_counts)+'\n'+'Percent KEGG annotated genes: '+str(Percentage_KO_annotated)+'\n'+'Total COG annotated genes: '+str(COG_gene_counts)+'\n'+'Percentage COG annotated genes: '+str(Percentage_COG_annotated)+'\n')
