import pandas as pd
import sys
import os
from Bio import SeqIO
import numpy as np
from functools import reduce

input_directory = sys.argv[1] 
base_name = sys.argv[2] 

#First let's transform each count table into a dataframe
DNA_file = [file for file in os.listdir(input_directory) if file.endswith("DNA_count.txt")][0]
DNA_count_df = pd.read_csv(DNA_file, sep="\t", names=['Geneid','Chr','Start', 'End', 'Strand','Gene_length','DNA_count'])

print("DNA counts converted to database")

#Next, let's get the corresponding RNA counts. We only need the gene ID and the count.
RNA_file = [file for file in os.listdir(input_directory) if file.endswith("RNA_count.txt")][0]
tmpRNA_count_df = pd.read_csv(RNA_file, sep="\t", names=['Geneid','Chr','Start', 'End', 'Strand','Gene_length','RNA_count'])

RNA_count_df = tmpRNA_count_df.drop(['Chr','Start', 'End', 'Strand','Gene_length'], axis=1)

print("RNA counts converted to database")

#Now, let's merge that data
DNA_and_RNA_df = pd.merge(DNA_count_df,RNA_count_df,on='Geneid', how='inner')

print("DNA and RNA databases have been successfully merged")

#Next, we'll start calculating the cov weighted and unweighted TPMs
#First, let's get the coverage/abundance of each gene. This is the DNA counts.
DNA_and_RNA_df['RPK'] = DNA_and_RNA_df['DNA_count']/DNA_and_RNA_df['Gene_length']*1000 #Correct for gene length and generate new column
RPK_sum =DNA_and_RNA_df["RPK"].sum() #Calculate sum of all length-corrected gene counts
DNA_scaling_factor = RPK_sum/1000000 #Calculate scaling factor
DNA_and_RNA_df["RPM"] = DNA_and_RNA_df["RPK"]/DNA_scaling_factor #Calculate RPM (Reads per million) for each gene. This corrects for bacterial sequencing depth.

print("DNA counts converted to RPM")

#Now let's do the same thing for the transcripts
DNA_and_RNA_df["TPK"] = DNA_and_RNA_df["RNA_count"]/DNA_and_RNA_df["Gene_length"]*1000 #Correct for gene length and generate new column
TPK_sum = DNA_and_RNA_df["TPK"].sum() #Calculate sum of all length-corrected transcript counts
RNA_scaling_factor = TPK_sum/1000000 #Calculate scaling factor
DNA_and_RNA_df["TPM"] = DNA_and_RNA_df["TPK"]/RNA_scaling_factor #Calculate RPM (Reads per million) for each gene. This corrects for bacterial sequencing depth.

print("RNA counts converted to TPM")

#Create new column with expression (beware /0 errors)
def divide_by_zero_check():
        try:
                DNA_and_RNA_df["Expression"] = DNA_and_RNA_df["TPM"]/DNA_and_RNA_df["RPM"]
        except ZeroDivisionError:
                DNA_and_RNA_df["Expression"] = 0

divide_by_zero_check()

#If RPM is zero, make expression value null
DNA_and_RNA_df.replace([np.inf, -np.inf], np.nan, inplace=True)

#Calculate relative expression for all genes
total_bacterial_expression = DNA_and_RNA_df["Expression"].sum()
DNA_and_RNA_df["Rel_Expression"] = DNA_and_RNA_df["Expression"]/total_bacterial_expression*100

#Let's print this out to a file for a sanity check 
DNA_and_RNA_df.to_csv(input_directory+'/'+base_name+'_TPM_interim.tab', sep='\t', index=False)

print("TPM_interim file generated")

#Next, let's attach some metadata so we know which genes are part of predicted biosynthetic gene clusters
#Note: We got the summary of our data using antismash_summary.py and then edited this file with a bash oneliner to make it tab separated.
anti_df_init = pd.read_csv(input_directory+'/'+base_name+'_antismash_summary.txt', sep='\t', header=0)
anti_df_init = anti_df_init.rename(columns={'Spades_Node': 'Chr'})
anti_df = anti_df_init.apply(lambda x: x.str.strip() if x.dtype == "object" else x) #Strip whitespace from df

#Function: Is range 1 a subset of range 2?
def range_subset(range1, range2):
        if not range1:
                return True  # empty range is subset of anything
        if not range2:
                return False  # non-empty range can't be subset of empty range
        if len(range1) > 1 and range1.step % range2.step:
                return False  # must have a single value or integer multiple step
        return range1.start in range2 and range1[-1] in range2

#Next, let's identify which genes in the DNA_and_RNA_df are part of BGCs
gene_bgcs = [] #List where we'll store dictionaries of BGC info for each geneid
BGC_nodes = anti_df['Chr'].tolist()
for index, row in DNA_and_RNA_df.iterrows(): #Iterate through rows (Effectively assess every gene)
    BGC_dec ='' #Blank decision
    BGC_type = '' #blank BGC type
    if row['Chr'] not in BGC_nodes: #If contig has no BGCs just list dict entry with No and N/A
        gene_bgc = {"Geneid":row['Geneid'], 'BGC?': 'No', 'Type': 'N/A','Cluster': 'N/A','BGC_name':'N/A', 'Chr':'N/A'} #Create geneid dict
        gene_bgcs.append(gene_bgc) #Add dict to list
    else:
        gene_start = int(row['Start'])
        gene_end = int(row['End'])
        gene_range = range(gene_start, gene_end+1)
        for i,r in anti_df.iterrows():
            if r['Chr'] == row['Chr']: #If the BGC contig is the same as the gene contig continue
                BGC_start = int(r['BGC_start'])
                BGC_end = int(r['BGC_end'])
                BGC_range = range(BGC_start, BGC_end+1)
                if range_subset(gene_range, BGC_range) == True:
                    BGC_dec ='Yes'
                    BGC_cluster = int(r['Cluster'])
                    BGC_contig = r['Chr']
                    BGC_type = r['Type']
                    BGC_name = r['BGC']
                    break
                else:
                    BGC_dec ='No'
                    BGC_cluster = 'N/A'
                    BGC_contig = 'N/A'
                    BGC_type = 'N/A'
                    BGC_name = 'N/A'
        gene_bgc = {"Geneid":row['Geneid'], 'BGC?': BGC_dec, 'Type': BGC_type, 'Cluster':BGC_cluster,'BGC_name':BGC_name, 'Chr':BGC_contig}
        gene_bgcs.append(gene_bgc)

print("Antismash metadata captured")

#Convert list of dicts to dataframe
is_bgc_df = pd.DataFrame(gene_bgcs)

#Merge dfs, remove extra Chr column and rename remaining Chr column
tmp_master_df1 = pd.merge(DNA_and_RNA_df,is_bgc_df,on=['Geneid'], how='outer')
tmp_master_df = tmp_master_df1.drop(['Chr_y'], axis=1)
tmp_master_df = tmp_master_df.rename(columns={'Chr_x': 'Chr'})

#Add kingdom classification data (should all be bacteria!)
taxonomy_df = pd.read_csv(input_directory+'/'+base_name+'.taxonomy.tsv',sep = '\t', header = 0)
taxonomy_df.rename({'contig': 'Chr'}, axis=1, inplace=True)
tax_df = taxonomy_df[['Chr','superkingdom']]
pre_master_df = pd.merge(tmp_master_df,tax_df,on='Chr', how='inner') 

#Generate new df with info of which contigs belong to which bin
bins_df = pd.read_csv(input_directory+'/contigs_per_bin.tab', sep='\t', names = ['Bin','Chr'])
sample_bin_df = bins_df[bins_df['Bin'].str.contains(base_name)] 

#Merge binning data and expression data
master_df = pd.merge(pre_master_df,sample_bin_df,on='Chr', how='outer') 

#Replace N/A with "Unbinned" for all contigs/genes that were not binned
master_df.Bin = master_df.Bin.fillna('Unbinned')

#Get contig lengths and remove contigs where contig length is NaN
master_df[['Contig_Length']] = master_df.Chr.str.split('_', expand=True)[3]

master_df.to_csv(input_directory+'/'+base_name+'_TPM_master.tab', sep='\t', index=False)

print("Master file successfully generated!")

#Getting expression stats

#Summed expression per bin
expression_summ_df = master_df[['Bin','Rel_Expression']]
expression_matrix_sum = expression_summ_df.pivot_table(index='Bin', values='Rel_Expression', aggfunc=np.sum).reset_index()
expression_matrix_mean = expression_summ_df.pivot_table(index='Bin', values='Rel_Expression', aggfunc=np.mean).reset_index()
expression_matrix_median = expression_summ_df.pivot_table(index='Bin', values='Rel_Expression', aggfunc=np.median).reset_index()

#Merge the tables
exp_tables = [expression_matrix_sum,expression_matrix_mean,expression_matrix_median]
expression_matrix_df = reduce(lambda  left,right: pd.merge(left,right,on=['Bin'],how='outer'), exp_tables)

#Rename headers
expression_matrix_df  = expression_matrix_df.rename(columns={'Rel_Expression_x': 'Sum_of_expression'})
expression_matrix_df  = expression_matrix_df.rename(columns={'Rel_Expression_y': 'Average_expression'})
expression_matrix_df  = expression_matrix_df.rename(columns={'Rel_Expression': 'Median_expression'})

#Write it out to file
expression_matrix_df.to_csv(input_directory+'/'+base_name+'_Bin_expression_summary.tab', sep='\t', index=False)

#Get and write out a subset of the master file with only BGC genes
BGC_subset = master_df.loc[master_df['BGC?'] == 'Yes']
BGC_subset.to_csv(input_directory+'/'+base_name+'_detailed_BGC_expression.tab', sep='\t', index=False)

#Getting BGC expression stats
bgc_expression_summ_df = master_df[['BGC_name','Rel_Expression']]

#Remove rows that are not BGCs
bgc_expression_summ_df = bgc_expression_summ_df[bgc_expression_summ_df.BGC_name != "N/A"]

bgc_expression_matrix_sum = bgc_expression_summ_df.pivot_table(index='BGC_name', values='Rel_Expression', aggfunc=np.sum).reset_index()
bgc_expression_matrix_mean = bgc_expression_summ_df.pivot_table(index='BGC_name', values='Rel_Expression', aggfunc=np.mean).reset_index()
bgc_expression_matrix_median = bgc_expression_summ_df.pivot_table(index='BGC_name', values='Rel_Expression', aggfunc=np.median).reset_index()

#Merge the tables
exp_tables = [bgc_expression_matrix_sum,bgc_expression_matrix_mean,bgc_expression_matrix_median]
bgc_expression_matrix_df = reduce(lambda  left,right: pd.merge(left,right,on=['BGC_name'],how='outer'), exp_tables)

#Rename headers
bgc_expression_matrix_df  = bgc_expression_matrix_df.rename(columns={'Rel_Expression_x': 'Sum_of_expression'})
bgc_expression_matrix_df  = bgc_expression_matrix_df.rename(columns={'Rel_Expression_y': 'Average_expression'})
bgc_expression_matrix_df  = bgc_expression_matrix_df.rename(columns={'Rel_Expression': 'Median_expression'})

#Add binning info
BGC_bins_df = master_df[['BGC_name','Bin','Type','Cluster']]
bgc_expression_matrix_df_bins = pd.merge(bgc_expression_matrix_df,BGC_bins_df,on='BGC_name',how='inner')
bgc_expression_matrix_df_bins.drop_duplicates(keep='first',inplace=True)

#Write it out to file
bgc_expression_matrix_df_bins.to_csv(input_directory+'/'+base_name+'_BGC_expression_summary.tab', sep='\t', index=False)

#Adding in KEGG annotations
KEGG_file = [file for file in os.listdir(input_directory) if file.endswith("detailed_reliable_tabbed_final")][0]
KEGG_df = pd.read_csv(input_directory+'/'+KEGG_file, sep="\t", header=0)
tmp_KO = pd.merge(master_df,KEGG_df, on='Geneid',how = 'inner')

#Create matrix of expression per KO and output to file
KO_bact_minimized = tmp_KO[['KO','Rel_Expression']]
sample_KO_matrix_df = KO_bact_minimized.pivot_table(index='KO', values='Rel_Expression', aggfunc=np.sum).reset_index() 

#Read in control file
KEGG_control = pd.read_csv(input_directory+'/KEGG_control.txt', sep="\t", header=0)
#Merge to get expression
KEGG_control_df = pd.merge(KEGG_control,sample_KO_matrix_df,on='KO',how='outer')

#Get expression stats per pathway
pathway_expression_summ_df = KEGG_control_df[['Pathway','Rel_Expression']]
pathway_expression_matrix_sum = pathway_expression_summ_df.pivot_table(index='Pathway', values='Rel_Expression', aggfunc=np.sum).reset_index()
pathway_expression_matrix_mean = pathway_expression_summ_df.pivot_table(index='Pathway', values='Rel_Expression', aggfunc=np.mean).reset_index()
pathway_expression_matrix_median = pathway_expression_summ_df.pivot_table(index='Pathway', values='Rel_Expression', aggfunc=np.median).reset_index()

#Merge the tables
exp_tables = [pathway_expression_matrix_sum,pathway_expression_matrix_mean,pathway_expression_matrix_median]
pathway_expression_matrix_df = reduce(lambda  left,right: pd.merge(left,right,on=['Pathway'],how='outer'), exp_tables)

#Rename headers
pathway_expression_matrix_df  = pathway_expression_matrix_df.rename(columns={'Rel_Expression_x': 'Sum_of_expression'})
pathway_expression_matrix_df  = pathway_expression_matrix_df.rename(columns={'Rel_Expression_y': 'Average_expression'})
pathway_expression_matrix_df  = pathway_expression_matrix_df.rename(columns={'Rel_Expression': 'Median_expression'})

#Write it out to file
pathway_expression_matrix_df.to_csv(input_directory+'/'+base_name+'_sample_KO_pathway_expression_summary.tab', sep='\t', index=False)

#Redo all of the above but get stats per bin
KO_bin_minimized = tmp_KO[['KO','Bin','Rel_Expression']]
KO_bin_expression = KO_bin_minimized.pivot_table(columns='Bin', index='KO', values='Rel_Expression',aggfunc=np.sum).reset_index()
#Convert NaN values to 0 expression
KO_bin_expression = KO_bin_expression.fillna(0)
#print out detailed matrix
KO_bin_expression.to_csv(input_directory+'/'+base_name+'_detailed_bin_KO_expression.tab', sep='\t', index=False)

#Tack on pathway info
KEGG_bin_control_df = KEGG_control[['Pathway','KO']]
KEGG_bin_pathway_merged = pd.merge(KEGG_bin_control_df,KO_bin_expression,on='KO',how='outer')

#Get average expression per pathway per bin
bin_pathway_mean_df = KEGG_bin_pathway_merged.groupby('Pathway', as_index=False).mean()
bin_pathway_mean_df.to_csv(input_directory+'/'+base_name+'_mean_pathway_expression_per_bin.tab', sep='\t', index=False)

#Get median expression per pathway per bin
bin_pathway_median_df = KEGG_bin_pathway_merged.groupby('Pathway', as_index=False).median()
bin_pathway_median_df.to_csv(input_directory+'/'+base_name+'_median_pathway_expression_per_bin.tab', sep='\t', index=False)

#Adding in COG annotations
COG_file = [file for file in os.listdir(input_directory) if file.endswith("_eggNOG_output.emapper.annotations")][0]
COG_df_init = pd.read_csv(input_directory+'/'+COG_file, skiprows=4, sep="\t", header=0)
COG_df_init = COG_df_init.rename(columns={'#query': 'Query'})
COG_df_init.drop(COG_df_init.tail(3).index,inplace = True)
COG_df_init[['Geneid']] = COG_df_init.Query.str.split('_', expand=True)[1]+"_"+COG_df_init.Query.str.split('_', expand=True)[6]
COG_df = COG_df_init[['Geneid','COG_category']]
COG_df  = COG_df .replace('-', np.nan) #Replace not annotated with Nan

#Merge with master and then extract expression per gene COG annot
tmp_COG = pd.merge(master_df,COG_df, on='Geneid',how = 'inner')
COG_mini1 = tmp_COG[['COG_category','Rel_Expression']]
COG_mini1 = COG_mini1.fillna(0)
exploded_cats_df = COG_mini1['COG_category'].apply(list).explode() #Generate new df where category string has been split into new rows
COG_mini=pd.merge(COG_mini1, exploded_cats_df, left_index=True, right_index=True) #Add back in expression data
COG_mini = COG_mini.drop(['COG_category_x'], axis=1)
COG_mini = COG_mini.rename(columns={'COG_category_y': 'COG_category'})

#Get sample expression per COG category
COG_expression_matrix_sum = COG_mini.pivot_table(index='COG_category', values='Rel_Expression', aggfunc=np.sum).reset_index()
COG_expression_matrix_mean = COG_mini.pivot_table(index='COG_category', values='Rel_Expression', aggfunc=np.mean).reset_index()
COG_expression_matrix_median = COG_mini.pivot_table(index='COG_category', values='Rel_Expression', aggfunc=np.median).reset_index()

#Merge the tables
exp_tables = [COG_expression_matrix_sum,COG_expression_matrix_mean,COG_expression_matrix_median]
COG_expression_matrix_df = reduce(lambda  left,right: pd.merge(left,right,on=['COG_category'],how='outer'), exp_tables)

#Rename headers
COG_expression_matrix_df  = COG_expression_matrix_df.rename(columns={'Rel_Expression_x': 'Sum_of_expression'})
COG_expression_matrix_df  = COG_expression_matrix_df.rename(columns={'Rel_Expression_y': 'Average_expression'})
COG_expression_matrix_df  = COG_expression_matrix_df.rename(columns={'Rel_Expression': 'Median_expression'})

#Write it out to file
COG_expression_matrix_df.to_csv(input_directory+'/'+base_name+'_sample_COG_expression_summary.tab', sep='\t', index=False)

#Redo all of the above but get stats per bin
COG_bins1 = tmp_COG[['COG_category','Bin','Rel_Expression']]
COG_bins1 = COG_bins1.fillna(0)
exploded_cats_bins_df = COG_bins1['COG_category'].apply(list).explode() #Generate new df where category string has been split into new rows
COG_bins = pd.merge(COG_bins1, exploded_cats_bins_df, left_index=True, right_index=True) #Add back in expression data
COG_bins = COG_bins.drop(['COG_category_x'], axis=1)
COG_bins = COG_bins.rename(columns={'COG_category_y': 'COG_category'})

#Get the summaries and print it out
COG_bin_expression_mean = COG_bins.pivot_table(columns='Bin', index='COG_category', values='Rel_Expression',aggfunc=np.mean).reset_index()
COG_bin_expression_mean.to_csv(input_directory+'/'+base_name+'_binned_COG_expression_mean.tab', sep='\t', index=False)

COG_bin_expression_median = COG_bins.pivot_table(columns='Bin', index='COG_category', values='Rel_Expression',aggfunc=np.median).reset_index()
COG_bin_expression_median.to_csv(input_directory+'/'+base_name+'_binned_COG_expression_median.tab', sep='\t', index=False)

#Counting genes annotated against KEGG and COG databases
all_gene_counts = master_df.groupby('Bin')['Geneid'].nunique() #Count genes per bin
kegg_gene_counts =tmp_KO.groupby('Bin')['KO'].count() #Count genes with KEGG annotations
COG_gene_counts = tmp_COG.groupby('Bin')['COG_category'].count() #Count genes with COG annotations

#merge the dfs
perc_annots_pre = pd.merge(all_gene_counts,kegg_gene_counts, on='Bin',how = 'outer')
perc_annots =pd.merge(perc_annots_pre,COG_gene_counts,on='Bin',how = 'outer')

#Calculate % annotated genes per bin
perc_annots['Percentage_KO_annotated'] = perc_annots['KO']/perc_annots['Geneid']*100
perc_annots['Percentage_COG_annotated'] = perc_annots['COG_category']/perc_annots['Geneid']*100

#print it out to file
perc_annots.to_csv(input_directory+'/'+base_name+'_percent_annots_per_bin.tab', sep='\t', index=True)
