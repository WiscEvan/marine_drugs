import pandas as pd
import sys
import os
from Bio import SeqIO

input_directory = sys.argv[1] #path to all files that we need. This includes all count_files, gff files, fna files, fasta files, COG_summary file and antismash summary file
count_files = [file for file in os.listdir(input_directory) if file.endswith("_count.txt")]

#Calculate TPM ber bin
for file_name in count_files:
	file = input_directory+'/'+file_name
	bam_file = file_name.replace("_count.txt",".bam")
	output_file = file.replace("_count.txt","_TPM.txt")
	counts = pd.read_csv(file, sep='\t', skiprows=0, header=1)
	counts["RPK"] = counts[bam_file]/counts["Length"]
	RPK_sum =counts["RPK"].sum()
	scaling_factor = RPK_sum/1000000
	counts["TPM"] = counts["RPK"]/scaling_factor
	final_count = counts.sort_values(by=["TPM"], ascending =False)
	final_count.to_csv(output_file, sep="\t", index=False)

#Next let's concatenate all TPM files together
TPM_files = [file for file in os.listdir(input_directory) if file.endswith("_TPM.txt")]
TPM_df = pd.DataFrame(columns=['Geneid', 'Start', 'End', 'Strand', 'Length', 'Count','RPK','TPM'])
for TPM_file in TPM_files:
	with open (input_directory+'/'+TPM_file, 'r') as input:
		for i in range(1):
			next(input) #Skip the first two lines of the file
		for line in input:
			geneid, chr, start, end, strand, length, count, rpk, tpm = line.strip('\n').split('\t')
			start_int = int(start)
			end_int = int(end)
			length_flt = float(length)
			count_int = int(count)
			#Now let's add each line of the count file to our master dataframe "counts".
			TPM_df = TPM_df.append({'Geneid':geneid, 'Start':start_int, 'End':end_int, 'Strand':strand, 'Length':length_flt, 'Count':count_int, 'RPK':rpk,'TPM':tpm},ignore_index=True)

#Next we're going to add prokka annotations and which contig each gene belongs to
gff_files = [file for file in os.listdir(input_directory) if file.endswith(".gff")]
df1 = pd.DataFrame(columns=['Geneid','Annotation','Chr']) #New dataframe to put all our info in
for file in gff_files:
	with open(input_directory+'/'+file,'r') as input:
		for line in input:
			if not line.startswith('##'):
				if 'Prodigal' in line: #gff files double up on the info - this prevents duplication
					info_list = line.strip('\n').split(";") #Pull geneid and annotation here
					gene_id_index = [idx for idx, s in enumerate(info_list) if 'Parent' in s][0] #Find index of Parent annotation
					geneid_string = str(info_list[gene_id_index]) #This pulls the fule string of "Parent=FL2014_3_0126_00001_gene"
					geneid = geneid_string.split("=")[1] #This just gets at the ID itself
					annot_id_index = [idx for idx, s in enumerate(info_list) if 'product' in s][0] #Find index of Parent annotation
					annot_string = str(info_list[annot_id_index]) #This pulls the full string of "Parent=FL2014_3_0126_00001_gene"
					annot = annot_string.split("=")[1]
					contig = line.strip('\n').split("\t")[0]
					df1 = df1.append({'Geneid': geneid, 'Annotation':annot, 'Chr':contig},ignore_index=True) #Add the data to our dataframe

#Next, let's find our SPAdes node that corresponds to the prokka contig
#Now, we want to bring in the antismash data, but in order to do that we need the correlating Spades NODE ID
dfy = pd.DataFrame(columns=['Chr'])
dfx = pd.DataFrame(columns=['Spades_Node'])
fna_files =[file for file in os.listdir(input_directory) if file.endswith(".fna")]
for fna in fna_files:
	fna_file = input_directory+"/"+fna
	prokka_contigs = [record.id for record in SeqIO.parse(fna_file, "fasta")]
	for contig in prokka_contigs:
		dfy = dfy.append({'Chr': contig}, ignore_index=True)
	filename, file_extension = os.path.splitext(fna)
	fasta = input_directory+'/'+filename+'.fasta'
	spades_contigs = [record.id for record in SeqIO.parse(fasta, "fasta")]
	for node in spades_contigs:
		dfx = dfx.append({'Spades_Node': node}, ignore_index=True)

df2 = pd.concat([dfy,dfx],axis=1,join='inner')

#Now let's add this info to the dataframe with all the annotations, using the Prokka contig number as a key
df3 = pd.merge(df1,df2,on='Chr', how='inner')
#Ab=nd finally, bring in the orginal datafarme of TPM values
df4 = pd.merge(df3,TPM_df,on='Geneid', how='inner')

#Now let's add the COG categories to each gene using our COG_summary.txt file
df5 = pd.read_csv(input_directory+'/COG_summary.txt', sep='\t', header=0)
df6 = pd.merge(df4,df5,on='Geneid', how='outer')
df6.to_csv(input_directory+'/TPM_and_annots.txt', sep='\t', index=False)

#Let's include our data from antismash. We got the summary of our data using antismash_summary.py.
#That script found the NODE and range therein where each BGC is located.
#We can't directly compare gene IDs between prokkka and antiSMASH.
#So we're going to see if each gene on a contig falls into the range of a detected BGC.
#Is it perfect? No. But this approach will be able to tell us if two BGCs on the same contig are expressed at different levels.

#Function: Is range 1 a subset of range 2?
def range_subset(range1, range2):
	if not range1:
		return True  # empty range is subset of anything
	if not range2:
		return False  # non-empty range can't be subset of empty range
	if len(range1) > 1 and range1.step % range2.step:
		return False  # must have a single value or integer multiple step
	return range1.start in range2 and range1[-1] in range2

anti_df = pd.read_csv(input_directory+'/antismash_summary.txt', sep='\t', header=0)
BGC_nodes  = anti_df['Spades_Node'].tolist()
list_of_BGC_nodes = []
for i in BGC_nodes:
    j = i.replace(' ','')
    list_of_BGC_nodes.append(j)
#Let's iterate through each geneide (row) and see if that Spades_contig is in our antismash summary
gene_bgcs = []
for index, row in df6.iterrows():
	spades_node = row['Spades_Node']
	if spades_node in list_of_BGC_nodes:
		gene_start = int(row['Start'])
		gene_end = int(row['End'])
		gene_range = range(gene_start, gene_end+1)
		print('Gene_contig: '+str(spades_node))
		print('Gene_range: '+str(gene_range))
		matching_BGC_contigs = anti_df.loc[anti_df['Spades_Node'] == ' '+spades_node]
		for i, r in matching_BGC_contigs.iterrows():
			print('BGC_node'+str(r['Spades_Node']))
			BGC_start = int(r['BGC_start'])
			BGC_end = int(r['BGC_end'])
			BGC_range = range(BGC_start, BGC_end+1)
			print('BGC_range: '+str(gene_range))
			is_bgc= ''
			print(str(range_subset(gene_range, BGC_range)))
			if range_subset(gene_range, BGC_range) == True:
				is_bgc = 'Yes'
				break
			else:
				is_bgc = 'No'
		gene_bgc = {"Geneid":row['Geneid'], "BGC?": is_bgc}
		gene_bgcs.append(gene_bgc)
	else:
		gene_bgc = {"Geneid":row['Geneid'], "BGC?": 'No'}
		gene_bgcs.append(gene_bgc)

df6.set_index('Geneid', inplace=True)
is_bgc_df = pd.DataFrame(gene_bgcs).set_index('Geneid')
master = pd.merge(df6,is_bgc_df,on='Geneid', how='outer')
master.to_csv(input_directory+'/final_gene_summary.txt', sep='\t', index=True)
