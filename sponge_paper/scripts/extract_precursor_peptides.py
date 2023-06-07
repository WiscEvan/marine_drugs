#Usage python extract_precursor_peptides.py /path/to/antismash/gbk_files

from Bio import SeqIO
import pandas as pd
import sys
import os
import numpy as np

input_directory = sys.argv[1]

gbk_files = [file for file in os.listdir(input_directory) if file.endswith(".gbk")]

output_string =''

for gbk_file in gbk_files:
	with open(input_directory+"/"+gbk_file, "r") as input:
		for i, line in enumerate(input):
			if i == 11 and 'NODE' in line:
				node = str((line.strip('\n').split('::'))[1])
			elif i ==13 and 'Orig' in line:
				bgc_start = int((line.strip('\n').split('::'))[1])
			elif i == 14 and 'Orig' in line:
				bgc_end = int((line.strip('\n').split('::'))[1])
			elif i >14:
				break
	count = 0
	BGC = SeqIO.read(input_directory+"/"+gbk_file,'genbank')
	file_name = str(gbk_file)
	molecule = file_name.replace('.gbk','')
	for feat in BGC.features:
		if feat.type == 'CDS':
			if 'gene_functions' not in feat.qualifiers:
				None
			elif "TIGR03793" in str(feat.qualifiers['gene_functions'][0]):
				descrip = feat.qualifiers['locus_tag'][0]
				gene_prot = feat.qualifiers['translation'][0]
				new_start = bgc_start + feat.location.start+1
				new_end = bgc_start + feat.location.end
				output_string = output_string+">"+molecule+"_"+str(new_start)+"_"+str(new_end)+"\n"+gene_prot+"\n"

if len(output_string) > 0:
	outputfile = "Extracted_precursor_peptides.faa"
	with open(input_directory+"/"+outputfile, "w") as text_file:
		text_file.write(output_string)
else:
	pass

