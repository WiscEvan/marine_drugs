#/usr/bin/env python

import pandas as pd

# GAF 2.1 format: http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/
# Column	Content	Required?	Cardinality	Example	 
# 1	DB	required	1	UniProtKB	 
# 2	DB Object ID	required	1	P12345	 
# 3	DB Object Symbol	required	1	PHO3	 
# 4	Qualifier	optional	0 or greater	NOT	 
# 5	GO ID	required	1	GO:0003993	 
# 6	DB:Reference (|DB:Reference)	required	1 or greater	PMID:2676709	 
# 7	Evidence Code	required	1	IMP	 
# 8	With (or) From	optional	0 or greater	GO:0000346	 
# 9	Aspect	required	1	F	 
# NOTE: Column 9 Refers to the GO ID (column 5) belongs;
# - one of P (biological process)
# - F (molecular function)
# - C (cellular component).
# 10	DB Object Name	optional	0 or 1	Toll-like receptor 4	 
# 11	DB Object Synonym (|Synonym)	optional	0 or greater	hToll	Tollbooth
# 12	DB Object Type	required	1	protein	 
# 13	Taxon(|taxon)	required	1 or 2	taxon:9606	 
# 14	Date	required	1	20090118	 
# 15	Assigned By	required	1	SGD	 
# 16	Annotation Extension	optional	0 or greater	part_of(CL:0000576)	 
# 17	Gene Product Form ID	optional	0 or 1	UniProtKB:P12345-2	 




