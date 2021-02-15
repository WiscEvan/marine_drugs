#Use python antismash_summary.py /path/to/gbk/files
import sys
import os
import pandas as pd
from Bio import SeqIO

def range_subset(range1, range2):
    """Whether range1 is a subset of range2."""
    if not range1:
        return True  # empty range is subset of anything
    if not range2:
        return False  # non-empty range can't be subset of empty range
    if len(range1) > 1 and range1.step % range2.step:
        return False  # must have a single value or integer multiple step
    return range1.start in range2 and range1[-1] in range2
#Use: range_subset(range(0, 1), range(0, 4))

input_directory = sys.argv[1]
df = pd.DataFrame(columns=['Cluster','Spades_Node', 'BGC_start','BGC_end'])
gbk_files = [file for file in os.listdir(input_directory) if file.endswith(".gbk")]
cluster_no = 0
for gbk in gbk_files:
        cluster_no += 1
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
        df = df.append({'Cluster': cluster_no, 'Spades_Node': node, 'BGC_start':start,'BGC_end': end}, ignore_index=True)

df.to_csv(input_directory+'/antismash_summary.txt', sep="\t", index=False)
