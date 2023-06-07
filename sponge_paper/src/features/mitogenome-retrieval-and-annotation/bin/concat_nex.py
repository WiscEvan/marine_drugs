#!/usr/bin/env python3


import argparse
import os
from Bio.Nexus import Nexus



def check_taxa(matrices):
    """Verify Nexus instances have the same taxa information.

    Checks that nexus instances in a list [(name, instance)...] have
    the same taxa, provides useful error if not and returns None if
    everything matches
    """
    first_taxa = matrices[0][1].taxlabels
    for name, matrix in matrices[1:]:
        first_only = [t for t in first_taxa if t not in matrix.taxlabels]
        new_only = [t for t in matrix.taxlabels if t not in first_taxa]
        if first_only:
            missing = ', '.join(first_only)
            msg = '%s taxa %s not in martix %s' % (nexi[0][0], missing, name)
            raise Nexus.NexusError(msg)
        elif new_only:
            missing = ', '.join(new_only)
            msg = '%s taxa %s not in all matrices'  % (name, missing)
            raise Nexus.NexusError(msg)
    return None # will only get here if it hasn't thrown an exception


def concat(file_list, same_taxa=True):
    """Combine multiple nexus data matrices in one partitioned file.

    By default this will only work if the same taxa are present in each file
    use same_taxa=False if you are not concerned by this
    """
    nexi = [(fname, Nexus.Nexus(fname)) for fname in file_list]
    if same_taxa:
        if not check_taxa(nexi):
            return Nexus.combine(nexi)
    else:
        return Nexus.combine(nexi)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nex", nargs="+", help="Path to nexus file", required=True)
    parser.add_argument("--out", help="Path to write combined nexus file", required=True)
    args = parser.parse_args()
    # Tuple[str, str]
    # NOTE: 0-th element in tuple is to define 'charset' for the particular nexus file. In this case,
    # we define the charset as simply the gene name
    # the combine function takes a list of tuples [(charset, nexus file)...],
    nexi =  [(os.path.basename(nex).replace(".bmge.nex", ""), Nexus.Nexus(nex)) for nex in args.nex]

    combined = Nexus.combine(nexi)
    combined.write_nexus_data(filename=open(args.out, 'w'))
    print(f"combined {len(nexi)} nexus files to {args.out}")

if __name__ == '__main__':
    main()