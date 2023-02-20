#!/usr/bin/env python

"""
Sort tRNA from all *.trna arwen output files.
command:
arwen_exe="${REPO}/marine_drugs/src/features/mitogenome-retrieval-and-annotation/bin/arwen"
- `${arwen_exe} -v -j -br -fo -gcmet -mtx -o ${mtdna/.fna/.arwen.trna} $mtdna`
# arwen 
# -v : -verbose
# -j : Display 4-base sequence on 3' end of astem regardless of predicted
# amino-acyl acceptor length.
# -br : Show secondary structure of tRNA gene primary sequence with round brackets.
# -fo : Print out primary sequence in fasta format only (no secondary structure).
# -gcmet : Use composite Metazoan mitochondrial genetic code.
# -mtx : Low scoring tRNA genes are not reported.
# -o : print output into <outfile>. If <outfile> exists, it is overwritten. 
# By default, output goes to STDOUT.

Filenames will be used to identify the organism's mitogenome.

We will need to appropriately identify all of the headers across files and place
these into their respective tRNA files for uploading to the forna visualizion
server.
"""


import argparse
import os

from glob import glob


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("indir")
    parser.add_argument("outdir")
    parser.add_argument("--glob", default="*.trna")
    args = parser.parse_args()
    search_string = os.path.join(args.indir, args.glob)
    trna_fps = glob(search_string)
    assert len(trna_fps) > 0, "FilesNotFound:\n\t\tsearch string: {}".format(
        search_string
    )
    trnas = {}
    for fp in trna_fps:
        org = os.path.splitext(os.path.basename(fp))[0]
        org_line = ">{}\n".format(org)
        with open(fp) as fh:
            for line in fh:
                if not line.startswith(">"):
                    trnas[trna][org] += line
                    continue
                trna = line.strip(">").strip()
                if trna not in trnas:
                    trnas.update({trna: {org: org_line}})
                else:
                    trnas[trna].update({org: org_line})
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    for trna in trnas:
        # First path handling for trna file
        trna_str = (
            trna.replace("-", "_")
            .replace("|", "_")
            .replace("(", "_")
            .replace(")", "")
            .replace("?", "qm")
        )
        outfpath = os.path.join(args.outdir, f"{trna_str}.seqs")
        # join all lines and write out to trna file
        outlines = "\n".join([trnas[trna][org] for org in trnas[trna]])
        outfile = open(outfpath, "w")
        outfile.write(outlines)
        outfile.close()
        print(f"Written: {os.path.basename(outfpath)}")


if __name__ == "__main__":
    main()
