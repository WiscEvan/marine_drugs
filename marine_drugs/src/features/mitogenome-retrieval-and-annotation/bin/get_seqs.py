#!/usr/bin/env python3

"""Extract annotated sequences from arwen output.txt
command:

- `./arwen  -v -j -br -fo -gcmet -mtx -o ${mtdna/.fna/.arwen.trna} $mtdna`
"""

import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("arwen-output", help="</path/to/arwen/output.txt>")
    parser.add_argument("outfpath", help="</path/to/output/seqs.txt>")
    args = parser.parse_args()
    with open(args.arwen_output) as fh:
        getseq = False
        getfolding = False
        outlines = ""
        for line in fh:
            if line.startswith("Primary sequence for"):
                mtrna = line.strip().split()[-1]
                outlines += ">{}\n".format(mtrna)
                getseq = True
            elif getseq:
                outlines += line.strip() + "\n"
                getseq = False
                getfolding = True
            elif getfolding:
                getfolding = False
                # See: https://stackoverflow.com/questions/5658369/how-to-input-a-regex-in-string-replace
                outlines += (
                    line.replace("t", ".")
                    .replace("A", ".")
                    .replace("d", ".")
                    .replace(" ", ".")
                )
    out = open(args.outfpath, "w")
    out.write(outlines)
    out.close()
    print(f"Written: {args.outfpath}")


if __name__ == "__main__":
    main()
