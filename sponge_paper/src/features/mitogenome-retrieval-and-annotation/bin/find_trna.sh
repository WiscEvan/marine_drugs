#!/usr/bin/env bash

REPO = "/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/sponge_paper"
export PATH="$REPO/sponge_paper/src/features/mitogenome-retrieval-and-annotation/bin":$PATH

mitogenomes="$REPO/sponge_paper/sponge_paper/data/interim/mitogenomes"
# Search each mitogenome for tRNA
for mtdna in `ls ${mitogenomes}/*.mtdna.fna`;do
    # -v: Verbose. Prints out search progress
    # -j: Display 4-base sequence on 3' end of astem regardless of predicted amino-acyl acceptor length.
    # -gcinvert: Use Invertebrate mitochondrial genetic code
    # -gcmet: Use composite Metazoan mitochondrial genetic code.
    # mtx: Low scoring tRNA genes are not reported.
    # o: output
    # -gc<num>: Use the GenBank transl_table = <num> genetic code.
    # -gcprot: Use Mold/Protozoan/Coelenterate mitochondrial genetic code. (equivalent to transl_table=4)
    arwen -v -j -br -gcprot -mtx -o ${mtdna/.fna/.gcprot.arwen.txt} $mtdna # (transl_table=4) -> this was used for hippospongia lachne mtDNA
    arwen -v -j -br -gcinvert -mtx -o ${mtdna/.fna/.gcinvert.arwen.txt} $mtdna # (transl_table=5)
    arwen -v -j -br -gcmet -mtx -o ${mtdna/.fna/.gcmet.arwen.txt} $mtdna
    # get_seqs.py ${mtdna/.fna/.arwen.out} ${mtdna/.fna/.trna}
    genbank_format.py --out ${mtdna/.fna/.gcprot.tRNA.gbk} ${mtdna/.fna/.gcprot.arwen.txt} ${mtdna}
    genbank_format.py --out ${mtdna/.fna/.gcinvert.tRNA.gbk} ${mtdna/.fna/.gcinvert.arwen.txt} ${mtdna}
    genbank_format.py --out ${mtdna/.fna/.gcmet.tRNA.gbk} ${mtdna/.fna/.gcmet.arwen.txt} ${mtdna}
done

# python sort_trna.py $curdir "${curdir}/tRNA_with_secondary_structure" --glob "*.seqs"
