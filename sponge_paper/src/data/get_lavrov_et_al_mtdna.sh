#!/usr/bin/env bash

OUTDIR="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/sponge_paper/sponge_paper/data/external/lavrov_et_al"
mkdir -p $OUTDIR

while read accession;do
    # curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${accession}&rettype=gb&retmode=txt">${OUTDIR}/$accession.gbk;
    # See: https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
    if [ ! -f ${OUTDIR}/${accession}.gbk ];then
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=gbwithparts&retmode=txt">${OUTDIR}/$accession.gbk;
        # curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=fasta&retmode=txt">${OUTDIR}/$accession.fna
        echo "downloaded to: ${OUTDIR}/${accession}.gbk"
    fi
done < /media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/sponge_paper/sponge_paper/data/external/lavrov_et_al_gbk_accessions.txt

# wget https://megasun.bch.umontreal.ca/People/lang/FMGP/Cantharellus.html

script="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/sponge_paper/sponge_paper/src/data/extract_seqs_from_gbk.py"
# gbk directory: lavrov_et_al
gbk_dir="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/sponge_paper/sponge_paper/data/external/lavrov_et_al"
gbk_orgs="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/sponge_paper/sponge_paper/data/external/lavrov_et_al/gbk_organism.tsv"
gbks=(`find $gbk_dir -name "*.gbk"`)
python $script \
    --gbk ${gbks[@]} \
    --gbk-orgs $gbk_orgs \
    --outdir $gbk_dir
