#!/usr/bin/env bash

OUTDIR="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/data/external/plese_et_al_2021"
mkdir -p $OUTDIR

while read accession;do
    # curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${accession}&rettype=gb&retmode=txt">${OUTDIR}/$accession.gbk;
    # See: https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
    if [ ! -f ${OUTDIR}/${accession}.gbk ];then
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=gbwithparts&report=gbwithparts&retmode=txt">${OUTDIR}/$accession.gbk;
        echo "downloaded to: ${OUTDIR}/${accession}.gbk"
    fi
    # curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=fasta&retmode=txt">${OUTDIR}/$accession.fna
done < /media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/data/external/plese_et_al_gbk_accessions.txt


wget -O ${OUTDIR}/stryphnus_fortis.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc102.txt
wget -O ${OUTDIR}/poecillastra_compressa.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc100.txt
wget -O ${OUTDIR}/cymbaxinella_damicornis.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc103.txt
wget -O ${OUTDIR}/phakellia_ventilabrum.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc106.txt
wget -O ${OUTDIR}/cliona_varians.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc96.txt
wget -O ${OUTDIR}/halisarca_caerulea.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc94.txt
wget -O ${OUTDIR}/dendrilla_antarctica.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc91.txt
wget -O ${OUTDIR}/halichondria_panicea.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc93.txt
wget -O ${OUTDIR}/phorbas_areolatus.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc95.txt
wget -O ${OUTDIR}/phorbas_tenacior.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc98.txt
wget -O ${OUTDIR}/geodia_atlantica.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc92.txt
wget -O ${OUTDIR}/spongia_officinalis.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc101.txt
wget -O ${OUTDIR}/ircinia_fasciculata.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc99.txt
wget -O ${OUTDIR}/kirkpatrickia_variolosa.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc105.txt
wget -O ${OUTDIR}/haliclona_tubifera.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc104.txt
wget -O ${OUTDIR}/stylissa_carteri.mtdna.cds.fna https://ars.els-cdn.com/content/image/1-s2.0-S1055790320302839-mmc90.txt

src="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/features/mitogenome-retrieval-and-annotation/bin/translate_mtdna_cds_fna.py"
for mtdna_fna in `find /media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/data/external/plese_et_al_2021 -name "*.mtdna.cds.fna"`;do
    python $src --fasta ${mtdna_fna} --out ${mtdna_fna/.fna/.faa} --code 4
done
