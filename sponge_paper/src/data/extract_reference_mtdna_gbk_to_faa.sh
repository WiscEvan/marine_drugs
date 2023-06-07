#!/usr/bin/env bash
# Retrieve CDS found within accession.gbk and write to accession.faa
# This is performed after running get_lavrov_et_al_mtdna.sh
# and manually curating the gbk_organisms.tsv files using grep and some editing
# i.e. grep "ORGANISM" *.gbk > gbk_organisms.tsv
# edit gbk_organisms.tsv to be tab-delimited and with headers Accession and Organism.
# NOTE: translate GBK accession using gbk_organisms.tsv 
# (file should be found in respective directory with headers:
# Accession and Organism columns.
script="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/sponge_paper/sponge_paper/src/data/extract_seqs_from_gbk.py"
# gbk directory: lavrov_et_al
gbk_dir="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/sponge_paper/sponge_paper/data/external/lavrov_et_al"
gbk_orgs="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/sponge_paper/sponge_paper/data/external/lavrov_et_al/gbk_organism.tsv"
gbks=(`find $gbk_dir -name "*.gbk"`)
python $script \
    --gbk ${gbks[@]} \
    --gbk-orgs $gbk_orgs \
    --outdir $gbk_dir