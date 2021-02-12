#!usr/bin/env bash


cd $HOME/Autometa


declare -a assemblies
assemblies=(FL2014_9 FL2015_43 FL2015_5 FL2015_9)
outdir="$HOME/marine_drugs/marine_drugs/data/interim/binning"
samples="$HOME/marine_drugs/marine_drugs/data/raw/assemblies"
for assembly in ${assemblies[@]};do
    # ../../data/raw/assemblies/FL2015_43.fasta
    sample="${samples}/${assembly}.fasta"
    nucls_out="${outdir}/${assembly}.orfs.fna"
    prots_out="${outdir}/${assembly}.orfs.faa"
    python -m autometa.common.external.prodigal --parallel --cpus 35 $sample $nucls_out $prots_out
done