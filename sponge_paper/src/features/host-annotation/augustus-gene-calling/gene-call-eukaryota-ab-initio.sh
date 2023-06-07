#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.augustus_gene_calling_ab_initio.err
#SBATCH --output=logs/%J.augustus_gene_calling_ab_initio.out

EUKARYA="/home/evan/sponge_paper/sponge_paper/data/interim/host-assembly/fastas"
OUTDIR="/home/evan/sponge_paper/sponge_paper/data/interim/host-annotation/gene-calling"
# species identifier for amphimedon queenslandica
SPECIES="amphimedon"

for eukaryota in `find $EUKARYA -name "eukaryota.fna"`;do
    INDIR=$(dirname $eukaryota)
    sponge=$(basename $(dirname $eukaryota))
    abinitio_output="${OUTDIR}/${sponge}.abinitio.gff"
    # 1. Perform ab initio predictions with Augustus
    if [ ! -f $abinitio_output ]
    then 
        docker run \
        --volume $INDIR:/input:ro \
        --user=$(id -u):$(id -g) \
        --rm \
        --detach=false \
        augustus:latest \
            augustus \
                --species=$SPECIES \
                --sample=${sponge} \
                --genemodel=partial \
                --protein=on \
                --codingseq=on \
                /input/eukaryota.fna > $abinitio_output
    else
        echo "${abinitio_output} already exists. skipping..."
    fi
done