#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=UNLIMITED
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.busco_abinitio_proteins_metazoa_eukaryota.err
#SBATCH --output=logs/%J.busco_abinitio_proteins_metazoa_eukaryota.out

# Path to directories containing eukaryota.fna fasta files.
PREDICTIONS_DIR="$HOME/marine_drugs/marine_drugs/data/interim/host-annotation/gene-calling/ab-initio"
OUTDIR="$HOME/marine_drugs/marine_drugs/data/interim/host-annotation/augustus-ab-initio-proteins-marker-assessment"
if [ ! -d $OUTDIR ];
then
    mkdir -p $OUTDIR
    echo "Created ${OUTDIR}"
fi
# AUGUSTUS predictions output directory
cpu=50

for prediction in `ls ${PREDICTIONS_DIR}/*.abinitio.predictions.aa`;do
    workdir=$(dirname $prediction)
    proteins=$(basename $prediction)
    sample=${proteins/.hints.predictions.aa/}
    # Search all eukaryota
    # docker run --rm -u $(id -u) -v $workdir:/busco_wd \
    #     ezlabgva/busco:v5.beta.1_cv1 \
    #     busco -m proteins \
    #         --auto-lineage-euk \
    #         --in /busco_wd/$proteins \
    #         --out ${sample}_busco_proteins_auto_lineage_euk \
    #         --cpu $cpu
    # mv ${workdir}/${sample}_busco_proteins_auto_lineage_euk ${OUTDIR}/.
    # Now search all metazoans
    docker run --rm -u $(id -u) -v $workdir:/busco_wd \
        ezlabgva/busco:v5.beta.1_cv1 \
        busco -m proteins \
            --lineage_dataset metazoa_odb10 \
            --in /busco_wd/$proteins \
            --out ${sample}_busco_proteins_metazoa_odb10 \
            --cpu $cpu
    mv ${workdir}/${sample}_busco_proteins_metazoa_odb10 ${OUTDIR}/.
done