#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=UNLIMITED
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.busco_abinitio.err
#SBATCH --output=logs/%J.busco_abinitio.out

# Path to directories containing eukaryota.fna fasta files.
EUKARYA="$HOME/marine_drugs/marine_drugs/data/interim/host-assembly/fastas"
# Ab initio output directory
cpu=30

for assembly in `ls ${EUKARYA}/*/eukaryota.fna`;do
    workdir=$(dirname $assembly)
    sample=$(basename $workdir)
    # Search all eukaryota
    docker run --rm -u $(id -u) -v $workdir:/busco_wd \
        ezlabgva/busco:v5.beta.1_cv1 \
        busco -m genome \
            --auto-lineage-euk \
            --in /busco_wd/eukaryota.fna \
            --out ${sample}_busco_genome_auto_lineage_euk \
            --cpu $cpu
    # Now search all metazoans
    docker run --rm -u $(id -u) -v $workdir:/busco_wd \
        ezlabgva/busco:v5.beta.1_cv1 \
        busco -m genome \
            --lineage_dataset metazoa_odb10 \
            --in /busco_wd/eukaryota.fna \
            --out ${sample}_busco_genome_metazoa_odb10 \
            --cpu $cpu
done