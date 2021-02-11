#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=UNLIMITED
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.busco_abinitio_metazoa_metagenomes.err
#SBATCH --output=logs/%J.busco_abinitio_metazoa_metagenomes.out

# Path to directories containing eukaryota.fna fasta files.
METAGENOMES="$HOME/marine_drugs/marine_drugs/data/interim/assemblies"
OUTDIR="$HOME/marine_drugs/marine_drugs/data/interim/host-annotation/metagenome-marker-assessment"

if [ ! -d $OUTDIR ];
then 
    mkdir -p $OUTDIR
    echo "Created ${OUTDIR}"
fi
# Ab initio output directory
cpu=50

for metagenome in `ls ${METAGENOMES}/*.fna`;do
    workdir=$(dirname $metagenome)
    # FL2014_3.filtered.fna
    sample=$(basename $metagenome)
    # Search all eukaryota
    docker run --rm -u $(id -u) -v $workdir:/busco_wd \
        ezlabgva/busco:v5.beta.1_cv1 \
        busco -m genome \
            --auto-lineage-euk \
            --in /busco_wd/${sample} \
            --out ${sample/.filtered.fna/}_busco_genome_auto_lineage_euk \
            --cpu $cpu
    echo mv "${workdir}/*busco*" "${OUTDIR}/."
    # Now search all metazoans
    echo docker run --rm -u $(id -u) -v $workdir:/busco_wd \
        ezlabgva/busco:v5.beta.1_cv1 \
        busco -m genome \
            --lineage_dataset metazoa_odb10 \
            --in /busco_wd/${sample} \
            --out ${sample/.filtered.fna/}_busco_genome_metazoa_odb10 \
            --cpu $cpu
    echo mv "${workdir}/*busco*" "${OUTDIR}/."
done