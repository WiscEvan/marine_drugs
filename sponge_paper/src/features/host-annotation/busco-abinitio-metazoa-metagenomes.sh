#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=UNLIMITED
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.busco_abinitio_metazoa_metagenomes.err
#SBATCH --output=logs/%J.busco_abinitio_metazoa_metagenomes.out

# Path to directories containing eukaryota.fna fasta files.
METAGENOMES="$HOME/sponge_paper/sponge_paper/data/interim/assemblies"
OUTDIR="$HOME/sponge_paper/sponge_paper/data/interim/host-annotation/busco-marker-assessment"

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
    # Now search all metazoans
    docker run --rm -u $(id -u) -v $workdir:/busco_wd \
        ezlabgva/busco:v5.beta.1_cv1 \
        busco -m genome \
            --lineage_dataset metazoa_odb10 \
            --in /busco_wd/${sample} \
            --out ${sample/.filtered.fna/}_busco_genome_metazoa_odb10 \
            --cpu $cpu
    mv "${workdir}/*busco*" "${OUTDIR}/."
done