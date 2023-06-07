#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=UNLIMITED
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.busco_abinitio_metazoa_references.err
#SBATCH --output=logs/%J.busco_abinitio_metazoa_references.out

# NOTE: Only Amphimedon Queenslandica has a reference protein set available
# Retrieve positive control from NCBI
# cd $HOME/sponge_paper/sponge_paper/data/external/
# Amphimedon Queenslandica
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/090/795/GCF_000090795.1_v1.0/GCF_000090795.1_v1.0_protein.faa.gz
# gzip -d GCF_000090795.1_v1.0_protein.faa.gz

# Path to directories containing eukaryota.fna fasta files.
OUTDIR="$HOME/sponge_paper/sponge_paper/data/interim/host-annotation/busco-marker-assessment"
if [ ! -d $OUTDIR ];
then 
    mkdir -p $OUTDIR
    echo "Created ${OUTDIR}"
fi
# Ab initio output directory
cpu=50
# Now perform search on positive control (amphimedon queenslandica NCBI proteins)
# full path to proteins: "${HOME}/sponge_paper/sponge_paper/data/external/GCF_000090795.1_v1.0_protein.faa"
workdir="${HOME}/sponge_paper/sponge_paper/data/external"
control="GCF_000090795.1_v1.0_protein.faa"
docker run --rm -u $(id -u) -v $workdir:/busco_wd \
        ezlabgva/busco:v5.beta.1_cv1 \
        busco -m proteins \
            --lineage_dataset metazoa_odb10 \
            --in /busco_wd/${control} \
            --out ${control/.faa/}_busco_proteins_metazoa_odb10 \
            --cpu $cpu
mv "${workdir}/${control/.fna/}_busco_proteins_metazoa_odb10" "${OUTDIR}/."