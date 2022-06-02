#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=UNLIMITED
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.busco_abinitio_metazoa_references.err
#SBATCH --output=logs/%J.busco_abinitio_metazoa_references.out

# Retrieve positive controls from NCBI
# cd $HOME/marine_drugs/marine_drugs/data/external/
# Amphimedon Queenslandica
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/090/795/GCF_000090795.1_v1.0/GCF_000090795.1_v1.0_genomic.fna.gz
# gzip -d GCF_000090795.1_v1.0_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/292/275/GCA_016292275.1_UQ_AmQuee_3/GCA_016292275.1_UQ_AmQuee_3_genomic.fna.gz
# gzip -d GCA_016292275.1_UQ_AmQuee_3_genomic.fna.gz
# Aplysina aerophoba
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/275/565/GCA_900275565.1_Aplysina16A_polished_assembly/GCA_900275565.1_Aplysina16A_polished_assembly_genomic.fna.gz
# gzip -d GCA_900275565.1_Aplysina16A_polished_assembly_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/275/575/GCA_900275575.1_Aplysina21_polished_assembly/GCA_900275575.1_Aplysina21_polished_assembly_genomic.fna.gz
# gzip -d GCA_900275575.1_Aplysina21_polished_assembly_genomic.fna.gz
# Ephydatia muelleri
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/339/895/GCA_013339895.1_Emu_genome_v1/GCA_013339895.1_Emu_genome_v1_genomic.fna.gz
# gzip -d GCA_013339895.1_Emu_genome_v1_genomic.fna.gz


# Path to directories containing eukaryota.fna fasta files.
OUTDIR="$HOME/marine_drugs/marine_drugs/data/interim/host-annotation/busco-marker-assessment"
if [ ! -d $OUTDIR ];
then 
    mkdir -p $OUTDIR
    echo "Created ${OUTDIR}"
fi
# Ab initio output directory
cpu=50
# Now perform search on positive control (amphimedon queenslandica NCBI genome)
# full path to genome: "${HOME}/marine_drugs/marine_drugs/data/external/GCF_000090795.1_v1.0_genomic.fna"
workdir="${HOME}/marine_drugs/marine_drugs/data/external"
controls=("GCF_000090795.1_v1.0_genomic.fna" "GCA_016292275.1_UQ_AmQuee_3_genomic.fna" "GCA_900275565.1_Aplysina16A_polished_assembly_genomic.fna" "GCA_900275575.1_Aplysina21_polished_assembly_genomic.fna" "GCA_013339895.1_Emu_genome_v1_genomic.fna")
for control in ${controls[@]};do
    docker run --rm -u $(id -u) -v $workdir:/busco_wd \
            ezlabgva/busco:v5.beta.1_cv1 \
            busco -m genome \
                --lineage_dataset metazoa_odb10 \
                --in /busco_wd/$control \
                --out ${control/.fna/}_busco_genome_metazoa_odb10 \
                --cpu $cpu
    mv ${workdir}/${control/.fna/}_busco_genome_metazoa_odb10 ${OUTDIR}/.
    # Search all eukaryota
    # docker run --rm -u $(id -u) -v $workdir:/busco_wd \
    #         ezlabgva/busco:v5.beta.1_cv1 \
    #         busco -m genome \
    #             --auto-lineage-euk \
    #             --in /busco_wd/$control \
    #             --out ${control/.fna/}_busco_genome_auto_lineage_euk \
    #             --cpu $cpu
    # mv ${workdir}/${control/.fna/}_busco_genome_auto_lineage_euk ${OUTDIR}/.
done