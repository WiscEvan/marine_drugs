#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.antismash_sample.err
#SBATCH --output=logs/%J.antismash_sample.out

# https://phdops.kblin.org/2015-running-antismash-standalone-from-docker.html
ANTISMASH="/home/evan/bin/run_antismash"
# Note: `run_antismash` wrapper has been altered to run the latest version of antismash docker image.

cpus=35
assemblies="$HOME/sponge_paper/sponge_paper/data/raw/assemblies"
sponge="sample"
fasta="${assemblies}/${sponge}.fasta"
if [ ! -f $fasta ]
then echo "could not find $fasta. skipping..."
else
    output="$HOME/sponge_paper/sponge_paper/data/interim/bgcs/${sponge}"
    zipfile="${output}/*.zip"
    if [ -f $zipfile ]
    then echo "${sponge} antiSMASH results already exist, skipping."
    else
        $ANTISMASH $fasta $output \
            --genefinding-tool prodigal-m \
            --fullhmmer \
            --cf-create-clusters \
            --clusterhmmer \
            --cb-general \
            --cb-subclusters \
            --cb-knownclusters \
            --pfam2go \
            --cpus $cpus
    fi
fi
