#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.antismash_metagenomes.err
#SBATCH --output=logs/%J.antismash_metagenomes.out

# https://phdops.kblin.org/2015-running-antismash-standalone-from-docker.html
ANTISMASH="/home/evan/bin/run_antismash"
# Note: `run_antismash` wrapper has been altered to run the latest version of antismash docker image.

# Split up sponges b/w deep-thought and kwanlab server.
# First 6 sponges
sponges=(FL2014_3 FL2014_9 FL2015_30 FL2015_34 FL2015_37 FL2015_42)
# Last 6 sponges
# sponges=(FL2015_43 FL2015_44 FL2015_4 FL2015_5 FL2015_8 FL2015_9)

cpus=35
assemblies="$HOME/marine_drugs/marine_drugs/data/raw/assemblies"
for sponge in ${sponges[@]};do
    fasta="${assemblies}/${sponge}.fasta"
    if [ ! -f $fasta ]
    then echo "could not find $fasta. skipping..."
    else
        output="$HOME/marine_drugs/marine_drugs/data/interim/bgcs/${sponge}"
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
done
