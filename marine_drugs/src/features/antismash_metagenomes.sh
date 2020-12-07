#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.antismash_metagenomes.err
#SBATCH --output=logs/%J.antismash_metagenomes.out


cpus=35
assemblies="$HOME/marine_drugs/marine_drugs/data/raw/assemblies"
for metagenome in `ls ${assemblies}/*.fasta`;do
    sponge=$(basename ${metagenome/.fasta/})
    output="$HOME/marine_drugs/marine_drugs/data/interim/bgcs/${sponge}"
    # https://phdops.kblin.org/2015-running-antismash-standalone-from-docker.html
    run_antismash $metagenome $output --genefinding-tool prodigal-m --cf-create-clusters --clusterhmmer --cb-general --cb-subclusters --cb-knownclusters --pfam2go --cpus $cpus
done
