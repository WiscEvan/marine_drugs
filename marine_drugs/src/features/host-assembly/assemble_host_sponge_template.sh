#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.star_alignment_sample.err
#SBATCH --output=logs/%J.star_alignment_sample.out

sponge="sample"
genome="$HOME/marine_drugs/marine_drugs/data/interim/host-assembly/fastas/${sponge}/eukaryota.fna"
all_reads="$HOME/marine_drugs/marine_drugs/data/interim/rnaseq"
fwd_reads="${all_reads}/${sponge}*PE1*.fastq*"
rev_reads="${all_reads}/${sponge}*PE2*.fastq.gz"
se_reads="${all_reads}/${sponge}*SE*.fastq.gz"
genomeDir="$HOME/marine_drugs/marine_drugs/data/interim/host-assembly/star-alignments/${sponge}"

threads=10

if [ ! -d $genomeDir ]
then mkdir -p $genomeDir
else echo "$genomeDir already exists."
fi

STAR \
    --runThreadN $threads \
    --runMode genomeGenerate \
    --genomeDir $genomeDir \
    --genomeFastaFiles $genome

# Note: We are not passing in single-end reads with the paired-end reads
# as this was not being recognized by STAR
STAR \
    --runThreadN $threads \
    --genomeDir $genomeDir \
	--outFileNamePrefix "${genomeDir}/${sponge}_" \
    --readFilesIn $fwd_reads $rev_reads \
    --readFilesCommand zcat \
    --alignIntronMax 100000 \
    --outSAMtype BAM SortedByCoordinate \
    --outWigType wiggle \
    --outWigStrand Unstranded

