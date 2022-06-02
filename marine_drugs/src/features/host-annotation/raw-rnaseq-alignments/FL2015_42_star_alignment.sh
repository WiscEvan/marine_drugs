#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.star_alignment_FL2015_42.err
#SBATCH --output=logs/%J.star_alignment_FL2015_42.out

sponge="FL2015_42"
genome="$HOME/marine_drugs/marine_drugs/data/interim/host-assembly/fastas/${sponge}/eukaryota.fna"
all_reads="$HOME/marine_drugs/marine_drugs/data/raw/rnaseq"
# Reads have FL2015-* rather than FL2015_* naming, so we replace the underscore ('_') with a dash ('-')
# We also need to truncate any FL2014 names to FL20
rna_name=${sponge/\1\4/}
rna_name=${rna_name/\_/\-}
prefix_num=$(echo $rna_name | cut -f1 -d"-")
if [ $prefix_num == "FL20" ]
then
    echo "Reads for FL2014 do not require padding. skipping..."
    suffix_num=$(echo $rna_name | cut -f2 -d"-")
else
    echo "Reads for FL2015 require padding. padding..."
    suffix_num=$(echo $rna_name | cut -f2 -d"-"  | xargs printf "%02d")
fi
rna_name="${prefix_num}-${suffix_num}"
fwd_reads="${all_reads}/${rna_name}_*R1*.fastq.gz"
rev_reads="${all_reads}/${rna_name}_*R2*.fastq.gz"
# se_reads="${all_reads}/${rna_name}_*SE*.fastq.gz"
genomeDir="$HOME/marine_drugs/marine_drugs/data/interim/host-assembly/star-alignments/raw-alignments/${sponge}"

threads=24

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

