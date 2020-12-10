#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.augustus_gene_calling.err
#SBATCH --output=logs/%J.augustus_gene_calling.out

EUKARYA="/home/evan/marine_drugs/marine_drugs/data/interim/host-assembly/fastas"
OUTDIR="/home/evan/marine_drugs/marine_drugs/data/interim/host-annotation/gene-calling"
# species identifier for amphimedon queenslandica
SPECIES="amphimedon"

for eukaryota in `find $EUKARYA -name "eukaryota.fna"`;do
    INDIR=$(dirname $eukaryota)
    sponge=$(basename $(dirname $eukaryota))
    output="${OUTDIR}/${sponge}.abinitio.txt"
    if [ ! -f $output ]
    then docker run \
        --volume $INDIR:/input:ro \
        --user=$(id -u):$(id -g) \
        --rm \
        --detach=false \
        augustus:latest \
            augustus --species=$SPECIES --sample=${sponge} --genemodel=partial /input/eukaryota.fna > $output
    else echo "${output} already exists. skipping..."
    fi
    
    # Aligned.sortedByCoord.out.bam written by STAR
    # Signal.Unique.str1.out.wig written by STAR
    # echo "chr2L   23513712" > chrom.sizes # size of chromosomes
    # mv Signal.Unique.str1.out.wig rnaseq.wig    # rename for convenience
    # mv Aligned.sortedByCoord.out.bam rnaseq.bam 
    # wigToBigWig rnaseq.wig chrom.sizes rnaseq.bw
    # bamtools index -in rnaseq.bam # index required by UCSC browser
    # # Now prepare files for visualization.
    # sponge_gff="$OUTDIR/${sponge}.gff"
    # echo "browser hide all" > customtrack
    # echo "track name=\"STAR RNA-Seq alignments\" type=bam visibility=4 bigDataUrl=http://hgwdev.cse.ucsc.edu/~mario/tutorial2015/results/rnaseq.bam" >> customtrack
    # echo "track name=\"STAR coverage\" type=bigWig visibility=2 bigDataUrl=http://hgwdev.cse.ucsc.edu/~mario/tutorial2015/results/rnaseq.bw" >> customtrack
    # echo "track name=\"AUGUSTUS $s abinitio\"  db=dm6 visibility=3" >> customtrack
    # grep -P "AUGUSTUS\t(CDS|exon)\t" $augustus_stdout >> customtrack # use only the relevant coordinate lines

done