#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.augustus_gene_calling_with_hints.err
#SBATCH --output=logs/%J.augustus_gene_calling_with_hints.out


# NOTE: Before running this script...
# Alignments must first be conducted using STAR as outlined below:
# ================================================================
# STAR \
#     --runThreadN $threads \
#     --runMode genomeGenerate \
#     --genomeDir $genomeDir \
#     --genomeFastaFiles $genome

# # Note: We are not passing in single-end reads with the paired-end reads
# # as this was not being recognized by STAR
# STAR \
#     --runThreadN $threads \
#     --genomeDir $genomeDir \
# 	  --outFileNamePrefix "${genomeDir}/${sponge}_" \
#     --readFilesIn $fwd_reads $rev_reads \
#     --readFilesCommand zcat \
#     --alignIntronMax 100000 \
#     --outSAMtype BAM SortedByCoordinate \
#     --outWigType wiggle \
#     --outWigStrand Unstranded
# ================================================================


# Path to directories containing eukaryota.fna fasta files.
EUKARYA="$HOME/sponge_paper/sponge_paper/data/interim/host-assembly/fastas"
# Ab initio output directory
OUTDIR="$HOME/sponge_paper/sponge_paper/data/interim/host-annotation/gene-calling"
# STAR Alignments directory
ALNDIR="$HOME/sponge_paper/sponge_paper/data/interim/host-annotation/star-alignments/seqycleaned-rnaseq"
# species identifier for amphimedon queenslandica
SPECIES="amphimedon"

for eukaryota in `find $EUKARYA -name "eukaryota.fna"`;do
    INDIR=$(dirname $eukaryota)
    sponge=$(basename $(dirname $eukaryota))
    if [ $sponge == "FL2015_8" ] || [ $sponge == "FL2015_30" ]
    then
        echo "Skipping ${sponge}. No RNA seq available to provide hints..." 
        continue
    fi
    # Aligned.sortedByCoord.out.bam written by STAR
    # Signal.Unique.str1.out.wig written by STAR
    # Generate hints for the introns.
    SPONGE_ALN_DIR="${ALNDIR}/${sponge}"
    introns_hints_filename="${sponge}.hints.introns.gff"

    # 2c. Create a directory that contains hints and concatenate
    HINTS_PRED_DIR="${OUTDIR}/${sponge}"
    if [ ! -d $HINTS_PRED_DIR ]
    then mkdir -p $HINTS_PRED_DIR
    fi
    
    # 2a. Generate intron hints for Augustus predictions from RNAseq alignments
    introns_hints="${HINTS_PRED_DIR}/${introns_hints_filename}"
    if [ ! -f $introns_hints ]
    then
        # --in follows convention from STAR parameter --outFileNamePrefix ${genomeDir}/${sponge}_'
        docker run \
            --volume $SPONGE_ALN_DIR:/sample:rw \
            --user=$(id -u):$(id -g) \
            --rm \
            --detach=false \
            augustus:latest \
                bin/bam2hints --intronsonly --in=/sample/${sponge}_Aligned.sortedByCoord.out.bam --out=/sample/${introns_hints_filename}
        mv "${SPONGE_ALN_DIR}/${introns_hints_filename}" $introns_hints
    else
        echo "${SPONGE_ALN_DIR}/${introns_hints_filename} exists... skipping bam2hints"
    fi
    # Now generate the exon part hints
    # Note: exon parts are handled differently than 'exons' in Augustus
    # For more information see the extrinsic/extrinsic.cfg file in the Augustus repository
    # So I am specifying 'exonparts' to avoid confusion
    exonparts_hints="${HINTS_PRED_DIR}/${sponge}.hints.exonparts.gff"

    # NOTE: sample name follows convention from STAR parameter --outFileNamePrefix ${genomeDir}/${sponge}_'
    # 2b. Generate exon hints for Augustus predictions from RNAseq alignments
    if [ ! -f $exonparts_hints ]
    then
        # NOTE: We must use /bin/bash so we can supply two commands within the container! 
        docker run \
            --volume $SPONGE_ALN_DIR:/sample:rw \
            --user=$(id -u):$(id -g) \
            --rm \
            --detach=false \
            augustus:latest /bin/bash -c "cat /sample/${sponge}_Signal.Unique.str1.out.wig \
                | \
                ./scripts/wig2hints.pl \
                    --width=10 \
                    --margin=10 \
                    --minthresh=2 \
                    --minscore=4 \
                    --prune=0.1 \
                    --src=W \
                    --type=ep \
                    --radius=4.5 \
                    --pri=4 \
                    --strand='.'" > $exonparts_hints
    else
        echo "${exonparts_hints} exists... skipping wig2hints exonparts generation"
    fi

    HINTS="${HINTS_PRED_DIR}/${sponge}.hints.gff"
    cat $introns_hints $exonparts_hints > $HINTS
    if [ ! -f $HINTS ]
    then
        echo "Missing hints file for ${sponge}: hints -> $HINTS"
        continue
    fi
    # 3. Perform Augustus predictions with the generated hints.
    HINTS_PREDICTIONS="${OUTDIR}/${sponge}.hints.predictions.gff"
    
    # NOTE: Copy extrinsic from Augusutus path below and write to $EXTRINSIC_CONFIG
    # AUGUSTUS_CONFIG="$HOME/Augustus/config/extrinsic/extrinsic.M.RM.E.W.cfg"
    EXTRINSIC_CONFIG="${OUTDIR}/${sponge}.extrinsic.cfg"

    # Augustus dev recommends switching the UTR flag to "on" because
    # RNA-Seq covers UTR as well. 
    # With --UTR=off, the exonpart hints would be (mis)interpreted to be hints for coding parts of exons.
    # See: http://augustus.gobics.de/binaries/readme.rnaseq.html
    if [ ! -f $HINTS_PREDICTIONS ]
    then
        docker run \
            --rm \
            --detach=false \
            --user=$(id -u):$(id -g) \
            --volume $INDIR:/fasta:ro \
            --volume $OUTDIR:/predictions:rw \
            augustus:latest \
                augustus \
                    --species=$SPECIES \
                    --alternatives-from-evidence=true \
                    --extrinsicCfgFile=/predictions/$(basename $EXTRINSIC_CONFIG) \
                    --allow_hinted_splicesites=atac \
                    --softmasking=on \
                    --protein=on \
                    --codingseq=on \
                    --UTR=on \
                    --hintsfile=/predictions/${sponge}/$(basename $HINTS) \
                    /fasta/eukaryota.fna > $HINTS_PREDICTIONS
    else
        echo "${HINTS_PREDICTIONS} exists... skipping"
    fi
done