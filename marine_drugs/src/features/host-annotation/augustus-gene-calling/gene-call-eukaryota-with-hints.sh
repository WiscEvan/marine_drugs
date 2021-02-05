#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.augustus_gene_calling.err
#SBATCH --output=logs/%J.augustus_gene_calling.out

# Path to directories containing eukaryota.fna fasta files.
EUKARYA="$HOME/marine_drugs/marine_drugs/data/interim/host-assembly/fastas"
# Ab initio output directory
OUTDIR="$HOME/marine_drugs/marine_drugs/data/interim/host-annotation/gene-calling"
# STAR Alignments directory
ALNDIR="$HOME/marine_drugs/marine_drugs/data/interim/host-annotation/star-alignments/seqycleaned-rnaseq/"
# species identifier for amphimedon queenslandica
SPECIES="amphimedon"

for eukaryota in `find $EUKARYA -name "eukaryota.fna"`;do
    INDIR=$(dirname $eukaryota)
    sponge=$(basename $(dirname $eukaryota))
    abinitio_output="${OUTDIR}/${sponge}.abinitio.gff"
    # 1. Perform ab initio predictions with Augustus
    if [ ! -f $abinitio_output ]
    then 
        docker run \
        --volume $INDIR:/input:ro \
        --user=$(id -u):$(id -g) \
        --rm \
        --detach=false \
        augustus:latest \
            augustus \
                --species=$SPECIES \
                --sample=${sponge} \
                --genemodel=partial \
                --protein=on \
                --codingseq=on \
                /input/eukaryota.fna > $abinitio_output
    else echo "${abinitio_output} already exists. skipping..."
    fi
    
    if [ $sponge == "FL2015_8" ] || [ $sponge == "FL2015_30" ]
    then
        echo "Skipping ${sponge}. No RNA seq available to provide hints..." 
        continue
    fi
    # Aligned.sortedByCoord.out.bam written by STAR
    # Signal.Unique.str1.out.wig written by STAR
    # Generate hints for the introns.
    SPONGE_ALN_DIR="${ALNDIR}/${sponge}"
    intron_hints_filename="${sponge}.hints.introns.gff"

    # 2a. Generate intron hints for Augustus predictions from RNAseq alignments
    docker run \
        --volume $SPONGE_ALN_DIR:/sample:rw \
        --user=$(id -u):$(id -g) \
        --rm \
        --detach=false \
        augustus:latest \
            bin/bam2hints --intronsonly --in=/sample/${sponge}_Aligned.sortedByCoord.out.bam --out=/sample/${intron_hints_filename}
    
    # Now generate the exon part hints
    # Note: exon parts are handled differently than 'exons' in Augustus
    # For more information see the extrinsic/extrinsic.cfg file in the Augustus repository
    # So I am specifying 'exonparts' to avoid confusion
    ep_hints="${SPONGE_ALN_DIR}/${sponge}.hints.exonparts.gff"

    # 2b. Generate exon hints for Augustus predictions from RNAseq alignments
    docker run \
        --volume $SPONGE_ALN_DIR:/sample:rw \
        --user=$(id -u):$(id -g) \
        --rm \
        --detach=false \
        augustus:latest \
            cat /sample/${sponge}_Signal.Unique.str1.out.wig \
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
                --strand="." > $ep_hints

    # 2c. Create a directory that contains hints and concatenate
    HINTS_PRED_DIR="${OUTDIR}/${sponge}"
    if [ ! -d $HINTS_PRED_DIR ]
    then mkdir -p $HINTS_PRED_DIR
    fi
    HINTS="${HINTS_PRED_DIR}/${sponge}.hints.gff"
    cat "${SPONGE_ALN_DIR}/${intron_hints_filename}" $ep_hints > $HINTS
    if [ ! -f $HINTS ]
    then
        echo "Missing hints file for ${sponge}: hints -> $HINTS"
        continue
    fi

    # 3. Perform Augustus predictions with the generated hints.
    HINTS_PREDICTIONS="${OUTDIR}/${sponge}.hints.predictions.gff"
    # AUGUSTUS_CONFIG="$HOME/Augustus/config/extrinsic/extrinsic.M.RM.E.W.cfg"
    # EXTRINSIC_CONFIG="${OUTDIR}/${sponge}.extrinsic.cfg"

    # Augustus dev recommends switching the UTR flag to "on" because
    # RNA-Seq covers UTR as well. 
    # With --UTR=off, the exonpart hints would be (mis)interpreted to be hints for coding parts of exons.
    # See: http://augustus.gobics.de/binaries/readme.rnaseq.html
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
                --extrinsicCfgFile=extrinsic.cfg \
                --allow_hinted_splicesites=atac \
                --softmasking=on \
                --protein=on \
                --codingseq=on \
                --UTR=on \
                --hintsfile=/predictions/$HINTS \
                /fasta/eukaryota.fna  > $HINTS_PREDICTIONS
                # --extrinsicCfgFile=$EXTRINSIC_CONFIG  > $HINTS_PREDICTIONS
done