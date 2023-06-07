#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.augustus_extract_fastas.err
#SBATCH --output=logs/%J.augustus_extract_fastas.out


# Path to directory containing all augustus predictions following wildcard *.hints.predictions.gff
PREDICTION_DIR="$HOME/sponge_paper/sponge_paper/data/interim/host-annotation/gene-calling"

# Path to getAnnoFasta.pl from Augustus/scripts
SCRIPT="$HOME/Augustus/scripts/getAnnoFasta.pl"

# Must follow file tree structure ${FASTAS}/<sponge>/eukaryota.fna
FASTAS="$HOME/sponge_paper/sponge_paper/data/interim/host-assembly/fastas"

for prediction in `ls ${PREDICTION_DIR}/*.hints.predictions.gff`;do
    # prediction -> e.g. FL2015_37.hints.predictions.gff
    sponge=$(basename ${prediction/.hints.predictions.gff/})
    seqfile="${FASTAS}/${sponge}/eukaryota.fna"
    $SCRIPT --seqfile=$seqfile $prediction
done