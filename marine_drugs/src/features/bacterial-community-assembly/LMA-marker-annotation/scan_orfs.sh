#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.sponge_hmmscan_taxon_specific_markers.err
#SBATCH --output=logs/%J.sponge_hmmscan_taxon_specific_markers.out



lma_sponges=($(cat lma_sponges.txt))

DATA="${HOME}/marine_drugs/marine_drugs/data/interim/binning"
OUTDIR="${HOME}/marine_drugs/marine_drugs/data/interim/LMA-marker-annotation"
CPUS=8


if [ ! -d $OUTDIR ]
then mkdir -p $OUTDIR
fi

EXTERNAL="${HOME}/marine_drugs/marine_drugs/data/external"

for sponge in ${lma_sponges[@]};do
    # hmmscan <hmmdb> <seqfile>
    seqfile="${DATA}/${sponge}.orfs.faa"
    for hmmdb in `ls ${EXTERNAL}/*/*.hmm`;do
        if [[ $hmmdb == *"Bacteria"* ]]
        then
            echo "Skipping $hmmdb"
        else
            hmmdb_basename=$(basename ${hmmdb/.hmm/})
            tblout="${OUTDIR}/${sponge}_${hmmdb_basename}.hmmscan.tsv"
            hmmscan --cpu $CPUS --tblout $tblout $hmmdb $seqfile
        fi
    done
done