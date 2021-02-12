#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.star_alignment_script_generation.err
#SBATCH --output=logs/%J.star_alignment_script_generation.out

template="$HOME/marine_drugs/marine_drugs/src/features/host-assembly/raw-alignments/star-alignment-template.sh"
outdir="$HOME/marine_drugs/marine_drugs/src/features/host-assembly/raw-alignments"

sponges_eukaryotic_contigs=`ls -d /home/evan/marine_drugs/marine_drugs/data/interim/host-assembly/fastas/*/`
for sample in $sponges_eukaryotic_contigs;do
    sponge=$(basename $sample)
    # FL2015-8 does not have RNAseq data.
    if [ $sponge == "FL2015_8" ] || [ $sponge == "FL2015_30" ]
    then echo "$sponge does not have RNA seq data. Skipping..."
    else
        script="${outdir}/${sponge}_star_alignment.sh"
        # Copy $template to sponge_start_alignment.sh
        cp ${template} ${script}
        # Now edit sample inputs and outputs filepaths
        sed -i "s,sample,${sponge},g" ${script}
        echo "Wrote: ${script}"
    fi
done
