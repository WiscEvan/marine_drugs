#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.bacteria_binning_script_generation.err
#SBATCH --output=logs/%J.bacteria_binning_script_generation.out

template="$HOME/marine_drugs/marine_drugs/src/features/sponge_binning_template.sh"
features="$HOME/marine_drugs/marine_drugs/src/features"
metagenomes=`ls $HOME/marine_drugs/marine_drugs/data/raw/assemblies/*.fasta`

for mg in $metagenomes;do
    sponge=$(basename ${mg/.fasta/})
    script="${features}/${sponge}_bacteria_binning.sh"
    # Copy sponge_binning_template.sh to FL2015_4_bacteria_binning.sh
    cp ${template} ${script}
    
    sed -i "s,sample,${sponge},g" ${script}

    echo "Wrote: ${script}"
done
