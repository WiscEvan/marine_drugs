#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.bacteria_binning_script_generation.err
#SBATCH --output=logs/%J.bacteria_binning_script_generation.out

template="$HOME/sponge_paper/sponge_paper/src/features/sponge_binning_template.sh"
features="$HOME/sponge_paper/sponge_paper/src/features"
metagenomes=`ls $HOME/sponge_paper/sponge_paper/data/raw/assemblies/*.fasta`

for mg in $metagenomes;do
    sponge=$(basename ${mg/.fasta/})
    script="${features}/${sponge}_bacteria_binning.sh"
    # Copy sponge_binning_template.sh to FL2015_4_bacteria_binning.sh
    cp ${template} ${script}
    
    sed -i "s,sample,${sponge},g" ${script}

    echo "Wrote: ${script}"
done
