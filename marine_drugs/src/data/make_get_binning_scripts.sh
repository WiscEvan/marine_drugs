#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.taxonomy_script_generation.err
#SBATCH --output=logs/%J.taxonomy_script_generation.out

template="$HOME/marine_drugs/marine_drugs/src/features/binning_template.mf"
features="$HOME/marine_drugs/marine_drugs/src/features"
metagenomes=`ls $HOME/marine_drugs/marine_drugs/data/raw/assemblies/*.fasta`
interim="$HOME/marine_drugs/marine_drugs/data/interim"
filtered_assemblies="${interim}/assemblies"
binnings="${interim}/binning"

for mg in $metagenomes;do
    metagenome=$(basename ${mg/.fasta/})
    makeflow="${features}/binning_${metagenome}.mf"
    # Copy template.mf to binning_FL2015_4.mf
    cp ${template} ${makeflow}
    
    # KMERS
    kmers="${binnings}/${metagenome}.kmers."
    sed -i "s,kmers\.,${kmers},g" ${makeflow}
    # COVERAGES
    coverages="${binnings}/${metagenome}.coverages."
    sed -i "s,coverages\.,${coverages},g" ${makeflow}
    # TAXONOMY
    taxonomy="${binnings}/${metagenome}.taxonomy."
    sed -i "s,taxonomy\.,${taxonomy},g" ${makeflow}
    # BACTERIA MARKERS
    markers="${binnings}/${metagenome}.bacteria.markers."
    sed -i "s,bacteria\.markers\.,${markers},g" ${makeflow}
    # ARCHAEA MARKERS
    markers="${binnings}/${metagenome}.archaea.markers."
    sed -i "s,archaea\.markers\.,${markers},g" ${makeflow}
    # BACTERIA BINNING
    binning="${binnings}/${metagenome}.bacteria.binning."
    sed -i "s,bacteria\.binning\.,${binning},g" ${makeflow}
    # ARCHAEA BINNING
    binning="${binnings}/${metagenome}.archaea.binning."
    sed -i "s,archaea\.binning\.,${binning},g" ${makeflow}
    # BACTERIA UNCLUSTERED RECRUITMENT
    recruitment="${binnings}/${metagenome}.bacteria.recruitment."
    sed -i "s,bacteria.unclustered_recruitment\.,${recruitment},g" ${makeflow}
    # ARCHAEA UNCLUSTERED RECRUITMENT
    recruitment="${binnings}/${metagenome}.archaea.recruitment."
    sed -i "s,archaea.unclustered_recruitment\.,${recruitment},g" ${makeflow}

    echo "Wrote: ${makeflow}"
done
