#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.taxonomy_script_generation.err
#SBATCH --output=logs/%J.taxonomy_script_generation.out

template="$HOME/sponge_paper/sponge_paper/src/features/taxonomy_template.sh"
features="$HOME/sponge_paper/sponge_paper/src/features"
metagenomes=`ls $HOME/sponge_paper/sponge_paper/data/raw/assemblies/*.fasta`
interim="$HOME/sponge_paper/sponge_paper/data/interim"
filtered_assemblies="${interim}/assemblies"
binnings="${interim}/binning"

for mg in $metagenomes;do
    metagenome=$(basename ${mg/.fasta/})
    get_taxon_script="${features}/get_${metagenome}_taxonomy.sh"
    # Copy template.mf to FL2015_4.sh
    cp ${template} ${get_taxon_script}
    # edit FL2015_4.mf inputs and outputs filepaths
    sed -i "s,sample,${metagenome},g" ${get_taxon_script}
    
    # FILTERED
    filtered_assembly="${filtered_assemblies}/${metagenome}.filtered.fna"
    sed -i "s,metagenome.filtered.fna,${filtered_assembly},g" ${get_taxon_script}
    # ORFS
    orfs="${binnings}/${metagenome}.orfs."
    sed -i "s,orfs\.,${orfs},g" ${get_taxon_script}
    # TAXONOMY
    taxonomy="${binnings}/${metagenome}.taxonomy."
    sed -i "s,taxonomy\.,${taxonomy},g" ${get_taxon_script}
    # BLASTP
    blastp="${binnings}/${metagenome}.blastp."
    sed -i "s,blastp\.,${blastp},g" ${get_taxon_script}
    # LCA
    lca="${binnings}/${metagenome}.lca."
    sed -i "s,lca\.,${lca},g" ${get_taxon_script}
    echo "Wrote: ${get_taxon_script}"
done


# autometa-taxonomy --ncbi /mnt/autometa_databases/ --assembly metagenome.filtered.fna --nucl-orfs orfs.fna --prot-orfs orfs.faa --blast blastp.tsv --lca lca.tsv --method majority_vote taxonomy.tsv
