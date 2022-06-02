#!/usr/bin/env bash
# Construct submission and shell files for HTCondor for blastp analysis.

shell_template="$HOME/marine_drugs/marine_drugs/src/features/condor_blastp_template.sh"
sub_template="$HOME/marine_drugs/marine_drugs/src/features/condor_blastp_template.sub"
binning="$HOME/marine_drugs/marine_drugs/data/interim/binning"
features="$HOME/marine_drugs/marine_drugs/src/features"
diamond="$HOME/marine_drugs/marine_drugs/src/features/diamond"

echo "Writing submission files to ${features}"

for filepath in `ls ${binning}/*.orfs.faa`;do
    sample=$(basename ${filepath/.orfs.faa/})
    # 1. Copy templates to sample.sub and sample.sh
    # e.g. blastp_FL2014_3.sh
    shellfile="${features}/condor_blastp_${sample}.sh"
    cp $shell_template $shellfile
    
    # e.g. blastp_FL2014_3.sub
    subfile="${features}/condor_blastp_${sample}.sub"
    cp $sub_template $subfile
    
    # 2. Alter file paths and parameters to match sample
    # 2.1 ORFs filepaths preparation
    # 2.1.1 ORFs filepath in shell executable
    # We just want the file name here because this
    # will be the path on the execute node
    filename=$(basename $filepath)
    sed -i "s,orfs.faa,${filename},g" $shellfile
    # 2.1.1 ORFs filepath in condor submission file
    sed -i "s,orfs.faa,${filepath},g" $subfile
    # 2.3 Replace blastp output name to match sample
    # replace blastp.tsv in shell file
    sed -i "s,blastp\.,${sample}.blastp.,g" $shellfile
    # 2.4 Rename <sample> placeholder with sample name
    # replace sponge in shell file
    sed -i "s,sample,${sample},g" $shellfile
    # 2.6
    # replace executable in sub file.
    sed -i "s,sample.sh,${shellfile},g" $subfile
    # replace job title in sub file.
    sed -i "s,sample,${sample},g" $subfile

    # Edit path to diamond executable
    sed -i "s,diamond,${diamond},g" $subfile

    echo "condor executable: ${shellfile}"
    echo "condor_submit ${subfile}"
done

echo "done"