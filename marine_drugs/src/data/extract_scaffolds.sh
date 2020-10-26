#!/usr/bin/env bash

tarballs_directory="/media/external2/Sponge_assemblies"
outdir="$HOME/marine_drugs/marine_drugs/data/raw/assemblies"
for src_scaffolds in `grep --no-filename "_Spades_output/scaffolds.fasta" ${outdir}/*_spades_tarball_listing.txt`;do
    # Example directory name: FL2015_30_Spades_output/scaffolds.fasta
    sponge_dirname=$(dirname ${src_scaffolds})
    sponge=${sponge_dirname/_Spades_output/}
    # Get full path to tarball
    tarball_name="$(dirname ${src_scaffolds}).tar.gz"
    tarball_path="${tarballs_directory}/${tarball_name}"
    # Now finish destination path
    dest_scaffolds="$(basename ${src_scaffolds})"
    # print paths for housekeeping
    echo "sponge: ${sponge} -> tarball path: ${tarball_path} -> output path: ${outdir}/${sponge}.fasta"
    # Good idea to check naming transforms with command below before running script.
    # tar -t -v --strip-components=1 --show-transformed-names  --file=${tarball_path} ${src_scaffolds}
    if [ ! -f ${outdir}/${sponge}.fasta ];
    then
        tar --keep-old-files --strip-components=1 --transform="s/scaffolds/${sponge}/" --file=${tarball_path} --directory ${outdir} --get ${src_scaffolds} &
    else
        echo "${outdir}/${sponge}.fasta already exists. skipping..."
    fi
done

wait