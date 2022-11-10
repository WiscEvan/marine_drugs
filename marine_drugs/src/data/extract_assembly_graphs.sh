#!/usr/bin/env bash

tarballs_directory="/media/external2/Sponge_assemblies"
outdir="$HOME/marine_drugs/marine_drugs/data/raw/assemblies"
for src_asm_graph in `grep --no-filename "_Spades_output/assembly_graph_with_scaffolds.gfa" ${outdir}/*_spades_tarball_listing.txt`;do
    # Example directory name: FL2015_30_Spades_output/assembly_graph.fastg
    # Example directory name: FL2015_30_Spades_output/assembly_graph_with_scaffolds.gfa
    sponge_dirname=$(dirname ${src_asm_graph})
    sponge=${sponge_dirname/_Spades_output/}
    # Get full path to tarball
    tarball_name="$(dirname ${src_asm_graph}).tar.gz"
    tarball_path="${tarballs_directory}/${tarball_name}"
    # Now finish destination path
    dest_scaffolds="$(basename ${src_asm_graph})"
    # print paths for housekeeping
    echo "sponge: ${sponge} -> tarball path: ${tarball_path} -> output path: ${outdir}/${sponge}.assembly_graph_with_scaffolds.gfa"
    # Good idea to check naming transforms with command below before running script.
    # tar -t -v --strip-components=1 --show-transformed-names  --file=${tarball_path} ${src_asm_graph}
    if [ ! -f ${outdir}/${sponge}.assembly_graph_with_scaffolds.gfa ];
    then
        tar --keep-old-files --strip-components=1 --transform="s/assembly_graph_with_scaffolds/${sponge}.assembly_graph_with_scaffolds/" --file=${tarball_path} --directory ${outdir} --get ${src_asm_graph} &
    else
        echo "${outdir}/${sponge}.assembly_graph_with_scaffolds.gfa already exists. skipping..."
    fi
done

wait