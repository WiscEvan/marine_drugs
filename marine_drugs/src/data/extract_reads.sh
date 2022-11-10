#!/usr/bin/env bash

#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=30
#SBATCH --mem 10000
#SBATCH --error=%J.reads_extraction.err
#SBATCH --output=%J.reads_extraction.out

tarballs_directory="/media/external2/Sponge_assemblies"
outdir="$HOME/marine_drugs/marine_drugs/data/raw/assemblies"
for reads in `grep --no-filename "_Spades_output/corrected/.*.*P.*.fastq.gz" ${outdir}/*_spades_tarball_listing.txt`;do
    # Example directory name: FL2015_30_Spades_output/corrected/FL2015_30_1P.00.0_0.cor.fastq.gz
    reads_filename=$(basename ${reads})
    sponge_fwd_rm=${reads_filename/_1P.00.0_0.cor.fastq.gz/}
    sponge=${sponge_fwd_rm/_2P.00.0_0.cor.fastq.gz/}
    # Get full path to tarball
    tarball_name="$(dirname $(dirname ${reads})).tar.gz"
    tarball_path="${tarballs_directory}/${tarball_name}"
    # Now finish destination path
    dest_scaffolds="$(basename ${reads})"
    # print paths for housekeeping
    echo "sponge: ${sponge} -> tarball path: ${tarball_path} -> output path: ${outdir}/${reads_filename}"
    # Good idea to check naming transforms with command below before running script.
    # tar -t -v --strip-components=2 --show-transformed-names --file=${tarball_path} ${reads}
    if [ ! -f ${outdir}/${reads_filename} ];
    then
        tar --keep-old-files --strip-components=2 --file=${tarball_path} --directory ${outdir} --get ${reads} &
    else
        echo "${outdir}/${reads_filename} already exists. skipping..."
    fi
done

wait