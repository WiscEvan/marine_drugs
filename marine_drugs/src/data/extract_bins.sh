#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.extract_bins.err
#SBATCH --output=logs/%J.extract_bins.out


fastas="$HOME/marine_drugs/marine_drugs/data/interim/assemblies"
binning_dir="$HOME/marine_drugs/marine_drugs/data/interim/binning"
outdir="$HOME/marine_drugs/marine_drugs/data/processed/bins"
src="/home/evan/marine_drugs/marine_drugs/src/data"
binning_filepaths="${binning_dir}/final_binning_filepaths.tsv"

for fasta in `ls ${fastas}/*.filtered.fna`;do
    # Example fasta name: FL2015_30.filtered.fna
    sponge=$(basename ${fasta/.filtered.fna/})
    binning=`grep "${sponge}\." $binning_filepaths | cut -f2`
    bin_outdir="${outdir}/${sponge}"
    if [ -f $binning ]
    then python "${src}/extract_bins.py" --use-latest-refinement --binning $binning --fasta $fasta --outdir $bin_outdir
    else echo "${binning} does not exists. Skipping..."
    fi
done