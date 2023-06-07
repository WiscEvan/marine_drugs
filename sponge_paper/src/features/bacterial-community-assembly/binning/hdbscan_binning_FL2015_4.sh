#!/bin/bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.hdscan_binning_FL2015_4.err
#SBATCH --output=logs/%J.hdscan_binning_FL2015_4.out



# Shared Data
norm_kmers="/home/evan/sponge_paper/sponge_paper/data/interim/binning/FL2015_4.kmers.am_clr.tsv"
coverages="/home/evan/sponge_paper/sponge_paper/data/interim/binning/FL2015_4.coverages.tsv"
embedded_kmers="/home/evan/sponge_paper/sponge_paper/data/interim/binning/FL2015_4.kmers.am_clr.bhsne.tsv"
taxonomy="/home/evan/sponge_paper/sponge_paper/data/interim/binning/FL2015_4.taxonomy.tsv"

# Archaea specific
archaea_markers="/home/evan/sponge_paper/sponge_paper/data/interim/binning/FL2015_4.archaea.markers.tsv"
archaea_out="/home/evan/sponge_paper/sponge_paper/data/interim/binning/FL2015_4.archaea.binning.hdbscan.tsv"
 
/home/evan/miniconda3/envs/autometa/bin/autometa-binning $norm_kmers $coverages $archaea_markers $archea_out \
    --embedded-kmers $embedded_kmers \
    --taxonomy $taxonomy \
    --clustering-method hdbscan \
    --domain archaea

# Bacteria specific
bacteria_out="/home/evan/sponge_paper/sponge_paper/data/interim/binning/FL2015_4.bacteria.binning.hdbscan.tsv"
bacteria_markers="/home/evan/sponge_paper/sponge_paper/data/interim/binning/FL2015_4.bacteria.markers.tsv"

/home/evan/miniconda3/envs/autometa/bin/autometa-binning $norm_kmers $coverages $bacteria_markers $bacteria_out \
    --embedded-kmers $embedded_kmers \
    --taxonomy $taxonomy \
    --clustering-method hdbscan \
    --domain bacteria