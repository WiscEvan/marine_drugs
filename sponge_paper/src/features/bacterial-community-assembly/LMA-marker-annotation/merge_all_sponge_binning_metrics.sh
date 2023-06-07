#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=14-00:00:00
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=logs/%J.merge_all_sponge_binning_metrics.err
#SBATCH --output=logs/%J.merge_all_sponge_binning_metrics.out



HMMSCAN_FILES="${HOME}/sponge_paper/sponge_paper/data/interim/LMA-marker-annotation/*.hmmscan.tsv"
BIN_FILES="${HOME}/sponge_paper/sponge_paper/data/interim/binning/final_binning_filepaths.tsv"
SRC="${HOME}/sponge_paper/sponge_paper/src/features/bacterial-community-assembly/LMA-marker-annotation/merge_all_sponge_binning_metrics.py"

INDIR="/home/evan/sponge_paper/sponge_paper/data/interim/LMA-marker-annotation"
OUTDIR="/home/evan/sponge_paper/sponge_paper/data/processed/LMA-marker-annotation"
if [ ! -d $OUTDIR ]
then mkdir -p $OUTDIR
fi

outfiles=()
for rank in phylum class;do
    out="${OUTDIR}/LMA_sponges_${rank}_markers_cluster_metrics.tsv"
    python $SRC --indir $INDIR --rank $rank --out $out
    outfile+=($out)
done

main_outfile="${OUTDIR}/lma_sponges_taxon_specific_markers_cluster_metrics.tsv"
python -c """
import pandas as pd;
dfs = []
for fp in ${outfiles[@]}:
    df = pd.read_csv(${out}, sep='\t')
    dfs.append(df)
df = pd.concat(dfs)
df.to_csv($main_outfile, sep='\t', index=False, header=True)
"""

