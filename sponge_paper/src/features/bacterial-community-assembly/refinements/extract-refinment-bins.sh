#!/usr/bin/env bash

# NOTE: First you will need to clone Automappa.
# git clone https://github.com/WiscEvan/Automappa.git
extract="$HOME/Automappa/scripts/extract_refinement_bins.py"

# Directory of refinements output by Automappa
refinements=`ls $HOME/sponge_paper/sponge_paper/data/interim/refined/second_round_of_refinements/*.csv`
# all_metagenomes=`ls ~/sponge_paper/sponge_paper/data/raw/assemblies/*.fasta`
METAGENOMES="$HOME/sponge_paper/sponge_paper/data/raw/assemblies"
OUTDIR="$HOME/sponge_paper/sponge_paper/data/interim/refined/second_round_of_refinements"

for refinement in $refinements;do
    # Now update arguments
    filename=$(basename ${refinement})
    sponge=${filename/_refinements.csv/}
    metagenome="${METAGENOMES}/${sponge}.fasta"
    # Now we extract our refined bins to their respective output directory which we'll title by sponge.
    output="${OUTDIR}/${sponge}"
    bin_seq_counts="${output}/bin_sequence_counts.tsv"
    # echo "$sponge metagenome: ${metagenome}"
    python $extract --input $refinement --fasta $metagenome --output $output > $bin_seq_counts
    echo "Wrote counts to: ${bin_seq_counts}"
done