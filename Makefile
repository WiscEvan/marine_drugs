.PHONY:
	help env

help : 
	@echo "Please inspect Makefile for list of commands"

env:
	@echo "Constructing conda environment: name=sponges"
	conda create -n sponges -c bioconda -c plotly -c conda-forge python=3.8 --file=requirements.txt --yes

# The wildcard matching needs to be fixed... But generally this is workflow for reproducibility
sponge_paper/data/assemblies/assemblies_stats.tsv: sponge_paper/data/raw/assemblies/*/*.fasta
	# @metagenomes=`ls ../data/assemblies/*/*.fasta`
	python sponge_paper/src/features/stat_metagenome.py $^ $@

sponge_paper/data/assemblies/%/scaffolds.fasta: sponge_paper/src/data/extract_scaffolds.sh sponge_paper/data/raw/assemblies/%_tarball_listing.txt
	bash sponge_paper/src/data/extract_scaffolds.sh

sponge_paper/src/data/list_assemblies_tarballs.sh: /media/external2/Sponge_assemblies/*.tar.gz
	for f in `ls /media/external2/Sponge_assemblies/*.tar.gz`;do echo "tar -tf $f > sponge_paper/data/raw/assemblies/$(basename ${f/.tar.gz/_tarball_listing.txt})";done > $@