.PHONY:
	help env

help : 
	@echo "Please inspect Makefile for list of commands"

env:
	@echo "Constructing conda environment: name=sponges"
	conda create -n sponges -c bioconda -c plotly -c conda-forge python=3.8 --file=requirements.txt --yes

# The wildcard matching needs to be fixed... But generally this is workflow for reproducibility
marine_drugs/data/assemblies/assemblies_stats.tsv: marine_drugs/data/raw/assemblies/*/*.fasta
	# @metagenomes=`ls ../data/assemblies/*/*.fasta`
	python marine_drugs/src/features/stat_metagenome.py $^ $@

marine_drugs/data/assemblies/%/scaffolds.fasta: marine_drugs/src/data/extract_scaffolds.sh marine_drugs/data/raw/assemblies/%_tarball_listing.txt
	bash marine_drugs/src/data/extract_scaffolds.sh

marine_drugs/src/data/list_assemblies_tarballs.sh: /media/external2/Sponge_assemblies/*.tar.gz
	for f in `ls /media/external2/Sponge_assemblies/*.tar.gz`;do echo "tar -tf $f > marine_drugs/data/raw/assemblies/$(basename ${f/.tar.gz/_tarball_listing.txt})";done > $@