#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process EXTRACT {
  tag "Extracting ${fasta.simpleName}"

  input:
    path fasta

  output:
    path "${fasta.simpleName}.fna"

  """
  TODO: 
  """
}

process ALIGN {
  tag "filtering metagenome ${metagenome.simpleName}"
  container = 'jason-c-kwan/autometa:dev'
  publishDir params.interim, pattern: "${metagenome.simpleName}.filtered.fna"

  input:
    path gene

  output:
    path "${metagenome.simpleName}.filtered.fna"

  """
  TODO: 
  """
}

process TRIM {
  tag "combining mtDNA genes nexus files"

  input:
    path genes

  output:
    path "mtdna.combined.nexus"

  """
  TODO: 
  """
}

process MAKE_NEXI {
  tag "combining mtDNA genes nexus files"

  input:
    path trimmed_alignments

  output:
    path "mtdna.combined.nexus"

  """
  TODO: 
  """
}

process COMBINE_NEXI {
  tag "combining mtDNA genes nexus files"

  input:
    path nexi

  output:
    path "mtdna.combined.nexus"

  """
  concat_nexus.py --input $nexi --output mtdna.combined.nexus
  """
}

process ML_TREE {
  tag "combining mtDNA genes nexus files"

  input:
    path nexi

  output:
    path "maximum_likelihood.tre"

  """
  TODO: 
  """
}

process BI_TREE {
  tag "combining mtDNA genes nexus files"

  input:
    path nexi

  output:
    path "bayesian_inference.tre"

  """
  TODO: 
  """
}