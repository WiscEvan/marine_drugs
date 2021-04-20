#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reads = "./data/raw/*.fastq"
params.database = "./data/external/"
params.interim = "./data/interim"
params.processed = "./data/processed"
params.cpus = 2

log.info """
 =======================================================================
 reads              : ${params.reads}
 database           : ${params.database}
 interim            : ${params.interim}
 processed          : ${params.processed}
 cpus               : ${params.cpus}
 -----------------------------------------------------------------------
"""


process CONCAT_PAIRS {
    tag "concatenating ${readsName}"
    publishDir params.interim, pattern: '*.fq.gz', mode:'symlink'

    input:
      tuple val(readsName), path(reads)
    
    output:
      path "${readsName}.fq.gz"

    script:
    """
    zcat $reads > ${readsName}.fq
    gzip ${readsName}.fq
    """
}

// Profile community presence/absence and abundance of microbial pathways with humann
// HUMAnN is a pipeline for efficiently and accurately profiling the presence/absence
// and abundance of microbial pathways in a community
process HUMANN {
  tag "functional profiling ${reads.simpleName}"
  container = 'biobakery/humann:latest'
  cpus params.cpus
  memory { reads.size() < 10.GB ? 64.GB : 124.GB }
  publishDir params.interim, pattern: '*.tsv', mode:'symlink'
  
  input:
    path reads
    path database
  
  output:
    path "*"
    path "*_genefamilies.tsv", emit: gene_families
    path "*_pathabundance.tsv", emit: path_abundance
    path "*_pathcoverage.tsv", emit: path_coverage

  script:
  """
  humann_config --update database_folders nucleotide ${database}/chocophlan
  humann_config --update database_folders protein ${database}/uniref
  humann_config --update database_folders utility_mapping ${database}/utility_mapping
  # Now run HUMANN with the configured databases from S3
  # Providing --bowtie-options '--very-sensitive-local' is equivalent to 'soft-clipping'
  humann \
    --threads ${task.cpus} \
    --input $reads \
    --output . \
    --bowtie-options='--very-sensitive-local' \
    --metaphlan-options="--nproc ${task.cpus}" \
    --remove-temp-output \
    --resume
  """
}

process RENORM_PATHABUNDANCES {
  tag "renormalizing ${pathabundances.simpleName} to copies per million"
  container = 'biobakery/humann:latest'
  publishDir params.interim, pattern: "*.tsv", mode: 'symlink'

  input:
    path pathabundances

  output:
    path "${pathabundances.simpleName}_cpm.tsv"

  """
  humann_renorm_table  \
    --input $pathabundances \
    --units cpm \
    --mode community \
    --output ${pathabundances.simpleName}_cpm.tsv
  """
}

process JOIN_PATHABUNDANCE_PROFILES {
  tag 'joining community pathabundance profiles'
  container = 'biobakery/humann:latest'
  publishDir params.processed, pattern: 'master_community_pathabundance_profiles.tsv', mode:'copy'

  input:
    path pathabundance_profiles

  output:
    path "master_community_pathabundance_profiles.tsv"

  """
  mkdir ./data
  mv $pathabundance_profiles ./data/.
  humann_join_tables \
    --input ./data \
    --file_name "pathabundance" \
    --output master_community_pathabundance_profiles.tsv
  """
}

process JOIN_PATHCOVERAGE_PROFILES {
  tag 'joining community pathcoverage profiles'
  container = 'biobakery/humann:latest'
  publishDir params.processed, pattern: 'master_community_pathcoverage_profiles.tsv', mode:'copy'

  input:
    path pathcoverage_profiles

  output:
    path "master_community_pathcoverage_profiles.tsv"

  """
  mkdir ./data
  mv $pathcoverage_profiles ./data/.
  humann_join_tables \
    --input ./data \
    --file_name "pathcoverage" \
    --output master_community_pathcoverage_profiles.tsv
  """
}

process RENORM_GENEFAMILIES {
  tag "renormalizing ${genefamilies.simpleName} to copies per million"
  container = 'biobakery/humann:latest'
  publishDir params.interim, pattern: "*.tsv", mode: 'symlink'

  input:
    path genefamilies

  output:
    path "${genefamilies.simpleName}_cpm.tsv"

  """
  humann_renorm_table  \
    --input $genefamilies \
    --units cpm \
    --mode community \
    --output ${genefamilies.simpleName}_cpm.tsv
  """
}

process JOIN_GENEFAMILIES_PROFILES {
  tag 'joining community genefamilies profiles'
  container = 'biobakery/humann:latest'
  publishDir params.processed, pattern: 'master_community_genefamilies_profiles.tsv', mode:'copy'

  input:
    path genefamilies_profiles

  output:
    path "master_community_genefamilies_profiles.tsv"

  """
  mkdir ./data
  mv $genefamilies_profiles ./data/.
  humann_join_tables \
    --input ./data \
    --file_name "genefamilies" \
    --output master_community_genefamilies_profiles.tsv
  """
}

workflow {
  Channel
    .fromPath( params.reads, checkIfExists: true, type: 'file')
    .set{ reads_ch }

  Channel
    .value(params.database)
    .set{database}
  
  //CONCAT_PAIRS( reads_ch )
  HUMANN( reads_ch, database )
  RENORM_PATHABUNDANCES( HUMANN.out.path_abundance )
  JOIN_PATHABUNDANCE_PROFILES( RENORM_PATHABUNDANCES.out.collect() )
  JOIN_PATHCOVERAGE_PROFILES( HUMANN.out.path_coverage.collect() )
  RENORM_GENEFAMILIES( HUMANN.out.gene_families )
  JOIN_GENEFAMILIES_PROFILES( RENORM_GENEFAMILIES.out.collect() )
}

/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}
