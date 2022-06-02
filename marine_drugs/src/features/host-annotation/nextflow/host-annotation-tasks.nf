#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.assembly = 'metagenome_assembly'
params.reads = 'reads'
params.outdir = 'output_directory'
params.species = 'amphimedon'
params.extrinsicCfgFile = 'extrinsic.cfg'
params.cpus = 1

process INTRON_HINTS {
    tag "Retrieving intron hints for ${bam.simpleName}"
    container 'augustus:latest'
    // publishDir params.outdir, pattern: "*.hints.introns.gff", mode: 'copy'
    
    input:
      path bam
    
    output:
      path "${bam.simpleName}.hints.introns.gff"
    
    script:
      """
      # Get intron hints
      bam2hints --intronsonly --in=${bam} --out=${bam.simpleName}.hints.introns.gff
      """
}

process EXONPARTS_HINTS {
    tag "Retrieving exonparts hints for ${wig.simpleName}"
    container 'augustus:latest'
    // publishDir params.outdir, pattern: "*.hints.exonparts.gff", mode: 'copy'
    
    input:
      path wig
    
    output:
      path "${wig.simpleName}.hints.exonparts.gff"
    
    script:
      """
      cat ${wig} \
          | wig2hints.pl \
              --width=10 \
              --margin=10 \
              --minthresh=2 \
              --minscore=4 \
              --prune=0.1 \
              --src=W \
              --type=ep \
              --radius=4.5 \
              --pri=4 \
              --strand="." > ${wig.simpleName}.hints.exonparts.gff
      """
}

process RETRIEVE_HINTS {
    tag "Concatenating hints ${introns.simpleName} and ${exonparts.simpleName}"
    
    publishDir params.outdir, pattern: "hints.gff", mode: 'copy'

    input:
      path introns
      path exonparts
    
    output:
      path "${introns.simpleName}.hints.gff"
    
    script:
      """
      cat $introns $exonparts > ${introns.simpleName}.hints.gff
      """
}

process PREDICT_GENES {
    tag "Augustus gene prediction w/hints ${assembly.simpleName}"
    container 'augustus:latest'
    // containerOptions "--volume ${params.containerVolume}:/input:ro"
    publishDir params.outdir, pattern: "*.predictions.gff", mode: 'copy'
    
    input:
      path assembly
      path hints
      path extrinsicCfgFile
    
    output:
      path "${assembly.simpleName}.predictions.gff"
    
    script:
      """  
      augustus \
        --species=${params.species} \
        --sample=${assembly.simpleName} \
        --alternatives-from-evidence=true \
        --extrinsicCfgFile=${extrinsicCfgFile} \
        --allow_hinted_splicesites=atac \
        --softmasking=on \
        --protein=on \
        --codingseq=on \
        --UTR=on \
        --hintsfile=${hints} \
        ${assembly} > ${assembly.simpleName}.predictions.gff
      """
}

process EXTRACT_ORFS {
    tag "Extracting ORFS for ${predictions.simpleName}"
    container 'augustus:latest'
    publishDir params.outdir, pattern: "*.{codingseq,aa,cdsexons,mrna}", mode: 'copy'
    
    input:
      path assembly
      path predictions
    
    output:
      path "*.aa", emit: orfs
      path "*.codingseq", emit: codingseq
      path "*.cdsexons", emit: cdsexons
      path "*.mrna", emit: mrna
    
    script:
      """
      # script will generate *.codingseq, *.aa, *.cdsexons, *.mrna
      getAnnoFasta.pl --seqfile=${assembly} ${predictions}
      """
}

process ANNOTATE_MARKERS {
    tag "BUSCO marker prediction on ${proteins.simpleName}"
    cpus params.cpus
    container 'ezlabgva/busco:v5.beta.1_cv1'
    containerOptions '-u $(id -u)'
    publishDir params.outdir, pattern: "*_markers", mode: 'copy'
    
    input:
      path proteins
    
    output:
      path "${proteins.simpleName}_markers"
    
    script:
      """  
      busco \
        --mode proteins \
        --lineage_dataset metazoa_odb10 \
        --cpu ${task.cpus} \
        --in ${proteins} \
        --out ${proteins.simpleName}_markers
      """
}

workflow ANNOTATE_HOST {
  take:
    assembly
    bam
    wig
    extrinsicCfgFile

  main:
    INTRON_HINTS(bam)
    EXONPARTS_HINTS(wig)
    RETRIEVE_HINTS(INTRON_HINTS.out, EXONPARTS_HINTS.out)
    PREDICT_GENES(assembly, RETRIEVE_HINTS.out, extrinsicCfgFile)
    EXTRACT_ORFS(assembly, PREDICT_GENES.out)
    ANNOTATE_MARKERS(EXTRACT_ORFS.out.orfs)

  emit:
    predictions = PREDICT_GENES.out
    orfs = EXTRACT_ORFS.out.orfs
    mrna = EXTRACT_ORFS.out.mrna
    codingseq = EXTRACT_ORFS.out.codingseq
    cdsexons = EXTRACT_ORFS.out.cdsexons
    markers = ANNOTATE_MARKERS.out
}