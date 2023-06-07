#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.assembly = 'metagenome_assembly'
params.fwd_reads = 'fwd_reads'
params.rev_reads = 'rev_reads'
params.se_reads = 'se_reads'
params.outdir = 'output_directory'
params.sample_species = 'sample species name'
params.busco_db = 'metazoa'
params.augustus_species = 'amphimedon'
params.extrinsicCfgFile = 'extrinsic.cfg'
params.cpus = 1
params.memory = '50GB'

process CLEAN {
    tag "cleaning ${assembly.simpleName}"
    container 'nextgenusfs/funannotate:nextflow'
    
    input:
      path assembly
    
    output:
      path "${assembly.simpleName}.cleaned.fa"
    
    script:
      """
      funannotate clean -i $assembly --minlen 1000 -o ${assembly.simpleName}.cleaned.fa
      """
}

process SORT {
    tag "sorting scaffolds by length ${assembly.simpleName}"
    container 'nextgenusfs/funannotate:nextflow'
    
    input:
      path assembly
    
    output:
      path "${assembly.simpleName}.sorted.fa"
    
    script:
      """
      funannotate sort -i $assembly -b scaffold -o ${assembly.simpleName}.sorted.fa
      """
}

process MASK {
    tag "Soft-masking repetitive elements ${assembly.simpleName}"
    container 'nextgenusfs/funannotate:nextflow'
    cpus params.cpus
    
    input:
      path assembly
    
    output:
      path "${assembly.simpleName}.masked.fa"
    
    script:
      """
      funannotate mask -i $assembly --cpus ${task.cpus} -o ${assembly.simpleName}.masked.fa
      """
}

process TRAIN {
    tag "training RNAseq against ${assembly.simpleName}"
    container 'nextgenusfs/funannotate:nextflow'
    cpus params.cpus
    
    input:
      path assembly
      tuple val(fwdReadsName), path(fwd_reads)
      tuple val(revReadsName), path(rev_reads)
      tuple val(seReadsName), path(se_reads)
    
    output:
      path "${assembly.simpleName}"
    
    script:
      """
      funannotate train \
        --input $assembly \
        --out ${assembly.simpleName} \
        --left $fwd_reads \
        --right $rev_reads \
        --single $se_reads \
        --stranded RF \
        --jaccard_clip \
        --species "${params.sample_species}" \
        --cpus ${task.cpus}
      """
}

process PREDICT {
    tag "FUNANNOTATE prediction on ${assembly.simpleName}"
    container 'nextgenusfs/funannotate:nextflow'
    cpus params.cpus
    
    input:
      path assembly
      path trainDir
    
    output:
      path "${assembly.simpleName}"
    
    script:
      """
      funannotate predict \
        --input $assembly \
        --out $trainDir \
        --species "${params.sample_species}" \
        --busco_db ${params.busco_db} \
        --augustus_species ${params.augustus_species} \
        --busco_seed_species ${params.augustus_species} \
        --organism other \
        --cpus ${task.cpus}
      """
}

process UPDATE {
    tag "FUNANNOTATE prediction on ${predictDir.simpleName}"
    container 'nextgenusfs/funannotate:nextflow'
    cpus params.cpus
    memory params.memory
    
    input:
      path predictDir
    
    output:
      path "${predictDir}"
    
    script:
      """
      funannotate update \
        --input $predictDir \
        --cpus ${task.cpus} \
        --memory ${task.memory}
      """
}

process IPRSCAN {
    tag "InterProScan5 on ${predictDir.simpleName}"
    container 'nextgenusfs/funannotate:nextflow'
    cpus params.cpus
    
    input:
      path predictDir
    
    output:
      path "${predictDir.simpleName}.iprscan.xml"
    
    script:
      """
      funannotate iprscan \
        --input $predictDir \
        --method local \
        --cpus ${task.cpus} \
        --out ${predictDir.simpleName}.iprscan.xml
      """
}

process ANNOTATE {
    tag "functional annotation on ${predictDir.simpleName}"
    container 'nextgenusfs/funannotate:nextflow'
    cpus params.cpus
    memory params.memory
    publishDir params.outdir, mode: 'move'
    
    input:
      path predictDir
      path iprscan
    
    output:
      path "${predictDir.simpleName}_annotations"
    
    script:
      """
      funannotate annotate \
        --input $predictDir \
        --cpus ${task.cpus} \
        --busco_db ${params.busco_db} \
        --out ${predictDir.simpleName}_annotations \
        --iprscan $iprscan
      """
}

workflow FUNANNOTATE {
    take:
        assembly
        fwd_reads
        rev_reads
        se_reads
    main:
        CLEAN(assembly)
        SORT(CLEAN.out)
        MASK(SORT.out)
        TRAIN(MASK.out, fwd_reads, rev_reads, se_reads)
        PREDICT(MASK.out, TRAIN.out)
        UPDATE(PREDICT.out)
        IPRSCAN(UPDATE.out)
        ANNOTATE(UPDATE.out, IPRSCAN.out)
}