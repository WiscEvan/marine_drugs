#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.cpus = 1

process GENOME_GENERATE {
    tag "Making STAR database for ${assembly.simpleName}"
    cpus params.cpus

    input:
        path assembly

    output:
        path "${assembly.simpleName}", type: 'dir'
    
    """
    STAR \
        --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir ${assembly.simpleName} \
        --limitGenomeGenerateRAM 1811573025717 \
        --genomeFastaFiles ${assembly}
    """
}

process ALIGN_RNA {
    tag "STAR alignment on ${genomeDir.simpleName}"
    cpus params.cpus

    input:
        path genomeDir
        tuple val(readsName), path(reads)

    output:
        path "${genomeDir.simpleName}.Aligned.sortedByCoord.out.bam", emit: bam
        path "${genomeDir.simpleName}.Signal.Unique.str1.out.wig", emit: wig
    
    """
    # Note: We are not passing in single-end reads with the paired-end reads
    # as this was not being recognized by STAR
    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${genomeDir} \
        --outFileNamePrefix ${genomeDir.simpleName}. \
        --readFilesIn ${reads} \
        --readFilesCommand zcat \
        --alignIntronMax 100000 \
        --outSAMtype BAM SortedByCoordinate \
        --outWigType wiggle \
        --outWigStrand Unstranded
    """
}

workflow RNASEQ_ALIGNMENT {
    take:
        assembly
        reads

    main:
        GENOME_GENERATE(assembly)
        ALIGN_RNA(GENOME_GENERATE.out, reads)

    emit:
        bam = ALIGN_RNA.out.bam
        wig = ALIGN_RNA.out.wig
}