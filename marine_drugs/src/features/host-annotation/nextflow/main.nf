#!/usr/bin/env nextflow

nextflow.enable.dsl=2


params.assembly = '/home/evan/marine_drugs/marine_drugs/data/interim/assemblies/FL2015_5.filtered.fna'
params.reads = '/home/evan/marine_drugs/marine_drugs/data/interim/rnaseq/FL2015-5_RNA_seqyclean_PE{1,2}.fastq.gz'
params.outdir = '/home/evan/marine_drugs/marine_drugs/data/processed/sponge-markers'
params.species = 'amphimedon'
params.extrinsicCfgFile = '/home/evan/marine_drugs/marine_drugs/src/features/host-annotation/augustus-gene-calling/extrinsic.cfg'
params.cpus = 50

include { RNASEQ_ALIGNMENT } from './rnaseq-tasks.nf'
include { ANNOTATE_HOST } from './host-annotation-tasks.nf'

log.info """
 Sponge-annotation - Automated host-marker analysis
 =====================================================
 projectDir           : ${workflow.projectDir}
 -----------------------------------------------------
 Data
 =====================================================
 assembly             : ${params.assembly}
 reads                : ${params.reads}
 outdir               : ${params.outdir}
 -----------------------------------------------------
 Parameters
 =====================================================
 cpus                 : ${params.cpus}
 -----------------------------------------------------
 Databases
 =====================================================
 extrinsicCfgFile     : ${params.extrinsicCfgFile}
 species              : ${params.species}
 -----------------------------------------------------
"""

workflow {
    Channel
        .fromPath( params.assembly, checkIfExists: true, type: 'file')
        .set{ assembly_ch }
    Channel
        .fromFilePairs( params.reads, checkIfExists: true, type: 'file')
        .set{ reads_ch }

    RNASEQ_ALIGNMENT(assembly_ch, reads_ch)
    ANNOTATE_HOST(assembly_ch, RNASEQ_ALIGNMENT.out.bam, RNASEQ_ALIGNMENT.out.wig)

}