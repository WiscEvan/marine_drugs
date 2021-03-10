#!/usr/bin/env nextflow

nextflow.enable.dsl=2


params.assembly = '/home/evan/marine_drugs/marine_drugs/data/interim/assemblies/FL2015_5.filtered.fna'
params.reads = '/media/bigdrive2/evan/marine_drugs/marine_drugs/data/interim/rnaseq/FL2015-5_RNA_seqyclean_PE{1,2}.fastq.gz'
params.extrinsicCfgFile = '/media/bigdrive2/evan/marine_drugs/marine_drugs/src/features/host-annotation/augustus-gene-calling/extrinsic.cfg'
params.cpus = 50
params.species = 'amphimedon'

include { RNASEQ_ALIGNMENT } from './rnaseq-tasks.nf'
include { ANNOTATE_HOST } from './host-annotation-tasks.nf'

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