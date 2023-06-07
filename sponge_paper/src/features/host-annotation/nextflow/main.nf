#!/usr/bin/env nextflow

nextflow.enable.dsl=2


params.assembly = "$HOME/sponge_paper/sponge_paper/data/interim/assemblies/FL2015_9.filtered.fna"
params.reads = "$HOME/sponge_paper/sponge_paper/data/interim/rnaseq/FL2015-9_RNA_seqyclean_PE{1,2}.fastq.gz"
params.outdir = "$HOME/sponge_paper/sponge_paper/data/processed/sponge-markers"
params.species = null
if ( !params.species || params.species instanceof Boolean ) 
    error """
    No species identifier provided! 
    Run `augustus --species=help` for list of available species.
    Then pass in with `nextflow run --species=<selected-identifier> ...`
    """
params.extrinsicCfgFile = "$HOME/sponge_paper/sponge_paper/src/features/host-annotation/nextflow/extrinsic.cfg"
params.cpus = 65

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
    Channel
        .fromPath( params.extrinsicCfgFile, checkIfExists: true, type: 'file')
        .set{ extrinsic_config_ch }

    RNASEQ_ALIGNMENT(assembly_ch, reads_ch)
    ANNOTATE_HOST(assembly_ch, RNASEQ_ALIGNMENT.out.bam, RNASEQ_ALIGNMENT.out.wig, extrinsic_config_ch)

}