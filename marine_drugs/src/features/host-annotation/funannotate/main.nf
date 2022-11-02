#!/usr/bin/env nextflow

nextflow.enable.dsl=2


params.sample_species = null
if ( !params.sample_species || params.sample_species instanceof Boolean ) 
    error """
    No sample species identifier provided! What sponge sample is this?
    """
params.extrinsicCfgFile = "$HOME/marine_drugs/marine_drugs/src/features/host-annotation/nextflow/extrinsic.cfg"
params.cpus = 65
params.busco_db = 'metazoa'

include { FUNANNOTATE } from './funannotate-tasks.nf'

log.info """
 Sponge-annotation - Automated host-functional annotation
 =====================================================
 projectDir           : ${workflow.projectDir}
 -----------------------------------------------------
 Data & Input Parameters
 =====================================================
 assembly             : ${params.assembly}
 fwd_reads            : ${params.fwd_reads}
 rev_reads            : ${params.rev_reads}
 se_reads             : ${params.se_reads}
 sample_species       : ${params.sample_species}
 outdir               : ${params.outdir}
 -----------------------------------------------------
 Runtime Parameters
 =====================================================
 cpus                 : ${params.cpus}
 memory               : ${params.memory}
 -----------------------------------------------------
 Databases
 =====================================================
 extrinsicCfgFile     : ${params.extrinsicCfgFile}
 augustus_species     : ${params.augustus_species}
 busco_db             : ${params.busco_db}
 -----------------------------------------------------
"""

workflow {
    Channel
        .fromPath( params.assembly, checkIfExists: true, type: 'file')
        .set{ assembly_ch }
    Channel
        .fromPath( params.fwd_reads, checkIfExists: true, type: 'file')
        .collect()
        .set{ fwd_reads_ch }
    Channel
        .fromPath( params.rev_reads, checkIfExists: true, type: 'file')
        .collect()
        .set{ rev_reads_ch }
    Channel
        .fromPath( params.se_reads, checkIfExists: true, type: 'file')
        .collect()
        .set{ se_reads_ch }

    FUNANNOTATE(assembly_ch, fwd_reads_ch, rev_reads_ch, se_reads_ch)

}