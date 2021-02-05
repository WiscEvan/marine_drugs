#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.assembly = 'metagenome_assembly'
// See 
params.species = 'amphimedon'
params.reads = 'reads'
params.extrinsicCfgFile = 'extrinsic.cfg'
params.cpus = 1
params.threads = 1

process predict_genes_abinitio {
    tag "predicting genes (ab initio) for ${assembly.simpleName}"
    container 'augustus:latest'
    containerOptions "--volume $INDIR:/input:ro"
    publishDir "host-annotation", mode: 'copy'
    
    input:
      path x
    
    output:
      path "*.abinitio.*"
    
    script:
      """
      augustus \
        --species=$params.species \
        --sample=${assembly.simpleName} \
        --genemodel=partial \
        --protein=on \
        --codingseq=on \
        $assembly > ${assembly.simpleName}.abinitio.gff
      getAnnoFasta.pl ${assembly.simpleName}.abinitio.gff > ${assembly.simpleName}.abinitio.faa
      """
}

process ALIGN_READS {
    tag "predicting genes (ab initio) for ${assembly.simpleName}"
    container 'augustus:latest'
    containerOptions "--volume $INDIR:/input:ro"
    publishDir "host-annotation", mode: 'copy'
    
    input:
      path x
    
    output:
      path "*.abinitio.*"
    
    script:
    """  
    mkdir gindex
    STAR \
    --runThreadN $params.threads \
    --runMode genomeGenerate \
    --genomeFastaFiles $assembly
    
    # Note: We are not passing in single-end reads with the paired-end reads
    # as this was not being recognized by STAR
    
    STAR \
        --runThreadN $threads \
        --outFileNamePrefix ${assembly.simpleName} \
        --readFilesIn $fwd_reads $rev_reads \
        --readFilesCommand zcat \
        --alignIntronMax 100000 \
        --outSAMtype BAM SortedByCoordinate \
        --outWigType wiggle \
        --outWigStrand Unstranded
    """
}
process retrieve_hints {
    tag "predicting genes (ab initio) for ${assembly.simpleName}"
    container 'augustus:latest'
    containerOptions "--volume $INDIR:/input:ro"
    publishDir "host-annotation", mode: 'copy'
    
    input:
    path bam
    path wig
    // Regex for Augustus output:
    // bam: "*_Aligned.sortedByCoord.out.bam"
    // wig: *_Signal.Unique.str1.out.wig
    
    output:
    path 'hints.*.gff'
    
    script:
    """
    # Get intron hints
    bin/bam2hints --intronsonly --in=$bam --out=hints.introns.gff
    # Get exon hints
    # Input: *_Signal.Unique.str1.out.wig
    cat $wig \
        | \
        ./scripts/wig2hints.pl \
            --width=10 \
            --margin=10 \
            --minthresh=2 \
            --minscore=4 \
            --prune=0.1 \
            --src=W \
            --type=ep \
            --radius=4.5 \
            --pri=4 \
            --strand="." > hints.exons.gff
    cat "${SPONGE_ALN_DIR}/${intron_hints_filename}" $ep_hints > $HINTS
      """
}

process predict_genes_with_hints {
    tag "predicting genes (ab initio) for ${assembly.simpleName}"
    container 'augustus:latest'
    containerOptions "--volume $INDIR:/input:ro"
    publishDir "host-annotation", mode: 'copy'
    
    input:
      path x
    
    output:
      path 
    
    script:
      """  
      augustus --species=$SPECIES --sample=${sponge} --genemodel=partial /input/eukaryota.fna > $abinitio_output
      """
}

process generate_gene_tracks {
    tag "predicting genes (ab initio) for ${assembly.simpleName}"
    container 'augustus:latest'
    containerOptions "--volume $INDIR:/input:ro"
    publishDir "host-annotation", mode: 'copy'
    
    input:
      path x
    
    output:
      path 
    
    script:
      """  
      augustus --species=$SPECIES --sample=${sponge} --genemodel=partial /input/eukaryota.fna > $abinitio_output
      """
}

process BUSCO_GENOME {
    tag "predicting genes (ab initio) for ${assembly.simpleName}"
    container 'augustus:latest'
    containerOptions "--volume $INDIR:/input:ro"
    publishDir "host-annotation", mode: 'copy'
    
    input:
      path assembly
    
    output:
      path "${assembly.simpleName}_busco_genome"
    
    script:
      """  
      busco -m genome --auto-lineage-euk -i $assembly -o ${assembly.simpleName}_busco_genome
      """
}

process BUSCO_PROTEIN {
    tag "predicting genes (ab initio) for ${assembly.simpleName}"
    container 'augustus:latest'
    containerOptions "--volume $INDIR:/input:ro"
    publishDir "host-annotation", mode: 'copy'
    
    input:
      path assembly
    
    output:
      path "${assembly.simpleName}_busco_genome"
    
    script:
      """  
      busco \
        --mode protein \
        --auto-lineage-euk \
        --cpu $params.cpus \
        
        --in $proteins \
        --out ${proteins.simpleName}_busco_protein
      """
}
workflow ANNOTATE_HOST {
    Channel
        .fromPath( params.assembly, checkIfExists: true, type: 'file')
        .set{ assembly_ch }
    Channel
        .fromPath( params.reads, checkIfExists: true, type: 'file')
        .set{ reads_ch }

    ALIGN_READS( params.assembly; reads_ch)   
    predict_genes_abinitio( assembly_ch )
    retrieve_hints( reads_ch )
    predict_genes_with_hints( )
    generate_gene_tracks( )
}