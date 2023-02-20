#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EXTRACT; ALIGN; TRIM; COMBINE; ML_TREE; BI_TREE } from './tasks.nf'

workflow MTDNA {
  take:
    mtdna

  main:
    // Perform various annotations on provided metagenome
    // retrieve_mtdna.py
    // find_trna.sh (arwen)
    // find_trna.sh (genbank_format.py)
    // <<separately align mitogenome genes>>
    // <<concatenate sorted mitogenome genes>>
    // sort_trna.py (arwen outdir -> )
    // format_headers.py || get_seqs.py
    EXTRACT(mtdna)
    ALIGN(EXTRACT.out)
    TRIM(ALIGN.out)
    COMBINE( TRIM.out.collect() )
    // Perform taxon assignment with filtered metagenome
    ML_TREE(COMBINE.out)
    // Now perform binning with all of our annotations.
    BI_TREE(COMBINE.out)

  emit:
    binning = BINNING.out.binning
    recruitment = UNCLUSTERED_RECRUITMENT.out
    all_trees = ML_TREE.out | mix(BI_TREE.out) | collect
}
