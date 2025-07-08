#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Pipeline for diagnostic rnaseq analyses
*/

// ----------------Workflow---------------- //

include { RNA_DIAGNOSTIC } from './workflows/local/rna.nf'

workflow {

  // Channel for source file (txt file with sample ID and file name for each fastq file (2 line for paired-end samples))
  Channel
    .fromPath(params.source_file, checkIfExists: true)
    .set{ source_file }

  // Run workflow
  RNA_DIAGNOSTIC(source_file)

}