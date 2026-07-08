#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Pipeline for diagnostic rnaseq analyses
*/

// ----------------Workflow---------------- //

include { RNA_DIAGNOSTIC } from './workflows/local/rna.nf'
include { CONTROL_COHORT_PREP } from './workflows/local/prep_control_cohort.nf'

workflow {

  // Channel for source file: txt file with sample ID and file name for each fastq file
  // Use 2 lines with the same ID for paired-end samples
  Channel
    .fromPath(params.source_file, checkIfExists: true)
    .set{ source_file }

  // Run workflow
  if (params.rna_workflow) {
    
    RNA_DIAGNOSTIC(source_file)

  }
  
  if (params.control_cohort_prep) {
    
    CONTROL_COHORT_PREP(source_file)

  }

}