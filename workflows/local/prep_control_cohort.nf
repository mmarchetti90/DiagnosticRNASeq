/*
Workflow for diagnostic rnaseq analyses
*/

// ----------------Workflow---------------- //

// Subworkflows
include { LOAD_RESOURCES } from '../../subworkflows/local/load_resources.nf'
include { PREPROCESSING } from '../../subworkflows/local/preprocessing.nf'
include { ABERRANT_SPLICING } from '../../subworkflows/local/aberrant_splicing.nf'

workflow CONTROL_COHORT_PREP {

  take:
  source_file

  main:
  // LOADING RESOURCES -------------------- //

  LOAD_RESOURCES(source_file)

  LOAD_RESOURCES.out.scripts_dir.set{ scripts_dir }
  LOAD_RESOURCES.out.genome_fasta.set{ genome_fasta }
  LOAD_RESOURCES.out.genome_fasta_index.set{ genome_fasta_index }
  LOAD_RESOURCES.out.genome_annotation.set{ genome_annotation }
  LOAD_RESOURCES.out.star_index.set{ star_index }
  LOAD_RESOURCES.out.control_junc_dir.set{ control_junc_dir }
  LOAD_RESOURCES.out.control_gene_counts_dir.set{ control_gene_counts_dir }
  LOAD_RESOURCES.out.raw_reads.set{ raw_reads }
  LOAD_RESOURCES.out.sample_ids.set{ sample_ids }

  // PREPROCESSING ------------------------ //

  PREPROCESSING(raw_reads, scripts_dir, genome_fasta, genome_annotation, star_index, control_gene_counts_dir)

  PREPROCESSING.out.indexed_bam.set{ indexed_bam }
  PREPROCESSING.out.merged_counts.set{ merged_counts }
  PREPROCESSING.out.control_gene_counts_ids.set{ control_gene_counts_ids }

  // ABERRANT SPLICING -------------------- //

  ABERRANT_SPLICING(scripts_dir, genome_annotation, sample_ids, indexed_bam, control_junc_dir)

  ABERRANT_SPLICING.out.parsed_leafcutter_stats.set{ parsed_leafcutter_stats }
  ABERRANT_SPLICING.out.parsed_leafcutter_data.set{ parsed_leafcutter_data }
  ABERRANT_SPLICING.out.spot_pvals.set{ spot_pvals }
  ABERRANT_SPLICING.out.spot_dists.set{ spot_dists }

}