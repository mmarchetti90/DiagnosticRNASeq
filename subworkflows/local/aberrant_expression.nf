/*
Aberrant expression with Outrider
*/

// ----------------Workflow---------------- //

include { Outrider } from '../../modules/local/outrider/outrider.nf'
include { ParseOutrider } from '../../modules/local/outrider/parse_outrider_data.nf'

workflow ABERRANT_EXPRESSION {

  take:
  scripts_dir
  genome_annotation
  sample_ids
  merged_counts
  control_gene_counts_ids

  main:
  // ABERRANT EXPRESSION ------------------ //

  // Run Outrider
  Outrider(scripts_dir, merged_counts, control_gene_counts_ids)

  // Filter Outrider output for individual samples
  ParseOutrider(scripts_dir, genome_annotation, Outrider.out.outrider_expr, sample_ids)

  emit:
  parsed_outrider_data = ParseOutrider.out.parsed_outrider_data

}