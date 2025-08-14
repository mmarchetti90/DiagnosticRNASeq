/*
Workflow for diagnostic rnaseq analyses
*/

// ----------------Workflow---------------- //

include { LOAD_RESOURCES } from '../../subworkflows/local/load_resources.nf'
include { PREPROCESSING } from '../../subworkflows/local/preprocessing.nf'
include { ABERRANT_SPLICING } from '../../subworkflows/local/aberrant_splicing.nf'
include { ABERRANT_EXPRESSION } from '../../subworkflows/local/aberrant_expression.nf'
include { VARIANT_CALLING } from '../../subworkflows/local/variant_calling.nf'
include { ALLELIC_IMBALANCE } from '../../subworkflows/local/allelic_imbalance.nf'
include { FIND_CANDIDATES } from '../../subworkflows/local/find_candidates.nf'

workflow RNA_DIAGNOSTIC {

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
  LOAD_RESOURCES.out.gatk_dict.set{ gatk_dict }
  LOAD_RESOURCES.out.control_junc_dir.set{ control_junc_dir }
  LOAD_RESOURCES.out.control_gene_counts_dir.set{ control_gene_counts_dir }
  LOAD_RESOURCES.out.hpo_obo.set{ hpo_obo }
  LOAD_RESOURCES.out.genes_to_phenotype.set{ genes_to_phenotype }
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

  // ABERRANT EXPRESSION ------------------ //

  ABERRANT_EXPRESSION(scripts_dir, genome_annotation, sample_ids, merged_counts, control_gene_counts_ids)

  ABERRANT_EXPRESSION.out.parsed_outrider_data.set{ parsed_outrider_data }

  // VARIANT CALLING ---------------------- //

  VARIANT_CALLING(genome_fasta, genome_fasta_index, gatk_dict, indexed_bam)

  VARIANT_CALLING.out.mrkdup_indexed_bam.set{ mrkdup_indexed_bam }
  VARIANT_CALLING.out.joint_rna_vcf.set{ joint_rna_vcf }

  // ALLELIC IMBALANCE -------------------- //

  ALLELIC_IMBALANCE(scripts_dir, genome_fasta, genome_fasta_index, genome_annotation, gatk_dict, sample_ids, mrkdup_indexed_bam, joint_rna_vcf)

  ALLELIC_IMBALANCE.out.ase_snp_stats.set{ ase_snp_stats }
  ALLELIC_IMBALANCE.out.ase_gene_stats.set{ ase_gene_stats }

  // FIND CANDIDATE GENES ----------------- //

  FIND_CANDIDATES(scripts_dir, hpo_obo, genes_to_phenotype, parsed_leafcutter_data, parsed_outrider_data, ase_gene_stats)

  FIND_CANDIDATES.out.tools_ranking.set{ tools_ranking }
  FIND_CANDIDATES.out.hpo_ranking_leafcutter.set{ hpo_ranking_leafcutter }
  FIND_CANDIDATES.out.hpo_ranking_outrider.set{ hpo_ranking_outrider }
  FIND_CANDIDATES.out.hpo_ranking_ase.set{ hpo_ranking_ase }
  FIND_CANDIDATES.out.hpo_ranking_integrated.set{ hpo_ranking_integrated }

}