/*
Allelic imbalance
*/

// ----------------Workflow---------------- //

include { IndexVcfGATK as IndexGenomicVcf } from '../../modules/local/indexing/index_vcf_gatk.nf'
include { SelectVariants } from '../../modules/local/variant_calling/select_variants.nf'
include { IndexVcfGATK as IndexSubVcf } from '../../modules/local/indexing/index_vcf_gatk.nf'
include { AseCounter } from '../../modules/local/allelic_imbalance/ase_counter.nf'
include { AseAnalysisPython } from '../../modules/local/allelic_imbalance/ase_analysis_python.nf'

workflow ALLELIC_IMBALANCE {

  take:
  scripts_dir
  genome_fasta
  genome_fasta_index
  genome_annotation
  gatk_dict
  sample_ids
  mrkdup_indexed_bam
  joint_rna_vcf

  main:
  // SELECT VARIANTS ---------------------- //

  // If a genomic VCF is provided, it will be preferred, otherwise the RNA VCF will be used instead
  if (new File("${params.genomic_vcf_path}").exists()) {

    Channel
      .fromPath("${params.genomic_vcf_path}")
      .map{ g -> tuple("genomic_joint_calling", file(g)) }
      .set{ joint_genomic_vcf }

    IndexGenomicVcf(joint_genomic_vcf)

    joint_genomic_vcf
      .map{ it[1] }
      .set{ joint_vcf }

    IndexGenomicVcf.out.vcf_index
      .map{ it[1] }
      .set{ joint_vcf_index }

  }
  else {

    joint_rna_vcf
      .map{ it[1] }
      .set{ joint_vcf }

    joint_rna_vcf
      .map{ it[2] }
      .set{ joint_vcf_index }

  }

  // Select variants
  SelectVariants(genome_fasta, genome_fasta_index, gatk_dict, joint_vcf, joint_vcf_index, sample_ids)

  // Index
  IndexSubVcf(SelectVariants.out.filtered_vcf)

  // Join BAM, VCF, and VCF index channels
  SelectVariants.out.filtered_vcf
    .join(IndexSubVcf.out.vcf_index, by: 0, remainder: false)
    .set{ subset_indexed_vcf }

  // ALLELIC IMBALANCE -------------------- //

  // Join BAM and VCF channels
  mrkdup_indexed_bam
    .join(subset_indexed_vcf, by: 0, remainder: false)
    .set{ ase_counter_input }

  // ASEReadCounter

  AseCounter(genome_fasta, genome_fasta_index, gatk_dict, ase_counter_input)

  // Run ASE analysis

  AseAnalysisPython(scripts_dir, genome_annotation, AseCounter.out.ase_counts)

  emit:
  ase_snp_stats = AseAnalysisPython.out.ase_snp_stats
  ase_gene_stats = AseAnalysisPython.out.ase_gene_stats

}