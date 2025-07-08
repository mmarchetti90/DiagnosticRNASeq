/*
Joint variant calling on RNA data
*/

// ----------------Workflow---------------- //

include { SplitNCigarReads } from '../../modules/local/variant_calling/split_n_cigar.nf'
include { ReplaceReadGroup } from '../../modules/local/variant_calling/replace_read_group.nf'
include { MarkDuplicates } from '../../modules/local/variant_calling/mark_duplicates.nf'
include { IndexBam } from '../../modules/local/indexing/index_bam.nf'
include { HaplotypeCaller } from '../../modules/local/variant_calling/haplotype_caller.nf'
include { IndexVcfGATK as IndexSingleGVCF } from '../../modules/local/indexing/index_vcf_gatk.nf'
include { CombineGVCFs } from '../../modules/local/variant_calling/combine_gvcfs.nf'
include { IndexVcfGATK as IndexCombinedGVCFs } from '../../modules/local/indexing/index_vcf_gatk.nf'
include { GenotypeGVCF } from '../../modules/local/variant_calling/genotype_gvcf.nf'
include { IndexVcfGATK as IndexJointVCF } from '../../modules/local/indexing/index_vcf_gatk.nf'

workflow VARIANT_CALLING {

  take:
  genome_fasta
  genome_fasta_index
  gatk_dict
  indexed_bam

  main:
  // VARIANT CALLING ---------------------- //

  // Fix CIGAR
  SplitNCigarReads(genome_fasta, genome_fasta_index, gatk_dict, indexed_bam)

  // Replace read groups
  ReplaceReadGroup(SplitNCigarReads.out.bam_files)

  // Mark duplicates
  MarkDuplicates(ReplaceReadGroup.out.bam_files)

  // Index mrkdup bam
  IndexBam(MarkDuplicates.out.bam_files)

  // Define HaplotypeCaller bam input
  MarkDuplicates.out.bam_files
    .join(IndexBam.out.bam_index, by: 0, remainder: false)
    .set{ mrkdup_indexed_bam }

  // HaplotypeCaller
  HaplotypeCaller(genome_fasta, genome_fasta_index, gatk_dict, mrkdup_indexed_bam)

  // Index GVCFs
  IndexSingleGVCF(HaplotypeCaller.out.gatk_gvcf)

  // Prep input for CombineGVCFs (sample IDs are discarded, only files are kept)
  HaplotypeCaller.out.gatk_gvcf
    .join(IndexSingleGVCF.out.vcf_index, by: 0, remainder: false)
    .map{ g -> tuple(file(g[1]), file(g[2])) }
    .collect()
    .set{ gvcf_files }

  // CombineGVCFs
  CombineGVCFs(genome_fasta, genome_fasta_index, gatk_dict, gvcf_files)

  // Index combined gvcf
  IndexCombinedGVCFs(CombineGVCFs.out.combined_gvcf)

  // Prep input for CombineGVCFs
  CombineGVCFs.out.combined_gvcf
    .join(IndexCombinedGVCFs.out.vcf_index, by: 0, remainder: false)
    .set{ indexed_combined_gvcf }

  // Genotyping
  GenotypeGVCF(genome_fasta, genome_fasta_index, gatk_dict, indexed_combined_gvcf)

  // Index VCF
  IndexJointVCF(GenotypeGVCF.out.gatk_vcf_unfilt)

  // Merge joint VCF and index
  GenotypeGVCF.out.gatk_vcf_unfilt
    .join(IndexJointVCF.out.vcf_index, by: 0, remainder: false)
    .set{ joint_rna_vcf }

  emit:
  mrkdup_indexed_bam
  joint_rna_vcf

}