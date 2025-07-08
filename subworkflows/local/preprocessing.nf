/*
Reads trimming, alignment, and generation of files for downstream analyses
*/

// ----------------Workflow---------------- //

include { ParseSourceFile } from '../../modules/local/parse_source/parse_source_file.nf'
include { TrimFastQ } from '../../modules/local/trimgalore/trimgalore.nf'
include { RunFastQC } from '../../modules/local/qc/run_fastqc.nf'
include { GenerateStarIndex } from '../../modules/local/star/generate_star_index.nf'
include { RunSTAR } from '../../modules/local/star/run_star.nf'
include { SortBam } from '../../modules/local/star/sort_bam.nf'
include { IndexBam } from '../../modules/local/indexing/index_bam.nf'
include { MergeCounts } from '../../modules/local/star/merge_counts.nf'
include { MakeIgvSjBed } from '../../modules/local/igv_sj_bed/igv_sj_bed.nf'
include { BamToBW } from '../../modules/local/bam_to_bw/bam_to_bw.nf'

workflow PREPROCESSING {

  take:
  raw_reads
  scripts_dir
  genome_fasta
  genome_annotation
  star_index
  control_gene_counts_dir

  main:
  // READS QC ----------------------------- //

  if (params.trim_reads == true) {
  
    // TRIMGALORE --------------------------- //

    // Trimming adapters
    TrimFastQ(raw_reads)

    star_input = TrimFastQ.out.trimmed_fastq_files

  }
  else {

    // FASTQC ------------------------------- //

    // Basic fastq files QC
    RunFastQC(raw_reads)

    star_input = raw_reads

  }

  // STAR ALIGNMENT ----------------------- //

  // Run STAR alignment
  RunSTAR(star_index, star_input)

  // Sort BAM
  SortBam(RunSTAR.out.bam_files)

  // Index BAM
  IndexBam(SortBam.out.bam_files)

  // Join BAM and BAI channels
  SortBam.out.bam_files
    .join(IndexBam.out.bam_index, by: 0, remainder: false)
    .set{ indexed_bam }

  // Merge gene count files
  MergeCounts(RunSTAR.out.gene_counts.collect(), control_gene_counts_dir)

  // GENERATING USEFUL FILES -------------- //

  // Create bed file of STAR SJ for IGV
  MakeIgvSjBed(scripts_dir, RunSTAR.out.star_sj)

  // Bam to bigWig
  BamToBW(indexed_bam)

  emit:
  indexed_bam
  merged_counts = MergeCounts.out.merged_counts
  control_gene_counts_ids = MergeCounts.out.control_gene_counts_ids

}