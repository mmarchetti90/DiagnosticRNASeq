/*
Workflow for diagnostic rnaseq analyses
*/

// ----------------Workflow---------------- //

include { ParseSourceFile } from '../../modules/local/parse_source/parse_source_file.nf'
include { RunFastQC } from '../../modules/local/qc/run_fastqc.nf'
include { GenerateStarIndex } from '../../modules/local/star/generate_star_index.nf'
include { RunSTAR } from '../../modules/local/star/run_star.nf'
include { MergeCounts } from '../../modules/local/star/merge_counts.nf'
include { BamToJunc } from '../../modules/local/leafcutter/bam_to_junc.nf'
include { IntronClustering } from '../../modules/local/leafcutter/intron_clustering.nf'
include { LeafcutterMD } from '../../modules/local/leafcutter/leafcutter_md.nf'
include { Spot } from '../../modules/local/spot/spot.nf'
include { Outrider } from '../../modules/local/outrider/outrider.nf'

workflow RNA_DIAGNOSTIC {

  take:
  source_file

  main:
  // LOADING RESOURCES -------------------- //

  // Channel for the directory containing the scripts used by the pipeline
  Channel
    .fromPath("${projectDir}/scripts")
    .set{ scripts_dir }

  // Channel for genome fasta
  Channel
    .fromPath("${params.genome_fasta_path}")
    .set{ genome_fasta }

  // Channel for genome annotation
  Channel
    .fromPath("${params.genome_annotation_path}")
    .set{ genome_annotation }

  // Creating channel for existing star index, or building de novo
  if (new File("${params.star_index_dir}/Genome").exists()) {

    Channel
    .fromPath("${params.star_index_dir}")
    .set{ star_index }

  }
  else {

    sjdboverhang = params.read_length - 1
    GenerateStarIndex(genome_fasta, genome_annotation, sjdboverhang)
    star_index = GenerateStarIndex.out.star_index

  }

  // Control cohort data
  Channel
    .fromPath("${params.control_cohort_junc_dir}")
    .ifEmpty("${projectDir}/modules")
    .set{ control_junc_dir }

  Channel
    .fromPath("${params.control_cohort_counts_dir}")
    .ifEmpty("${projectDir}/modules")
    .set{ control_gene_counts_dir }
  
  Channel
    .fromPath("${params.control_cohort_counts_ids}")
    .ifEmpty("${projectDir}/nextflow.config")
    .set{ control_gene_counts_ids }
  
  // Parsing source file to output lists fastq files
  ParseSourceFile(scripts_dir, source_file)

  // Creating raw_reads channel
  ParseSourceFile.out.reads_list
    .splitCsv(header: true, sep: '\t')
    .map{row -> tuple(row.SampleID, file(row.File1), file(row.File2))}
    .set{ raw_reads }

  // FASTQC ------------------------------- //

  // Checking reads quality with FasQC
  RunFastQC(raw_reads)

  // STAR ALIGNMENT ----------------------- //

  // Run STAR alignment
  RunSTAR(star_index, raw_reads)

  // Merge gene count files
  MergeCounts(RunSTAR.out.gene_counts.collect(), control_gene_counts_dir)

  // BAM TO JUNC -------------------------- //

  BamToJunc(scripts_dir, RunSTAR.out.bam_files)

  // INTRON CLUSTERING -------------------- //

  IntronClustering(scripts_dir, BamToJunc.out.junc_file.collect(), control_junc_dir)

  // LEAFCUTTER MD ------------------------ //

  LeafcutterMD(scripts_dir, IntronClustering.out.intron_counts)

  // SPOT --------------------------------- //

  Spot(scripts_dir, IntronClustering.out.intron_counts)

  // OUTRIDER ---------------------------- //

  Outrider(scripts_dir, MergeCounts.out.merged_counts, control_gene_counts_ids)

}
