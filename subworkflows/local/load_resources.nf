/*
Load needed resources
*/

// ----------------Workflow---------------- //

include { IndexFastaSamtools } from '../../modules/local/indexing/index_fasta_samtools.nf'
include { GenerateStarIndex } from '../../modules/local/star/generate_star_index.nf'
include { MakeGATKDict } from '../../modules/local/variant_calling/make_gatk_dict.nf'
include { ParseSourceFile } from '../../modules/local/parse_source/parse_source_file.nf'

workflow LOAD_RESOURCES {

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

  // Index fasta with samtools
  IndexFastaSamtools(genome_fasta)

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

  // Create GATK dictionary
  MakeGATKDict(genome_fasta)

  // Control cohort data
  Channel
    .fromPath("${params.control_cohort_junc_dir}")
    .ifEmpty("${projectDir}/modules")
    .set{ control_junc_dir }

  Channel
    .fromPath("${params.control_cohort_counts_dir}")
    .ifEmpty("${projectDir}/modules")
    .set{ control_gene_counts_dir }
  
  // Creating channel for HPO obo file
  Channel
    .fromPath("${params.hpo_obo_path}")
    .ifEmpty("mock.hpo.obo")
    .set{ hpo_obo }

  // Creating channel for genes to phenotype file
  Channel
    .fromPath("${params.genes_to_phenotype_path}")
    .ifEmpty("mock.genes_to_phenotype.txt")
    .set{ genes_to_phenotype }

  // Parsing source file to output lists fastq files
  ParseSourceFile(scripts_dir, source_file)

  // Creating raw_reads channel
  ParseSourceFile.out.reads_list
    .splitCsv(header: true, sep: '\t')
    .map{ row -> tuple(row.SampleID, file(row.File1), file(row.File2)) }
    .set{ raw_reads }

  // Channel for sample IDs
  raw_reads
    .map{ it[0] }
    .set{ sample_ids }

  emit:
  scripts_dir
  genome_fasta
  genome_fasta_index = IndexFastaSamtools.out.reference_fasta_index
  genome_annotation
  star_index
  gatk_dict = MakeGATKDict.out.gatk_dict
  control_junc_dir
  control_gene_counts_dir
  hpo_obo
  genes_to_phenotype
  raw_reads
  sample_ids

}