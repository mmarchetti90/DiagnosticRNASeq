process SplitNCigarReads {

  // Corrects the CIGAR of RNA bams
  
  label 'variantcalling'

  publishDir "${projectDir}/${params.main_output_dir}/${params.bam_dir}", mode: "copy", pattern: "*.{bam,bai}"

  input:
  each path(reference_fasta)
  each path(reference_index)
  each path(gatk_dict)
  tuple val(sample_id), path(bam), path(bai)

  output:
  tuple val(sample_id), path("${sample_id}_cigarfixed.bam"), path("${sample_id}_cigarfixed.bam.bai"), emit: bam_files

  """
  # Fix CIGAR
  gatk SplitNCigarReads \
  -R ${reference_fasta} \
  -I ${bam} \
  -O ${sample_id}_cigarfixed.bam \
  --read-index ${bai} \
  --create-output-bam-index

  # Rename index
  mv ${sample_id}_cigarfixed.bai ${sample_id}_cigarfixed.bam.bai
  """

}