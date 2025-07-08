process MarkDuplicates {

  // Mark duplicates
  
  label 'variantcalling'

  publishDir "${projectDir}/${params.bam_dir}", mode: "copy", pattern: "*_mrkdup.bam"
  publishDir "${projectDir}/${params.qc_dir}/${params.duplication_subdir}", mode: "copy", pattern: "*_marked_dup_metrics.txt"

  input:
  tuple val(sample_id), path(bam)

  output:
  tuple val(sample_id), path("${sample_id}_mrkdup.bam"), emit: bam_files

  """
  # Mark duplicates
  gatk MarkDuplicates \
  -I ${bam} \
  -O ${sample_id}_mrkdup.bam  \
  -M ${sample_id}_marked_dup_metrics.txt
  """

}