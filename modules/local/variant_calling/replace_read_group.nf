process ReplaceReadGroup {

  // Replace read groups for post-processing with GATK
  
  label 'variantcalling'

  //publishDir "${projectDir}/${params.main_output_dir}/${params.bam_dir}", mode: "copy", pattern: "*_rg.bam"

  input:
  tuple val(sample_id), path(bam), path(bai)

  output:
  tuple val(sample_id), path("${sample_id}_rg.bam"), emit: bam_files

  """
  # Add/replace read groups for post-processing with GATK
  picard AddOrReplaceReadGroups \
  I=${bam} \
  O=${sample_id}_rg.bam \
  RGID=${sample_id} \
  RGLB=${sample_id} \
  RGPL=illumina \
  RGPU=${sample_id} \
  RGSM=${sample_id}
  """

}