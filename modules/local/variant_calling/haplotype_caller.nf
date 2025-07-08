process HaplotypeCaller {

  // HaplotypeCaller with GATK
  
  label 'variantcalling'

  publishDir "${projectDir}/${params.variants_out}", mode: "copy", pattern: "*_gatk.g.vcf.gz"

  input:
  each path(reference_fasta)
  each path(reference_index)
  each path(gatk_dict)
  tuple val(sample_id), path(bam), path(bam_index)

  output:
  tuple val(sample_id), path("${sample_id}_gatk.g.vcf.gz"), emit: gatk_gvcf

  """
  if [ ${params.variants_only} == false ]
  then 
  
    # Call variants with GATK, output GVCF
    gatk HaplotypeCaller \
    -R ${reference_fasta} \
    -ploidy ${params.ploidy} \
    -I ${bam} \
    -ERC BP_RESOLUTION \
    --output-mode EMIT_ALL_CONFIDENT_SITES \
    -O ${sample_id}_gatk.g.vcf.gz
  
  else

    # Call variants with GATK, output GVCF
    gatk HaplotypeCaller \
    -R ${reference_fasta} \
    -ploidy ${params.ploidy} \
    -I ${bam} \
    -ERC GVCF \
    --output-mode EMIT_VARIANTS_ONLY \
    -O ${sample_id}_gatk.g.vcf.gz
  
  fi
  """

}