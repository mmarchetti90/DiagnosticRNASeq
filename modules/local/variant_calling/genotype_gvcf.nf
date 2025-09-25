process GenotypeGVCF {

  // Genotyping with GATK
  
  label 'variantcalling'

  publishDir "${projectDir}/${params.main_output_dir}/${params.variants_out}", mode: "copy", pattern: "*_gatk_unfilt.vcf.gz"

  input:
  each path(reference_fasta)
  each path(reference_index)
  each path(gatk_dict)
  tuple val(sample_id), path(gvcf), path(gvcf_index)

  output:
  tuple val(sample_id), path("${sample_id}_gatk_unfilt.vcf.gz"), emit: gatk_vcf_unfilt

  """
  if [ ${params.variants_only} == false ]
  then 
  
    # GVCF to VCF. Min base quality score is 10 by default. Including non-variant sites in order to differentiate between consensus call and no-call sites
    gatk GenotypeGVCFs \
    -R ${reference_fasta} \
    --variant ${gvcf} \
    -ploidy ${params.ploidy} \
    --include-non-variant-sites true \
    --output ${sample_id}_gatk_unfilt.vcf.gz
  
  else

    # GVCF to VCF. Min base quality score is 10 by default.
    gatk GenotypeGVCFs \
    -R ${reference_fasta} \
    --variant ${gvcf} \
    -ploidy ${params.ploidy} \
    --include-non-variant-sites false \
    --output ${sample_id}_gatk_unfilt.vcf.gz
  
  fi
  """

}