process IndexVcfGATK {

  // Indexing vcf file
  
  label 'variantcalling'

  publishDir "${projectDir}/${params.main_output_dir}/${params.variants_out}", mode: "copy", pattern: "*.tbi"

  input:
  tuple val(sample_id), path(vcf)

  output:
  tuple val(sample_id), path("*.tbi"), emit: vcf_index

  """
  # Index vcf
  gatk IndexFeatureFile -I ${vcf}
  """

}