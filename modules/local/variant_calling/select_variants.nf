process SelectVariants {

  // Select variants with GATK for ASE
  
  label 'variantcalling'

  publishDir "${projectDir}/${params.main_output_dir}/${params.variants_out}", mode: "copy", pattern: "*_select.vcf.gz"

  input:
  each path(reference_fasta)
  each path(reference_index)
  each path(gatk_dict)
  each path(vcf)
  each path(vcf_index)
  val sample_id

  output:
  tuple val(sample_id), path("${sample_id}_select.vcf.gz"), emit: filtered_vcf

  """
  gatk SelectVariants \
  -R ${reference_fasta} \
  -V ${vcf} \
  --restrict-alleles-to BIALLELIC \
  -select 'vc.getGenotype("${sample_id}").isHet()' \
  --select-type-to-include SNP \
  --sample-name ${sample_id} \
  -O ${sample_id}_tmp.vcf

  bcftools norm \
  --rm-dup all \
  -o ${sample_id}_select.vcf \
  -O v \
  ${sample_id}_tmp.vcf

  bgzip ${sample_id}_select.vcf
  """

}