process CombineGVCFs {

  // Combine GVCFs
  
  label 'variantcalling'

  publishDir "${projectDir}/${params.main_output_dir}/${params.variants_out}", mode: "copy", pattern: "joint_calling.g.vcf.gz"

  input:
  path reference_fasta
  path reference_index
  path gatk_dict
  path(gvcf_files)

  output:
  tuple val("joint_calling"), path("joint_calling.g.vcf.gz"), emit: combined_gvcf

  """
  # Collect GVCFs
  gvcf_files=" "
  for file in *.g.vcf.gz
  do

    gvcf_files="\${gvcf_files} --variant \${file}"

  done

  # Combine GVCFs
  gatk CombineGVCFs \
  -R ${reference_fasta} \
  \${gvcf_files} \
  -O joint_calling.g.vcf.gz
  """

}