process MakeGATKDict {

  // Generating GATK dictionary
  
  label 'variantcalling'

  publishDir "${projectDir}/${params.main_output_dir}/${params.resources_dir}/${params.gatk_dict_subdir}", mode: "copy", pattern: "*.dict"

  input:
  path reference_fasta

  output:
  path "*.dict", emit: gatk_dict

  """
  # Generate GATK dictionary
  gatk CreateSequenceDictionary -R ${reference_fasta}
  """

}