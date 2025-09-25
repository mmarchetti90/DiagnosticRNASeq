process IndexFastaSamtools {

  // Index reference fasta
  
  label 'samtools'

  publishDir "${projectDir}/${params.main_output_dir}/${params.resources_dir}/${params.fasta_index_subdir}", mode: "copy", pattern: "*.fai"

  input:
  path reference_fasta

  output:
  path "*.fai", emit: reference_fasta_index

  """
  # Index reference fasta
  samtools faidx ${reference_fasta}
  """

}