process IndexBam {

  // Index mapped bam file
  
  label 'samtools'

  publishDir "${projectDir}/${params.main_output_dir}/${params.bam_dir}", mode: "copy", pattern: "*.bai"

  input:
  tuple val(sample_id), path(bam)

  output:
  tuple val(sample_id), path("*.bai"), emit: bam_index

  """
  # Index bam
  samtools index \
  -b \
  -@ \$SLURM_CPUS_ON_NODE \
  ${bam}
  """

}