process SortBam {
  
  // Samtools sort

  label 'samtools'

  publishDir "${projectDir}/${params.main_output_dir}/${params.bam_dir}", mode: "copy", pattern: "*_Aligned.sortedByCoord.out.bam"

  input:
  tuple val(sample_id), path(bam)

  output:
  tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam"), emit: bam_files
  
  """
  # Sorting by coordinates
  # N.B. I use samtools instead of STAR --outSAMtype BAM SortedByCoordinate because STAR was sometimes runnin into memory issues on the cluster
  samtools sort ${bam} -o ${sample_id}_Aligned.sortedByCoord.out.bam -@ \$SLURM_CPUS_ON_NODE
  """

}