process RunFastQC {
  
  // FastQC

  label 'fastqc'

  publishDir "${projectDir}/${params.qc_dir}/${params.fastqc_subdir}", mode: "copy", pattern: "*_fastqc.{html,zip}"

  input:
  tuple val(sample_id), path(read1), path(read2)

  output:
  path "*_fastqc.{html,zip}", optional: true, emit: fastqc_reports

  """
  if [[ "${read2}" == "mock.fastq" ]]
  then

    # Single-end
    fastqc --outdir . -t \$SLURM_CPUS_ON_NODE ${read1}

  else

    # Paired-end
    fastqc --outdir . -t \$SLURM_CPUS_ON_NODE ${read1} ${read2}

  fi
  """

}