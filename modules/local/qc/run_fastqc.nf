process RunFastQC {
  
  label 'slurm'

  publishDir "${projectDir}/${params.fastqc_dir}", mode: "move", pattern: "*_fastqc.{html,zip}"

  input:
  tuple val(read_id), path(read1), path(read2)

  output:
  path "*_fastqc.{html,zip}", optional: true, emit: fastqc_reports

  """
  if [[ "${read2}" == "mock.fastq" ]]
  then

    # Single-end
    fastqc --outdir . -t \$SLURM_CPUS_ON_NODE ${read1}

  elif [[ "${workflow.profile}" != "cellranger" ]]
  then

    # Paired-end
    fastqc --outdir . -t \$SLURM_CPUS_ON_NODE ${read1} ${read2}

  else

    # Single-cell
    fastqc --outdir . -t \$SLURM_CPUS_ON_NODE ${read2}/*.fastq.gz

  fi
  """

}