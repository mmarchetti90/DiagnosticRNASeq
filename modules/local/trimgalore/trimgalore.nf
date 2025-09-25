process TrimFastQ {

  // Trim reads

  label 'trimgalore'

  publishDir "${projectDir}/${params.qc_dir}/${params.main_output_dir}/${params.trimming_subdir}", mode: "copy", pattern: "*_trimming_report.txt"
  publishDir "${projectDir}/${params.qc_dir}/${params.main_output_dir}/${params.fastqc_subdir}", mode: "copy", pattern: "*_fastqc.{html,zip}"

  input:
  tuple val(sample_id), path(read1), path(read2)

  output:
  path "*_trimming_report.txt", optional: true, emit: trimming_reports
  path "*_fastqc.{html,zip}", optional: true, emit: fastqc_reports
  tuple val(sample_id), path("{${sample_id}_val_1.fq.gz,${sample_id}_trimmed.fq.gz}"), path("{${sample_id}_val_2.fq.gz,mock.fastq}"), optional: true, emit: trimmed_fastq_files

  """
  if [[ "${read2}" == "mock.fastq" ]]
  then

      trim_galore \
      --cores 4 \
      --output_dir . \
      --basename ${sample_id} \
      --gzip \
      --fastqc \
      ${params.trimgalore_params} \
      ${read1}

      # Adding mock read2 output
      touch mock.fastq

  else

      trim_galore \
      --cores 4 \
      --output_dir . \
      --basename ${sample_id} \
      --gzip \
      --fastqc \
      ${params.trimgalore_params} \
      --paired \
      ${read1} ${read2}

  fi
  """

}