process RunSTAR {
  
  // STAR alignment

  label 'star'

  publishDir "${projectDir}/${params.main_output_dir}/${params.bam_dir}", mode: "copy", pattern: "*{Unmapped,SJ}*"
  publishDir "${projectDir}/${params.main_output_dir}/${params.gene_counts_dir}", mode: "copy", pattern: "*_ReadsPerGene.out.tab"
  publishDir "${projectDir}/${params.main_output_dir}/${params.qc_dir}/${params.star_subdir}", mode: "copy", pattern: "*Log*"

  input:
  each path(index)
  tuple val(sample_id), path(read1), path(read2)

  output:
  tuple val(sample_id), path("${sample_id}_Aligned.out.bam"), emit: bam_files
  tuple val(sample_id), path("${sample_id}_SJ.out.tab"), emit: star_sj
  path "${sample_id}_Unmapped.out.mate1"
  path "${sample_id}_ReadsPerGene.out.tab", emit: gene_counts
  path "${sample_id}_Log.{final.out,out,progress.out}", emit: star_reports
  
  """
  # STAR alignment
  temp="${read1}"

  if [[ \${temp: -2} == "gz" ]]
  then

    if [[ "${read2}" == "mock.fastq" ]]
    then

      STAR \
      --runThreadN \$SLURM_CPUS_ON_NODE \
      --runMode alignReads \
      --twopassMode Basic \
      --genomeDir ${index} \
      --readFilesIn ${read1} \
      --outFileNamePrefix ./${sample_id}_ \
      --outSAMtype BAM Unsorted \
      --outReadsUnmapped Fastx \
      --quantMode GeneCounts \
      --readFilesCommand zcat

    else

      STAR \
      --runThreadN \$SLURM_CPUS_ON_NODE \
      --runMode alignReads \
      --twopassMode Basic \
      --genomeDir ${index} \
      --readFilesIn ${read1} ${read2} \
      --outFileNamePrefix ./${sample_id}_ \
      --outSAMtype BAM Unsorted \
      --outReadsUnmapped Fastx \
      --quantMode GeneCounts \
      --readFilesCommand zcat

    fi

  else

    if [[ "${read2}" == "mock.fastq" ]]
    then

      STAR \
      --runThreadN \$SLURM_CPUS_ON_NODE \
      --runMode alignReads \
      --twopassMode Basic \
      --genomeDir ${index} \
      --readFilesIn ${read1} \
      --outFileNamePrefix ./${sample_id}_ \
      --outSAMtype BAM Unsorted \
      --outReadsUnmapped Fastx \
      --quantMode GeneCounts

    else

      STAR \
      --runThreadN \$SLURM_CPUS_ON_NODE \
      --runMode alignReads \
      --twopassMode Basic \
      --genomeDir ${index} \
      --readFilesIn ${read1} ${read2} \
      --outFileNamePrefix ./${sample_id}_ \
      --outSAMtype BAM Unsorted \
      --outReadsUnmapped Fastx \
      --quantMode GeneCounts

    fi

  fi
  """

}