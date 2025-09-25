process BamToBW {

  // BAM to bigWig using deepTools
  
  label 'deeptools'

  publishDir "${projectDir}/${params.main_output_dir}/${params.qc_dir}/${params.bw_subdir}", mode: "copy", pattern: "*.{bw,log}"

  input:
  tuple val(sample_id), path(bam), path(bai)

  output:
  tuple val(sample_id), path("${sample_id}.bw"), emit: bw
  path "${sample_id}_bw.log"

  """
  bamCoverage \
  -b ${bam} \
  -of bigwig \
  -o ${sample_id}.bw \
  &> ${sample_id}_bw.log
  """

}