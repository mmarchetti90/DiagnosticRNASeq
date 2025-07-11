process MakeIgvSjBed {

  // Create bed file of STAR SJ for IGV
  
  label 'variantcalling'

  publishDir "${projectDir}/${params.qc_dir}/${params.sj_bed_subdir}", mode: "copy", pattern: "*{.bed.gz,.bed.gz.tbi}"

  input:
  each path(scripts_dir)
  tuple val(sample_id), path(sj)

  output:
  tuple val(sample_id), path("${sample_id}_sj.bed.gz"), path("${sample_id}_sj.bed.gz.tbi"), emit: sj_bed

  """
  # Generate bed file
  python3 ${scripts_dir}/sj_to_bed/sj_to_bed.py \
  --file_path ${sj} \
  --filter > ${sample_id}_sj.bed

  # Bgzip
  bgzip -@ $SLURM_CPUS_ON_NODE ${sample_id}_sj.bed

  # Index
  tabix -f ${sample_id}_sj.bed.gz
  """

}