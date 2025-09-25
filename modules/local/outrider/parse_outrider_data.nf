process ParseOutrider {

  // Parses the Outrider output file
  
  label 'python'

  publishDir "${projectDir}/${params.main_output_dir}/${params.outrider_out}", mode: "copy", pattern: "*_aberrant_expression.tsv.gz"

  input:
  each path(scripts_dir)
  each path(reference_annotation)
  each path(raw_outrider_data)
  val sample_id

  output:
  tuple val(sample_id), path("${sample_id}_aberrant_expression.tsv.gz"), emit: parsed_outrider_data

  """
  python3 ${scripts_dir}/outrider/parse_outrider_data.py \
  --sample ${sample_id} \
  --outrider_data ${raw_outrider_data} \
  --gtf ${reference_annotation}
  """

}