process ToolsIntegration {

  // Combines the output files from aberrant splicing, aberrant expression, and allelic imbalance to rank genes
  
  label 'python'

  publishDir "${projectDir}/${params.data_integration_out}", mode: "copy", pattern: "*_gene_ranks.tsv.gz"

  input:
  each path(scripts_dir)
  tuple val(sample_id), path(aberrant_splicing_data), path(aberrant_expression_data), path(allelic_imbalance_data)

  output:
  tuple val(sample_id), path("${sample_id}_gene_ranks.tsv.gz"), emit: tools_ranking

  """
  python3 ${scripts_dir}/data_integration/tools_integration.py \
  --aberrant_splicing ${aberrant_splicing_data} \
  --aberrant_expression ${aberrant_expression_data} \
  --allelic_imbalance ${allelic_imbalance_data}

  mv gene_ranks.tsv.gz ${sample_id}_gene_ranks.tsv.gz
  """

}