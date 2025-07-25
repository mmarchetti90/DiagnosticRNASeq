process HpoInterpolation {

  // Overlaps ranked genes to HPO terms
  
  label 'python'

  publishDir "${projectDir}/${params.data_integration_out}", mode: "copy", pattern: "*_gene_ranks_with_hpo.tsv.gz"

  input:
  each path(scripts_dir)
  each path(hpo_obo)
  each path(genes_to_phenotype)
  each path(target_hpo_terms)
  tuple val(sample_id), path(ranks)

  output:
  tuple val(sample_id), path("${sample_id}_gene_ranks_with_hpo.tsv.gz"), emit: hpo_ranking

  """
  python3 ${scripts_dir}/data_integration/hpo_interpolation.py \
  --obo ${hpo_obo} \
  --hpo ${target_hpo_terms} \
  --genes2hpo ${genes_to_phenotype} \
  --genes_rank ${ranks}

  mv gene_ranks_with_hpo.tsv.gz ${sample_id}_gene_ranks_with_hpo.tsv.gz
  """

}