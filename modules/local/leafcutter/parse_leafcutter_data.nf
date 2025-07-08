process ParseLeafcutter {

  // Parses the LeafcutterMD output files
  
  label 'python'

  publishDir "${projectDir}/${params.leafcuttermd_out}", mode: "copy", pattern: "*_splicing_stats.tsv.gz"
  publishDir "${projectDir}/${params.leafcuttermd_out}", mode: "copy", pattern: "*_abberrant_splicing.tsv.gz"

  input:
  each path(scripts_dir)
  each path(reference_annotation)
  each path(leafcutter_pvals)
  each path(leafcutter_effects)
  each path(leafcutter_cluster_pvals)
  val sample_id

  output:
  tuple val(sample_id), path("${sample_id}_splicing_stats.tsv.gz"), emit: parsed_leafcutter_stats
  tuple val(sample_id), path("${sample_id}_abberrant_splicing.tsv.gz"), emit: parsed_leafcutter_data

  """
  python3 ${scripts_dir}/leafcutter/parse_leafcutter_data.py \
  --sample ${sample_id} \
  --intron_pvals ${leafcutter_pvals} \
  --intron_effect ${leafcutter_effects} \
  --cluster_pvals ${leafcutter_cluster_pvals} \
  --gtf ${reference_annotation}
  """

}