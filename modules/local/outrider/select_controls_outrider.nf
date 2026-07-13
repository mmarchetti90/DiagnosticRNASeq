process SubsetGeneCountsControls {

  // Selects controls closest to the probandds for Outrider
  
  label 'python'

  publishDir "${projectDir}/${params.main_output_dir}/${params.qc_dir}/${params.outrider_ctrl_select_subdir}", mode: 'copy', pattern: "pca.pkl"
  publishDir "${projectDir}/${params.main_output_dir}/${params.qc_dir}/${params.outrider_ctrl_select_subdir}", mode: 'copy', pattern: "*.png"
  publishDir "${projectDir}/${params.main_output_dir}/${params.qc_dir}/${params.outrider_ctrl_select_subdir}", mode: 'copy', pattern: "pca_transform.tsv"
  publishDir "${projectDir}/${params.main_output_dir}/${params.gene_counts_dir}", mode: 'copy', pattern: "MergedGeneCounts_Subset.tsv"

  input:
  path scripts_dir
  path merged_counts
  path control_gene_counts_ids

  output:
  path "pca_explained_variance.png", optional: true
  path "pca_transform.tsv", optional: true
  path "pca.pkl", optional: true
  path "controls_selection.png", optional: true
  path "MergedGeneCounts_Subset.tsv", emit: merged_counts_filtered
  path "{ctrl_ids,mock}.filtered.txt", emit: control_gene_counts_ids_filtered

  """
  python3 ${scripts_dir}/outrider/select_gene_counts_control_cohort.py \
  --counts ${merged_counts} \
  --controls ctrls.list \
  --min_reads ${params.outrider_ctrls_min_reads} \
  --max_ctrls ${params.outrider_ctrls_max}
  """

}