process SubsetJuncControls {

  // Selects controls closest to the probandds for LeafcutterMD
  
  label 'python'

  publishDir "${projectDir}/${params.main_output_dir}/${params.qc_dir}/${params.leafcutter_ctrl_select_subdir}", mode: 'copy', pattern: "pca.pkl"
  publishDir "${projectDir}/${params.main_output_dir}/${params.qc_dir}/${params.leafcutter_ctrl_select_subdir}", mode: 'copy', pattern: "*.png"
  publishDir "${projectDir}/${params.main_output_dir}/${params.qc_dir}/${params.leafcutter_ctrl_select_subdir}", mode: 'copy', pattern: "pca_transform.tsv"
  publishDir "${projectDir}/${params.main_output_dir}/${params.qc_dir}/${params.leafcutter_ctrl_select_subdir}", mode: 'copy', pattern: "junc_selected_control_manifest.txt"

  input:
  path scripts_dir
  path junc_files
  path control_junc_dir

  output:
  path "pca_explained_variance.png", optional: true
  path "pca_transform.tsv", optional: true
  path "pca.pkl", optional: true
  path "controls_selection.png", optional: true
  path "junc_selected_control_manifest.txt", emit: selected_junc_controls

  """
  # Create manifest of proband samples
  ls -1 *.junc > samples_manifest.txt

  # Create manifest of control samples
  ls -1 ${control_junc_dir}/*.junc > control_manifest.txt

  # Filter controls
  python3 ${scripts_dir}/leafcutter/select_junc_control_cohort.py \
  --samples samples_manifest.txt \
  --controls control_manifest.txt \
  --min_reads ${params.leafcutter_ctrls_min_reads} \
  --max_ctrls ${params.leafcutter_ctrls_max}
  """

}