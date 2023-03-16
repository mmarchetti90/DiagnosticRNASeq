process ParseSourceFile {

  label 'local'

  publishDir "${params.resources_dir}", mode: "copy", pattern: "*{txt,tsv}"
  
  input:
  path scripts_dir
  path metadata_file

  output:
  path "ReadsList.txt", emit: reads_list

  """
  python ${scripts_dir}/parse_source/parse_source_file.py --source_file ${metadata_file}
  """

}