process AseAnalysisPython {

  // Binomial test for allelic imbalance
  
  label 'python'

  publishDir "${projectDir}/${params.main_output_dir}/${params.allelic_imbalance_out}", mode: "copy", pattern: "*_ase_{snp,gene}_stats.tsv.gz"

  input:
  each path(scripts_dir)
  each path(reference_annotation)
  tuple val(sample_id), path(counts)

  output:
  tuple val(sample_id), path("${sample_id}_ase_snp_stats.tsv.gz"), emit: ase_snp_stats
  tuple val(sample_id), path("${sample_id}_ase_gene_stats.tsv.gz"), emit: ase_gene_stats

  """
  python3 ${scripts_dir}/allelic_imbalance/allelic_imbalance.py \
  --ase_counts ${counts} \
  --gtf ${reference_annotation} \
  --min_depth ${params.min_depth}

  mv ase_snp_stats.tsv.gz ${sample_id}_ase_snp_stats.tsv.gz
  mv ase_gene_stats.tsv.gz ${sample_id}_ase_gene_stats.tsv.gz
  """

}