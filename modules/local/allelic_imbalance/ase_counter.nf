process AseCounter {

  // Count allele reads
  
  label 'variantcalling'

  publishDir "${projectDir}/${params.allelic_imbalance_out}", mode: "copy", pattern: "*_ASE.tsv"

  input:
  each path(reference_fasta)
  each path(reference_index)
  each path(gatk_dict)
  tuple val(sample_id), path(bam), path(bam_index), path(vcf), path(vcf_index)

  output:
  tuple val(sample_id), path("${sample_id}_ASE.tsv"), emit: ase_counts

  """
  gatk ASEReadCounter \
  -R ${reference_fasta} \
  -I ${bam} \
  -V ${vcf} \
  -O ${sample_id}_ASE.tsv
  """

}