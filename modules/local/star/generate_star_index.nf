process GenerateStarIndex {

  label 'slurm'

  input:
  path genome_fasta
  path genome_annotation
  val overhang

  output:
  path "StarIndex", emit: star_index

  script:
  """
  mkdir StarIndex

  STAR \
  --runThreadN \$SLURM_CPUS_ON_NODE \
  --runMode genomeGenerate \
  --genomeDir ./StarIndex/ \
  --genomeFastaFiles ${genome_fasta} \
  --sjdbGTFfile ${genome_annotation} \
  --sjdbOverhang ${overhang} \
  --genomeSAindexNbases ${params.saindexnbases}
  """

}