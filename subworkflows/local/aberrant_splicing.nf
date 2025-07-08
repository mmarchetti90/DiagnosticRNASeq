/*
Aberrant splicing with LeafcutterMD and Spot
*/

// ----------------Workflow---------------- //

include { BamToJunc } from '../../modules/local/leafcutter/bam_to_junc.nf'
include { IntronClustering } from '../../modules/local/leafcutter/intron_clustering.nf'
include { LeafcutterMD } from '../../modules/local/leafcutter/leafcutter_md.nf'
include { ParseLeafcutter } from '../../modules/local/leafcutter/parse_leafcutter_data.nf'
include { Spot } from '../../modules/local/spot/spot.nf'

workflow ABERRANT_SPLICING {

  take:
  scripts_dir
  genome_annotation
  sample_ids
  bam_files
  control_junc_dir

  main:
  // ABERRANT SPLICING -------------------- //

  // BAM to JUNC
  BamToJunc(scripts_dir, bam_files)

  // Intron clustering
  IntronClustering(scripts_dir, BamToJunc.out.junc_file.collect(), control_junc_dir)

  // LeafcutterMD
  LeafcutterMD(scripts_dir, IntronClustering.out.intron_counts)

  // Parse LeafcutterMD output
  ParseLeafcutter(scripts_dir, genome_annotation, LeafcutterMD.out.leafcutter_pvals, LeafcutterMD.out.leafcutter_cl_pvals, LeafcutterMD.out.leafcutter_effsize, sample_ids)

  // SPOT
  Spot(scripts_dir, IntronClustering.out.intron_counts)

  emit:
  parsed_leafcutter_stats = ParseLeafcutter.out.parsed_leafcutter_stats
  parsed_leafcutter_data = ParseLeafcutter.out.parsed_leafcutter_data
  spot_pvals = Spot.out.spot_pvals
  spot_dists = Spot.out.spot_dists

}