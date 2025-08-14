/*
Merge analyses to find candidate genes
*/

// ----------------Workflow---------------- //

include { ToolsIntegration } from '../../modules/local/data_integration/tools_integration.nf'
include { HpoInterpolation as HpoInterpolationLeafcutter } from '../../modules/local/data_integration/hpo_interpolation.nf'
include { HpoInterpolation as HpoInterpolationOutrider } from '../../modules/local/data_integration/hpo_interpolation.nf'
include { HpoInterpolation as HpoInterpolationAse } from '../../modules/local/data_integration/hpo_interpolation.nf'
include { HpoInterpolation as HpoInterpolationIntegrated } from '../../modules/local/data_integration/hpo_interpolation.nf'

workflow FIND_CANDIDATES {

  take:
  scripts_dir
  hpo_obo
  genes_to_phenotype
  parsed_leafcutter_data
  parsed_outrider_data
  parsed_ase_data

  main:
  // JOIN TOOLS OUTPUTS ------------------- //
  
  parsed_leafcutter_data
    .join(parsed_outrider_data, by: 0, remainder: false)
    .join(parsed_ase_data, by: 0, remainder: false)
    .set{ tools_outputs }

  // TOOLS GENE RANKING ------------------- //

  ToolsIntegration(scripts_dir, tools_outputs)

  // HPO TERMS INTEGRATION ---------------- //

  if (new File("${params.target_hpo_terms_path}").exists()) {

    // Creating channel for target HPO terms
    Channel
      .fromPath("${params.target_hpo_terms_path}")
      .ifEmpty("mock.target_hpo_terms.txt")
      .set{ target_hpo_terms }

    // Aberrant splicing
    HpoInterpolationLeafcutter('aberrant_splicing', scripts_dir, hpo_obo, genes_to_phenotype, target_hpo_terms, parsed_leafcutter_data)

    hpo_ranking_leafcutter = HpoInterpolationLeafcutter.out.hpo_ranking

    // Aberrant expression
    HpoInterpolationOutrider('aberrant_expression', scripts_dir, hpo_obo, genes_to_phenotype, target_hpo_terms, parsed_outrider_data)

    hpo_ranking_outrider = HpoInterpolationOutrider.out.hpo_ranking

    // Allelic imbalance
    HpoInterpolationAse('allelic_imbalance', scripts_dir, hpo_obo, genes_to_phenotype, target_hpo_terms, parsed_ase_data)

    hpo_ranking_ase = HpoInterpolationAse.out.hpo_ranking

    // Integrated rankings
    HpoInterpolationIntegrated('integrated_data', scripts_dir, hpo_obo, genes_to_phenotype, target_hpo_terms, ToolsIntegration.out.tools_ranking)

    hpo_ranking_integrated = HpoInterpolationIntegrated.out.hpo_ranking

  }
  else {

    hpo_ranking_leafcutter = Channel.empty()

    hpo_ranking_outrider = Channel.empty()

    hpo_ranking_ase = Channel.empty()

    hpo_ranking_integrated = Channel.empty()

  }

  emit:
  tools_ranking = ToolsIntegration.out.tools_ranking
  hpo_ranking_leafcutter
  hpo_ranking_outrider
  hpo_ranking_ase
  hpo_ranking_integrated

}