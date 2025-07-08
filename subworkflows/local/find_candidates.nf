/*
Merge analyses to find candidate genes
*/

// ----------------Workflow---------------- //

include { ToolsIntegration } from '../../modules/local/data_integration/tools_integration.nf'
include { HpoInterpolation } from '../../modules/local/data_integration/hpo_interpolation.nf'

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

    HpoInterpolation(scripts_dir, hpo_obo, genes_to_phenotype, target_hpo_terms, ToolsIntegration.out.tools_ranking)

    hpo_ranking = HpoInterpolation.out.hpo_ranking

  }
  else {

    hpo_ranking = Channel.empty()

  }

  emit:
  tools_ranking = ToolsIntegration.out.tools_ranking
  hpo_ranking

}