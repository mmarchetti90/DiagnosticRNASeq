process Outrider {

    // Run Outrider

    label 'outrider'

    publishDir "${projectDir}/${params.main_output_dir}/${params.outrider_out}", mode: 'copy', pattern: "*.{tsv,rds}"

    input:
    path scripts_dir
    path counts
    path ctrl_ids

    output:
    path "aberrant_expression.tsv", emit: outrider_expr
    path "outrider_analysis.rds", optional: true

    """
    if [[ "${ctrl_ids}" == "mock.txt" ]]
    then

        Rscript ${scripts_dir}/outrider/outrider.R \
        --counts ${counts} \
        --mincounts ${params.min_counts} \
        --p_thr ${params.pval_threshold} \
        --threads \$SLURM_CPUS_ON_NODE

    else

        Rscript ${scripts_dir}/outrider/outrider.R \
        --counts ${counts} \
        --mincounts ${params.min_counts} \
        --p_thr ${params.pval_threshold} \
        --threads \$SLURM_CPUS_ON_NODE \
        --ctrl_ids_list ${ctrl_ids}

    fi
    """

}