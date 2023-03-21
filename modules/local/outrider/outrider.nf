process Outrider {

    label 'slurm'

    publishDir "${projectDir}/${params.outrider_out}", mode: 'copy', pattern: "*.{tsv,rds}"

    input:
    path scripts_dir
    path counts
    path ctrl_ids

    output:
    path "aberrant_expression.tsv"
    path "outrider_analysis.rds"

    """
    if [[ "${ctrl_ids}" == "nextflow.config" ]]
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
