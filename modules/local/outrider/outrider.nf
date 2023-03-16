process BamToJunc {

    label 'slurm'

    publishDir "${projectDir}/${params.outrider_out}", mode: 'copy', pattern: "*.{tsv,rds}"

    input:
    path scripts_dir
    path counts

    output:
    path "aberrant_expression.tsv"
    path "outrider_analysis.rds"

    """
    Rscript ${scripts_dir}/outrider/outrider.R --counts ${counts} --mincounts ${params.min_counts} --p_thr ${params.pval_threshold} --threads \$SLURM_CPUS_ON_NODE
    """

}