process LeafcutterMD {

    // Runs LeafcutterMD

    label 'leafcutter'

    publishDir "${projectDir}/${params.main_output_dir}/${params.leafcuttermd_out}", mode: "copy", pattern: "*.txt"

    input:
    path scripts_dir
    path intron_counts

    output:
    path "leafcutter_outlier_pVals.txt", emit: leafcutter_pvals
    path "leafcutter_outlier_clusterPvals.txt", emit: leafcutter_cl_pvals
    path "leafcutter_outlier_effSize.txt", emit: leafcutter_effsize

    """
    Rscript ${scripts_dir}/leafcutter/leafcutterMD.R --num_threads \$SLURM_CPUS_ON_NODE ${intron_counts}
    """

}