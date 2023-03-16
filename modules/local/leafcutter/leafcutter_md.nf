process LeafcutterMD {

    label 'slurm'

    publishDir "${projectDir}/${params.leafcuttermd_out}", mode: "copy", pattern: "*.txt"

    input:
    path scripts_dir
    path intron_counts

    output:
    path "leafcutter_outlier_pVals.txt"
    path "leafcutter_outlier_clusterPvals.txt"
    path "leafcutter_outlier_effSize.txt"

    """
    Rscript ${scripts_dir}/leafcutter/leafcutterMD.R --num_threads 4 ${intron_counts}
    """

}