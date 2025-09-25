process Spot {

    // Runs Spot

    label 'leafcutter'

    publishDir "${projectDir}/${params.main_output_dir}/${params.spot_out}", mode: "copy", pattern: "spot_*.txt"

    input:
    path scripts_dir
    path junc_files

    output:
    path "spot_empirical_pvalue.txt", emit: spot_pvals
    path "spot_md.txt", emit: spot_dists

    """
    python3 ${scripts_dir}/spot/spot.py --juncfile ${junc_files} --chains 4
    """

}