process Spot {

    label 'slurm'

    publishDir "${projectDir}/${params.spot_out}", mode: "copy", pattern: "spot_*.txt"

    input:
    path scripts_dir
    path junc_files

    output:
    path "spot_empirical_pvalue.txt"
    path "spot_md.txt"

    """
    python3 ${scripts_dir}/spot/spot.py --juncfile ${junc_files} --chains 4
    """

}
