process IntronClustering {

    label 'slurm'

    publishDir "${projectDir}/${params.intron_clusters}", mode: "copy", pattern: "intronclust_perind_numers.counts.gz"

    input:
    path scripts_dir
    path junc_files
    path control_junc_dir

    output:
    path "intronclust_perind_numers.counts.gz", emit: intron_counts

    """
    # Create list of junc files
    ls -1 *.junc > juncfiles.txt

    # Adding junc files from control cohort (if any)
    control_junc_num=\$(ls ${control_junc_dir} | grep ".junc" | wc -l)
    if (( \$control_junc_num > 0 ))
    then

        ls -1 ${control_junc_dir}/*.junc >> juncfiles.txt

    fi

    # intron clustering
    python3 ${scripts_dir}/leafcutter/leafcutter_cluster.py -j juncfiles.txt -m 50 -o intronclust -l 500000
    """

}