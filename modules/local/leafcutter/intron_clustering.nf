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

    # Adding junc files from control cohort (if any) and creating list of ctrl ids
    control_junc_num=\$(ls ${control_junc_dir} | grep ".junc" | wc -l)
    if (( \$control_junc_num > 0 ))
    then

        ls -1 ${control_junc_dir}/*.junc >> juncfiles.txt

    fi

    # Intron clustering
    python3 ${scripts_dir}/leafcutter/leafcutter_cluster.py -j juncfiles.txt -m 50 -o intronclust -l 500000

    # Batch correction
    mv intronclust_perind_numers.counts.gz raw_intronclust_perind_numers.counts.gz
    if (( \$control_junc_num > 0 ))
    then

        # Create list of ctrl_ids
        touch ctrl_ids.txt
        for ctrl in ${control_junc_dir}/*.junc
        do

            basename -s .junc \${ctrl} >> ctrl_ids.txt

        done

        # Batch correction
        Rscript ${scripts_dir}/leafcutter/intronclust_batch_correction.R --counts raw_intronclust_perind_numers.counts.gz --ctrl_ids_list ctrl_ids.txt
        mv corrected_intronclust_perind_numers.counts.gz intronclust_perind_numers.counts.gz

    fi
    """

}