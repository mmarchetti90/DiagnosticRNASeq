process MergeCounts {

    label 'slurm'

    publishDir "${projectDir}/${params.gene_counts_dir}", mode: 'copy', pattern: "MergedGeneCounts.tsv"

    input:
    path count_files
    path control_gene_counts_dir

    output:
    path "MergedGeneCounts.tsv", emit: merged_counts
    path "{ctrl_ids,mock}.txt", emit: control_gene_counts_ids

    """
    # Get info
    header_row=4
    file_pattern=_ReadsPerGene.out.tab
    element_type=GeneID
    output_file=MergedGeneCounts.tsv
    element_id_column=1

    if [[ "${params.strand}" == "unstranded" ]]
    then

        counts_column=2

    elif [[ "${params.strand}" == "stranded" ]]
    then

        counts_column=3

    else

        counts_column=4
            
    fi

    # N.B. This script assumes that all files processed have the same elements and that they are in the same order.

    # Extract element ids
    echo "Extracting element ids"
    ls -1 *\${file_pattern} | head -n 1 | xargs awk -v header_row=\${header_row} '(NR > header_row) { print \$1 }' > merged_temp.txt

    # Extracting counts and merging files
    echo "Merging count files"
    for file in *\${file_pattern}
    do

        awk -v header_row=\${header_row} -v counts_column=\${counts_column} '(NR > header_row) { print \$counts_column }' \$file > counts_temp.txt

        paste merged_temp.txt counts_temp.txt > new_merged_temp.txt

        rm counts_temp.txt
        rm merged_temp.txt
        mv new_merged_temp.txt merged_temp.txt

    done

    # Adding counts from control cohort and creating a list of control ids (if any)
    control_counts_num=\$(ls ${control_gene_counts_dir} | grep "\${file_pattern}" | wc -l)
    if (( \$control_counts_num > 0 ))
    then

        # Adding counts to counts matrix
        for file in ${control_gene_counts_dir}/*\${file_pattern}
        do

            awk -v header_row=\${header_row} -v counts_column=\${counts_column} '(NR > header_row) { print \$counts_column }' \$file > counts_temp.txt

            paste merged_temp.txt counts_temp.txt > new_merged_temp.txt

            rm counts_temp.txt
            rm merged_temp.txt
            mv new_merged_temp.txt merged_temp.txt

        done

        # Creating list of control ids
        touch ctrl_ids.txt
        for file in ${control_gene_counts_dir}/*\${file_pattern}
        do

            sample_name=\$(basename \$file)
            sample_name=\${sample_name%\$file_pattern*}
            echo \${sample_name} >> ctrl_ids.txt

        done
    
    else
    
        touch mock.txt

    fi

    # Removing duplicate lines. Duplicates happen if the GTF used for featureCounts had duplicated entries (common for exons, apparently)
    echo "Removing duplicated lines"
    sort merged_temp.txt | uniq >> merged_temp_unique.txt
    rm merged_temp.txt

    # Add header
    echo "Adding header"
    header=\${element_type}

    for file in *\${file_pattern}
    do

        sample_name=\${file%\$file_pattern*}
        header="\${header}\t\${sample_name}"

    done

    if (( \$control_counts_num > 0 ))
    then

        for file in ${control_gene_counts_dir}/*\${file_pattern}
        do

            sample_name=\$(basename \$file)
            sample_name=\${sample_name%\$file_pattern*}
            header="\${header}\t\${sample_name}"

        done

    fi

    touch \${output_file}
    echo -e "\${header}" >> \${output_file}
    cat merged_temp_unique.txt >> \${output_file}
    rm merged_temp_unique.txt
    """

}
