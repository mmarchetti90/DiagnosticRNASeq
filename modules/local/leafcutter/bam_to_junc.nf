process BamToJunc {

    label 'slurm'

    publishDir "${projectDir}/${params.junc_dir}", mode: "copy", pattern: "*.junc"

    input:
    each path(scripts_dir)
    tuple val(sample_id), path(bam_file)

    output:
    path "*junc", emit: junc_file

    """
    samtools view ${bam_file} | python3 ${scripts_dir}/leafcutter/filter_cs.py | ${scripts_dir}/leafcutter/sam2bed.pl --use-RNA-strand - ${sample_id}.bed
    ${scripts_dir}/leafcutter/bed2junc.pl ${sample_id}.bed ${sample_id}.junc
    #rm ${sample_id}.bed
    """

}