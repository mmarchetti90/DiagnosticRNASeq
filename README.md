# Diagnostic RNASeq analysis
## Dockerized nextflow workflow for diagnostic RNASeq analysis with leafcutterMD, SPOT, OUTRIDER, and ASE

/// --------------------------------------- //

### RUN COMMAND:

nextflow run [OPTIONS] --source_file "/path/to/source/file"

### OPTIONS: (See config file for more)

--source_file

		Path to source file, a tab-separated txt document with lines sample_id\tfile_name
		Note that paired-end files are reported in separate lines, but have same sample_id
		(REQUIRED)

--genome_fasta_path

		Path to genome fasta file for STAR genomeGenerate
		Not needed if STAR index path is provided

--genome_annotation_path

		Path to genome annotation gtf file for STAR genomeGenerate
		Not needed if STAR index path is provided

--star_index_dir

		Path to STAR index
		If not provided, a new index will be generated

--control_cohort_junc_dir

		Path to junc files generated from a control cohort
		Can be omitted, in which case only provided samples will be used

--control_cohort_counts_dir

		Path to genes count files generated from a control cohort
		Can be omitted, in which case only provided samples will be used

--hpo_obo_path

		Path to OBO file for HPO terms
		Get from https://hpo.jax.org/data/ontology

--genes_to_phenotype_path

		Path to genes to HPO term annotation file
		Get from https://hpo.jax.org/data/annotations

--target_hpo_terms_path

		Diagnostic HPO terms for proband

--genomic_vcf_path

		Path to vcf file generated from DNA data
		If absent, the vcf generated from RNA will be used for ASE

--read_length

		Used for STAR index generation
		(Default, 50)

--saindexnbases

		Used for STAR index generation
		(Default, 14)

--strand

		For STAR --quantMode counts parsing. Possible values are unstranded, stranded, reversestrand
		(Default, unstranded)

--min_counts

		For OUTRIDER pre-filtering
		(Default, 10)

--pval_threshold

		For results filtering
		(Default, 0.05)

--ploidy

		Ploidy for GATK variant calling
		(Default,  2)

--variants_only

    	If true, HaplotypeCaller will run in EMIT_VARIANTS_ONLY mode
    	If false, HaplotypeCaller will run in EMIT_ALL_CONFIDENT_SITES mode
    	(Default, true)

--min_depth

		Min depth of variants for ASE analysis
    	(Default, 20)

/// --------------------------------------- ///

### DEPENDENCIES:

Nextflow 20.10+

Singularity 3.8.5+

Python 3.9.7+ &

	os
	pandas
	sys

leafcutter

	https://github.com/davidaknowles/leafcutter
	Copy bed2junc.pl, filter_cs.py, leafcutterMD.R, and sam2bed.pl to ./scripts/leafcutter

SPOT

	https://github.com/BennyStrobes/SPOT

OUTRIDER (included in the docker image)

	https://github.com/gagneurlab/OUTRIDER

/// --------------------------------------- ///

### NOTES:

- Python code from leafcutterMD (leafcutter_cluster.py) and SPOT (spot.py) was updated for Python3 compatibility;

- OUTRIDER is included in the docker image;

- SPOT was updated to use pystan 3;

- Useful scripts to extract results can be found in scripts/utils

/// --------------------------------------- ///

### RECOMMENDATIONS:

- For LeafcutterMD, if a trio is present, the control cohort could be omitted, but still best to have it to account for cases the aberrant splicing is "inherited" from a parent;

- For HPO interpolation, it's best to only use key HPO terms rather than secondary clinical manifestations;

- The final ranking is dependent on HPO terms, so candidates from WGS analyses may score low if the phenotype overlap is only partial or if the HPO annotation is lacking.
