# Diagnostic RNASeq analysis
## Dockerized nextflow workflow for diagnostic RNASeq analysis with leafcutterMD, SPOT, and OUTRIDER

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
	Copy dirichlet_multinomial.stan and dirichlet_multinomial.py to ./scripts/spot

OUTRIDER (included in the docker image)

	https://github.com/gagneurlab/OUTRIDER

/// --------------------------------------- ///

### NOTES:

- Python code from leafcutterMD (leafcutter_cluster.py) and SPOT (spot.py) was updated for Python3 compatibility;

- OUTRIDER is included in the docker image;

- SPOT does not currently work since it still uses pystan 2, while the docker image contains pystan 3;

