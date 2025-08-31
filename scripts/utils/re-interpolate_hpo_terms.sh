#!/bin/bash

script=${1}
hpo_terms=${2}
hpo_obo=${3}
genes2hpo=${4}
analysis_dir=${5}
sample=${6}

# Init results dir

results_dir=reinterpolated_data
#rm -r ${results_dir}
mkdir -p ${results_dir}/6_data_integration

# Allelic imbalance

mkdir -p ${results_dir}/5_allelic_imbalance

genes_rank=${analysis_dir}/5_allelic_imbalance/${sample}_ase_gene_stats.tsv.gz

cp ${genes_rank} ${results_dir}/5_allelic_imbalance/

python ${script} --obo ${hpo_obo} --genes2hpo ${genes2hpo} --hpo ${hpo_terms} --genes_rank ${genes_rank} --rank_type allelic_imbalance

mv gene_ranks_with_hpo-allelic_imbalance.tsv.gz ${results_dir}/6_data_integration/${sample}_gene_ranks_with_hpo-allelic_imbalance.tsv.gz

# Leafcutter

mkdir -p ${results_dir}/5_leafcutter_md

genes_rank=${analysis_dir}/5_leafcutter_md/${sample}_aberrant_splicing.tsv.gz

cp ${genes_rank} ${results_dir}/5_leafcutter_md/

python ${script} --obo ${hpo_obo} --genes2hpo ${genes2hpo} --hpo ${hpo_terms} --genes_rank ${genes_rank} --rank_type aberrant_splicing

mv gene_ranks_with_hpo-aberrant_splicing.tsv.gz ${results_dir}/6_data_integration/${sample}_gene_ranks_with_hpo-aberrant_splicing.tsv.gz

# Outrider

mkdir -p ${results_dir}/5_outrider

genes_rank=${analysis_dir}/5_outrider/${sample}_aberrant_expression.tsv.gz

cp ${genes_rank} ${results_dir}/5_outrider/

python ${script} --obo ${hpo_obo} --genes2hpo ${genes2hpo} --hpo ${hpo_terms} --genes_rank ${genes_rank} --rank_type aberrant_expression

mv gene_ranks_with_hpo-aberrant_expression.tsv.gz ${results_dir}/6_data_integration/${sample}_gene_ranks_with_hpo-aberrant_expression.tsv.gz

# Integrated data

mkdir -p ${results_dir}/6_data_integration

genes_rank=${analysis_dir}/6_data_integration/${sample}_gene_ranks.tsv.gz

cp ${genes_rank} ${results_dir}/6_data_integration/

python ${script} --obo ${hpo_obo} --genes2hpo ${genes2hpo} --hpo ${hpo_terms} --genes_rank ${genes_rank} --rank_type integrated_data

mv gene_ranks_with_hpo-integrated_data.tsv.gz ${results_dir}/6_data_integration/${sample}_gene_ranks_with_hpo-integrated_data.tsv.gz
