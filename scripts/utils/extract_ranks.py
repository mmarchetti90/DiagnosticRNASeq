#!/usr/bin/env python3

"""
Script to conveniently extract ranks for a gene of interest

Parameters
----------
results_dir : str
    Path to results directory (can be .)
target_gene : str
    Gene symbol for which to extract ranks
"""

### ---------------------------------------- ###

### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd

from os import listdir
from sys import argv

### Parse args

main_dir = argv[argv.index('--results_dir') + 1]

target_gene = argv[argv.index('--target_gene') + 1]

### Extract rankings

rank_type_pval_columns = {'integrated_data' : 'final_pval',
                          'aberrant_splicing' : 'clusters_min_pval',
                          'aberrant_expression' : 'padj',
                          'allelic_imbalance' : 'compound_p'}

rank_type_effect_columns = {'integrated_data' : 'final_effect',
                            'aberrant_splicing' : 'intron_median_abs_effect',
                            'aberrant_expression' : 'log2FC',
                            'allelic_imbalance' : 'median_ase_score'}

subdirs_and_patterns = {'5_allelic_imbalance' : ['_ase_gene_stats.tsv.gz'],
                        '5_leafcutter_md' : ['_aberrant_splicing.tsv.gz'],
                        '5_outrider' : ['_aberrant_expression.tsv.gz'],
                        '6_data_integration' : ['_gene_ranks_with_hpo-aberrant_expression.tsv.gz',
                                                '_gene_ranks_with_hpo-aberrant_splicing.tsv.gz',
                                                '_gene_ranks_with_hpo-allelic_imbalance.tsv.gz',
                                                '_gene_ranks_with_hpo-integrated_data.tsv.gz',
                                                '_gene_ranks.tsv.gz']}

summary = {}

for subdir,patterns in subdirs_and_patterns.items():
    
    if subdir in listdir(main_dir):
    
        for file in listdir(f'{main_dir}/{subdir}'):
            
            if sum([file.endswith(p) for p in patterns]) == 0:
                
                continue
            
            try:
                
                # Extract sample id
                
                sample_id = file

                for p in patterns:
                    
                    sample_id = sample_id.replace(p, '')
                
                if sample_id not in summary.keys():
                    
                    summary[sample_id] = {'with_hpo' : {}, 'without_hpo' : {}}
                
                # Extract data type
                
                main_type = 'with_hpo' if 'with_hpo' in file else 'without_hpo'
                
                sub_type = ('aberrant_splicing' if 'aberrant_splicing' in file else
                            'aberrant_expression' if 'aberrant_expression' in file else
                            'allelic_imbalance' if 'allelic_imbalance' in file or 'ase' in file else
                            'integrated_data' if 'integrated_data' in file else
                            'integrated_data')
                
                pval_column = rank_type_pval_columns[sub_type]
                
                effect_column = rank_type_effect_columns[sub_type]
                
                # Extract data
                
                data = pd.read_csv(f'{main_dir}/{subdir}/{file}', sep='\t')
                
                if target_gene in data['gene_symbol'].values:
                
                    data.loc[:, '-pval'] = - data[pval_column].values
                    
                    if main_type == 'with_hpo':
                    
                        data.sort_values(['hpo_overlap', '-pval', effect_column], ascending=False, inplace=True)
                        
                    else:
                        
                        data.sort_values(['-pval', effect_column], ascending=False, inplace=True)
                    
                    data = data.reset_index(drop=True)
                    
                    rank = data.index[data['gene_symbol'] == target_gene][0] + 1
                    
                    summary[sample_id][main_type][sub_type] = rank
                
                else:
                    
                    summary[sample_id][main_type][sub_type] = np.nan
            
            except:
                
                pass

for sample_id,sample_summary in summary.items():

    sample_summary = pd.DataFrame(sample_summary)
    sample_summary.index.name = 'analysis'
    sample_summary = sample_summary.reset_index(drop=False)
    
    sample_summary.to_csv(f'{sample_id}_{target_gene}_ranks.tsv', sep='\t', index=False, header=True)
