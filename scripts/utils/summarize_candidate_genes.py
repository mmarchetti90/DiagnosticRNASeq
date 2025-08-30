#!/usr/bin/env python3

"""
Script to conveniently interpolate results with a list of target genes

Parameters
----------
results_dir : str
    Path to results directory (can be .)
target_genes : str
    Path to file with one target gene per line
"""

### ---------------------------------------- ###

### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt
from os import listdir
from sys import argv

### Parse args

main_dir = argv[argv.index('--results_dir') + 1]

target_genes_file = argv[argv.index('--target_genes') + 1]

target_genes = [gene for gene in open(target_genes_file, 'r').read().split('\n') if len(gene)]

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
                    
                    summary[sample_id] = {gene : {'with_hpo' : {}, 'without_hpo' : {}} for gene in target_genes}
                
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
                
                for target_gene in target_genes:
                
                    if target_gene in data['gene_symbol'].values:
                    
                        data.loc[:, '-pval'] = - data[pval_column].values
                        
                        if main_type == 'with_hpo':
                        
                            data.sort_values(['hpo_overlap', '-pval', effect_column], ascending=False, inplace=True)
                            
                        else:
                            
                            data.sort_values(['-pval', effect_column], ascending=False, inplace=True)
                        
                        data = data.reset_index(drop=True)
                        
                        rank = data.index[data['gene_symbol'] == target_gene][0] + 1
                        
                        summary[sample_id][target_gene][main_type][sub_type] = rank
                    
                    else:
                        
                        summary[sample_id][target_gene][main_type][sub_type] = np.nan
            
            except:
                
                pass

### Structure as heatmap

max_rank = 100

for sample_id,sample_summary in summary.items():
    
    heatmap_idx = target_genes
    heatmap_col = [f'{sub_type}_{main_type}' for main_type in ['with_hpo', 'without_hpo'] for sub_type in rank_type_pval_columns.keys()]
    
    heatmap_data = pd.DataFrame(np.zeros((len(heatmap_idx), len(heatmap_col))) - 1,
                                index=heatmap_idx,
                                columns=heatmap_col)
    
    for row in heatmap_idx:
        
        for col in heatmap_col:
            
            main_type = 'with_hpo' if 'with_hpo' in col else 'without_hpo'
            sub_type = col.replace('_with_hpo', '').replace('_without_hpo', '')
            
            rank = sample_summary[row][main_type][sub_type]
            
            heatmap_data.loc[row, col] = rank
    
    #heatmap_data = heatmap_data.fillna(-1)
    
    heatmap_data.to_csv(f'{sample_id}_ranks_summary.tsv', sep='\t', index=True, header=True)
    
    # Plot data
    
    plt.figure(figsize=(12, 12))
    sns.heatmap(heatmap_data, cmap='Blues_r', vmin=1, vmax=100, linecolor='white', linewidths=1)
    plt.xticks(ticks=np.arange(0.5, len(heatmap_col) + 0.5, 1),
               labels=[col.replace('_with_hpo', '\nwith_hpo').replace('_without_hpo', '\nwithout_hpo').replace('_', ' ') for col in heatmap_col],
               fontweight='bold')
    plt.yticks(fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{sample_id}_ranks_summary.png', dpi=300)
    plt.close()
    

