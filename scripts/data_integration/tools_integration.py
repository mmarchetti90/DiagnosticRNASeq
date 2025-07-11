#!/usr/bin/env python3

"""
This script integrates data from ABERRANT_SPLICING, ABERRANT_EXPRESSION, and ALLELIC_IMBALANCE
workflows to prioritize genes
"""

### ---------------------------------------- ###

def parse_args():
    
    # Aberrant splicing data
    
    if '--aberrant_splicing' in argv:
        
        aberrant_splicing_path = argv[argv.index('--aberrant_splicing') + 1]
        aberrant_splicing_data = pd.read_csv(aberrant_splicing_path, sep='\t', header=0)
    
    else:
        
        columns = [
            'gene_id', 'gene_symbol', 'gene_biotype', 'contig', 'strand',
            'n_clusters', 'n_clusters_p<0.05', 'n_clusters_min_pval', 'n_introns', 'n_introns_p<0.05',
            'n_introns_min_pval', 'intron_best_abs_effect', 'intron_median_abs_effect', 'intron_median_abs_effect_scaled'
        ]
        
        aberrant_splicing_data = pd.DataFrame({}, columns=columns)
    
    # Aberrant expression data
    
    if '--aberrant_expression' in argv:
        
        aberrant_expression_path = argv[argv.index('--aberrant_expression') + 1]
        aberrant_expression_data = pd.read_csv(aberrant_expression_path, sep='\t', header=0)
    
    else:
        
        columns = [
            'gene_id', 'gene_symbol', 'gene_biotype', 'contig', 'start',
            'end', 'strand', 'pval', 'padj', 'log2FC', 'log2FC_scaled'
        ]
        
        aberrant_expression_data = pd.DataFrame({}, columns=columns)
    
    # Allelic imbalance data
    
    if '--allelic_imbalance' in argv:
        
        allelic_imbalance_path = argv[argv.index('--allelic_imbalance') + 1]
        allelic_imbalance_data = pd.read_csv(allelic_imbalance_path, sep='\t', header=0)
    
    else:
        
        columns = [
            'gene_id', 'gene_symbol', 'gene_biotype', 'contig', 'start',
            'end', 'strand', 'n_snps', 'best_ase_score', 'median_ase_score',
            'median_ase_score_scaled', 'best_p', 'compound_p'
        ]
        
        allelic_imbalance_data = pd.DataFrame({}, columns=columns)
    
    return aberrant_splicing_data, aberrant_expression_data, allelic_imbalance_data

### ---------------------------------------- ###

def scale_data(vals):
    
    if len(vals) > 0:
        
        abs_max_val = np.max(np.abs(vals))
        
        abs_max_val = abs_max_val if abs_max_val > 0 else 1
        
        vals /= abs_max_val
        
    return vals

### ---------------------------------------- ###

def merge_data(a_s, a_e, a_i, cols, metric='final_pval'):
    
    a_s_col, a_e_col, a_i_col = cols
    
    common_cols = ['gene_id', 'gene_symbol', 'gene_biotype', 'contig', 'strand']

    all_dt = pd.merge(a_s[common_cols + [a_s_col]],
                      a_e[common_cols + [a_e_col]],
                      on=common_cols,
                      how='outer')

    all_dt = pd.merge(all_dt,
                      a_i[common_cols + [a_i_col]],
                      on=common_cols,
                      how='outer')

    all_dt[cols] = all_dt[cols].astype('float64')

    all_dt[cols] = all_dt[cols].fillna(1)

    all_dt[metric] = np.power(np.prod(all_dt[cols].abs(), axis=1), 1 / len(cols))

    all_dt.drop(cols, axis=1, inplace=True)
    
    return all_dt

### ------------------MAIN------------------ ###

try:
    
    import numpy as np
    import pandas as pd
    
    from sys import argv
    
except:
    
    print("One or more dependencies are not installed.\nAlso, make sure your terminal has been activated.")
    exit()

### Parse args

aberrant_splicing_data, aberrant_expression_data, allelic_imbalance_data = parse_args()

### Scale effect sizes

# Aberrant splicing

if aberrant_splicing_data.shape[0] > 0:
    
    aberrant_splicing_data.loc[:, 'intron_median_abs_effect_scaled'] = scale_data(aberrant_splicing_data['intron_median_abs_effect'].values.copy())

else:
    
    aberrant_splicing_data.loc[:, 'intron_median_abs_effect_scaled'] = aberrant_splicing_data['intron_median_abs_effect']

# Aberrant expression

if aberrant_expression_data.shape[0] > 0:
    
    aberrant_expression_data.loc[:, 'log2FC_scaled'] = scale_data(aberrant_expression_data['log2FC'].values.copy())

else:
    
    aberrant_expression_data.loc[:, 'log2FC_scaled'] = aberrant_expression_data['log2FC']

# Allelic imbalance

if allelic_imbalance_data.shape[0] > 0:
    
    allelic_imbalance_data.loc[:, 'median_ase_score_scaled'] = scale_data(allelic_imbalance_data['median_ase_score'].values.copy())

else:
    
    allelic_imbalance_data.loc[:, 'median_ase_score_scaled'] = allelic_imbalance_data['median_ase_score']

### Merge p values

pvals = merge_data(aberrant_splicing_data,
                   aberrant_expression_data,
                   allelic_imbalance_data,
                   ['n_clusters_min_pval', 'padj', 'compound_p'],
                   'final_pval')

### Merge effect size

eff_size = merge_data(aberrant_splicing_data,
                      aberrant_expression_data,
                      allelic_imbalance_data,
                      ['intron_median_abs_effect', 'log2FC', 'median_ase_score'],
                      'final_effect')

### Merge scaled effect size (used to break ties)

eff_size_scaled = merge_data(aberrant_splicing_data,
                             aberrant_expression_data,
                             allelic_imbalance_data,
                             ['intron_median_abs_effect_scaled', 'log2FC_scaled', 'median_ase_score_scaled'],
                             'final_effect_scaled')

### Merge

ranking = pd.merge(pvals,
                   eff_size,
                   on=['gene_id', 'gene_symbol', 'gene_biotype', 'contig', 'strand'],
                   how='inner')

ranking = pd.merge(ranking,
                   eff_size_scaled,
                   on=['gene_id', 'gene_symbol', 'gene_biotype', 'contig', 'strand'],
                   how='inner')

ranking.sort_values(by=['final_pval', 'final_effect_scaled'],
                    ascending=True,
                    inplace=True)

ranking.to_csv('gene_ranks.tsv.gz', sep='\t', index=False, header=True)
