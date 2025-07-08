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
            'n_introns_min_pval', 'intron_best_abs_effect', 'intron_median_abs_effect'
        ]
        
        aberrant_splicing_data = pd.DataFrame({}, columns=columns)
    
    # Aberrant expression data
    
    if '--aberrant_expression' in argv:
        
        aberrant_expression_path = argv[argv.index('--aberrant_expression') + 1]
        aberrant_expression_data = pd.read_csv(aberrant_expression_path, sep='\t', header=0)
    
    else:
        
        columns = [
            'gene_id', 'gene_symbol', 'gene_biotype', 'contig', 'start',
            'end', 'strand', 'pval', 'padj', 'log2FC'
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
            'best_p', 'compound_p'
        ]
        
        allelic_imbalance_data = pd.DataFrame({}, columns=columns)
    
    return aberrant_splicing_data, aberrant_expression_data, allelic_imbalance_data

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

    all_dt[cols] = all_dt[cols].fillna(1)

    all_dt[metric] = np.prod(all_dt[cols].abs(), axis=1)

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

### Merge p values

pvals = merge_data(aberrant_splicing_data,
                   aberrant_expression_data,
                   allelic_imbalance_data,
                   ['n_clusters_min_pval', 'padj', 'compound_p'],
                   'final_pval')

### Merge effect size (used to break ties)

eff_size = merge_data(aberrant_splicing_data,
                      aberrant_expression_data,
                      allelic_imbalance_data,
                      ['intron_median_abs_effect', 'log2FC', 'median_ase_score'],
                      'final_effect')

### Merge

ranking = pd.merge(pvals,
                   eff_size,
                   on=['gene_id', 'gene_symbol', 'gene_biotype', 'contig', 'strand'],
                   how='inner')

ranking.sort_values(by=['final_pval', 'final_effect'],
                    ascending=True,
                    inplace=True)

ranking.to_csv('gene_ranks.tsv.gz', sep='\t', index=False, header=True)
