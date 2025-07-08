#!/usr/bin/env python3

### ---------------------------------------- ###

def parse_argv():

    print(f'args = {argv}')

    # Sample ID
    
    sample_id = argv[argv.index('--sample') + 1]
    
    # Leafcutter files
    
    intron_pvals_path = argv[argv.index('--intron_pvals') + 1]
    intron_pvals = pd.read_csv(intron_pvals_path, sep='\t')
    
    intron_effect_path = argv[argv.index('--intron_effect') + 1]
    intron_effect = pd.read_csv(intron_effect_path, sep='\t')
    
    cluster_pvals_path = argv[argv.index('--cluster_pvals') + 1]
    cluster_pvals = pd.read_csv(cluster_pvals_path, sep='\t')
    
    # GTF annotation
    
    gtf_path = argv[argv.index('--gtf') + 1]
    gtf = load_gtf(gtf_path)
    
    return sample_id, intron_pvals, intron_effect, cluster_pvals, gtf

### ---------------------------------------- ###

def load_gtf(path, desired_biotypes=[], desired_chromosomes=[]):
    
    # Load GTF
    
    gtf_data = pd.read_csv(path, sep='\t', header=None, comment='#', dtype=str)
    gtf_data.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    ### Only keep genes and exons

    gtf_data = gtf_data.loc[gtf_data.feature.isin(['gene', 'exon']), ['seqname', 'start', 'end', 'strand', 'attribute']]
    
    # Get biotype and gene id

    biotypes, gene_ids, gene_symbols, transcript_ids, exon_ids = [], [], [], [], []
    for _,row in gtf_data.iterrows():
        
        info = row.values[-1]
        
        biotype = re.findall('gene_biotype "\w+";', info)[0]
        biotype = biotype.replace('gene_biotype ', '').replace(';', '').replace('"', '')
        
        biotypes.append(biotype)
        
        gene = re.findall('gene_id "\w+";', info)[0]
        gene = gene.replace('gene_id ', '').replace(';', '').replace('"', '')
        
        gene_ids.append(gene)
        
        if 'gene_name' in info:
            
            gene = info[info.index('gene_name "') + len('gene_name "'):]
            gene = gene[:gene.index('"')]
        
        else:
            
            gene = ''
        
        gene_symbols.append(gene)
        
        if 'transcript_id' in info:
            
            transcript = info[info.index('transcript_id "') + len('transcript_id "'):]
            transcript = transcript[:transcript.index('"')]
        
        else:
            
            transcript = ''
        
        transcript_ids.append(transcript)
        
        if 'exon_id' in info:
            
            exon = info[info.index('exon_id "') + len('exon_id "'):]
            exon = exon[:exon.index('"')]
        
        else:
            
            exon = ''
        
        exon_ids.append(exon)

    gtf_data['biotype'] = biotypes
    gtf_data['gene_id'] = gene_ids
    gtf_data['gene_symbol'] = gene_symbols
    gtf_data['transcript_id'] = transcript_ids
    gtf_data['exon_id'] = exon_ids
    
    # Filter based on biotype
    
    if len(desired_biotypes):

        gtf_data = gtf_data.loc[gtf_data.biotype.isin(desired_biotypes),]

    # Filter for desired chromosomes

    if len(desired_chromosomes):

        gtf_data = gtf_data.loc[gtf_data.seqname.isin(desired_chromosomes),]
    
    # Remove genes without gene_symbol
    
    #gtf_data = gtf_data.loc[gtf_data['gene_symbol'] != '',]
    
    # Fix dtypes
    
    gtf_data[['start', 'end']] = gtf_data[['start', 'end']].astype(int)
    
    return gtf_data

### ---------------------------------------- ###

def merge_data(sid, annot, i_pvals, i_effects, cl_pvals):
    
    # Filter data for sample of interest

    i_pvals = i_pvals[sid]
    i_pvals.name = 'intron_pval'
    
    i_effects = i_effects[sid]
    i_effects.name = 'intron_effect_size'
    
    cl_pvals = cl_pvals[sid]
    cl_pvals.name = 'cluster_pval'
    
    # Benjamini-Hochberg correction of cluster p values
    
    cl_pvals[:] = fdrcorrection(cl_pvals.values, alpha=0.05, is_sorted=False)[1]
    
    # Create intron to cluster series

    intron_to_cluster = pd.Series({intron : intron.split(':')[-1] for intron in i_pvals.index})
    intron_to_cluster.name = 'cluster'
    
    # Merge
    
    merged = pd.merge(i_pvals, i_effects, how='inner', left_index=True, right_index=True)
    merged = pd.merge(merged, intron_to_cluster, how='inner', left_index=True, right_index=True)
    merged = pd.merge(merged, cl_pvals, how='inner', left_on='cluster', right_index=True)
    
    # Reset index to add intron ID
    
    merged.index.name = 'intron_id'
    merged.reset_index(drop=False, inplace=True)
    
    # Add intron coordinates
    
    intron_coords = [(i.split(':')[0],
                      int(i.split(':')[1]),
                      int(i.split(':')[2]))
                     for i in merged['intron_id'].values]
    
    intron_coords = pd.DataFrame(intron_coords,
                                 columns=['seqname', 'start', 'end'],
                                 index=merged['intron_id'].values)
    
    merged = pd.merge(intron_coords, merged, how='inner', left_index=True, right_on='intron_id')
    
    # Find associated genes based on intersection of intron and exon coordinates (can be more than one gene)
    # N.B. Exon starts and ends have to be considered separately given the nature of introns
    # (intron start and end correspond to end of an exon and start of another, respectively)
    
    exon_left = pd.merge(annot, merged, how='inner', left_on=['seqname', 'end'], right_on=['seqname', 'start'], suffixes=['_exon', '_intron'])
    
    exon_right = pd.merge(annot, merged, how='inner', left_on=['seqname', 'start'], right_on=['seqname', 'end'], suffixes=['_exon', '_intron'])
    
    merged = pd.concat([exon_left, exon_right], axis=0)
    
    merged.drop_duplicates(inplace=True)
    
    merged.drop(['attribute', 'intron_id'], axis='columns', inplace=True)
    
    return merged

### ---------------------------------------- ###

def summarize_leafcutter(dt, p_thr=0.05):
    
    # Clusters stats
    
    dt_sub = dt[['gene_id', 'cluster', 'cluster_pval']].drop_duplicates()
    
    gene_clusters = dt_sub.groupby('gene_id')['cluster'].size()

    gene_clusters_significant = dt_sub.loc[dt_sub['cluster_pval'] < 0.05,].groupby('gene_id')['cluster'].size()
    
    gene_clusters_best_p = dt_sub.groupby('gene_id')['cluster_pval'].min()
    
    # Intron stats
    
    dt_sub = dt[['gene_id', 'seqname', 'start_intron', 'end_intron', 'intron_pval', 'intron_effect_size']].drop_duplicates()
    
    dt_sub.loc[:, 'intron_id'] = [':'.join([c, str(s), str(e)]) for _,(c,s,e) in dt_sub[['seqname', 'start_intron', 'end_intron']].iterrows()]
    
    dt_sub.loc[:, 'intron_effect_size'] = dt_sub.loc[:, 'intron_effect_size'].abs()
    
    gene_introns = dt_sub.groupby('gene_id')['intron_id'].size()
    
    gene_introns_significant = dt_sub.loc[dt_sub['intron_pval'] < 0.05,].groupby('gene_id')['intron_id'].size()
    
    gene_introns_best_p = dt_sub.groupby('gene_id')['intron_pval'].min()
    
    gene_introns_best_effect = dt_sub.groupby('gene_id')['intron_effect_size'].max()
    
    gene_introns_median_effect = dt_sub.groupby('gene_id')['intron_effect_size'].median()
    
    # Merge
    
    gene_data = dt[['gene_id', 'gene_symbol', 'biotype', 'seqname', 'strand']].drop_duplicates()
    gene_data.index = gene_data['gene_id'].values
    
    for db in [gene_clusters, gene_clusters_significant, gene_clusters_best_p, gene_introns, gene_introns_significant, gene_introns_best_p, gene_introns_best_effect, gene_introns_median_effect]:
    
        gene_data = pd.merge(gene_data, db, how='outer', left_index=True, right_index=True)
    
    gene_data.fillna(0, inplace=True)
    
    # Finalize format
    
    gene_data.columns = ['gene_id', 'gene_symbol', 'gene_biotype', 'contig', 'strand',
                         'n_clusters', 'n_clusters_p<0.05', 'n_clusters_min_pval',
                         'n_introns', 'n_introns_p<0.05', 'n_introns_min_pval', 'intron_best_abs_effect', 'intron_median_abs_effect']
    
    return gene_data

### ------------------MAIN------------------ ###

import pandas as pd
import re

from statsmodels.stats.multitest import fdrcorrection
from sys import argv

### Load data

sample_id, intron_pvals, intron_effect, cluster_pvals, gtf = parse_argv()

### Merge data

all_data = merge_data(sample_id, gtf, intron_pvals, intron_effect, cluster_pvals)

all_data.to_csv(sample_id + '_splicing_stats.tsv.gz', sep='\t', index=False, header=True)

### Summarize at gene level

gene_data = summarize_leafcutter(all_data)

gene_data.to_csv(sample_id + '_abberrant_splicing.tsv.gz', sep='\t', index=False, header=True)
