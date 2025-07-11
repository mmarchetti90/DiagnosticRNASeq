#!/usr/bin/env python3

### ---------------------------------------- ###

def parse_argv():

    print(f'args = {argv}')

    # Sample ID
    
    sample_id = argv[argv.index('--sample') + 1]
    
    # Outrider data
    
    outrider_path = argv[argv.index('--outrider_data') + 1]
    outrider_data = pd.read_csv(outrider_path, sep='\t')
    
    # GTF annotation
    
    gtf_path = argv[argv.index('--gtf') + 1]
    gtf = load_gtf(gtf_path)
    
    # Only keep gene info
    
    gene_starts = gtf.groupby(['gene_id', 'gene_symbol', 'biotype', 'seqname', 'strand'])['start'].min().reset_index(drop=False)
    gene_ends = gtf.groupby(['gene_id', 'gene_symbol', 'biotype', 'seqname', 'strand'])['end'].max().reset_index(drop=False)
    
    gtf = pd.merge(gene_starts, gene_ends, how='inner', on=['gene_id', 'gene_symbol', 'biotype', 'seqname', 'strand'])
    
    gtf = gtf[['gene_id', 'gene_symbol', 'biotype', 'seqname', 'start', 'end', 'strand']]
    gtf.columns = ['gene_id', 'gene_symbol', 'gene_biotype', 'contig', 'start', 'end', 'strand']
    
    return sample_id, outrider_data, gtf

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

### ------------------MAIN------------------ ###

import pandas as pd
import re

from sys import argv

### Load data

sample_id, outrider_data, gtf = parse_argv()

### Parse data

sample_data = outrider_data.loc[(outrider_data['sampleID'] == sample_id) &
                                (outrider_data['aberrant']),
                                ['geneID', 'pValue', 'padjust', 'l2fc']]

sample_data.columns = ['gene_id', 'pval', 'padj', 'log2FC']

sample_data = pd.merge(gtf, sample_data, on='gene_id', how='inner')

sample_data.to_csv(sample_id + '_aberrant_expression.tsv.gz', sep='\t', index=False, header=True)
