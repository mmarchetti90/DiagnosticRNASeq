#!/usr/bin/env python3

"""
Interpolating genes with target HPO terms
"""

### ---------------------------------------- ###

def parse_args():
    
    # HPO ontology
    
    obo_path = argv[argv.index('--obo') + 1]
    
    hpo_id_to_name, hpo_hierarchy = parse_obo(obo_path)
    
    # Target HPO terms
    
    target_hpo_terms_path = argv[argv.index('--hpo') + 1]
    
    target_hpo_terms = [term for term in open(target_hpo_terms_path, 'r').read().split('\n') if len(term)]
    
    # Genes to HPO
    
    genes2hpo_path = argv[argv.index('--genes2hpo') + 1]
    
    genes2hpo = pd.read_csv(genes2hpo_path, sep='\t')[['gene_symbol', 'hpo_id']]
    
    # Genes ranking
    
    genes_rank_path = argv[argv.index('--genes_rank') + 1]
    
    genes_rank = pd.read_csv(genes_rank_path, sep='\t')
    
    return hpo_id_to_name, hpo_hierarchy, target_hpo_terms, genes2hpo, genes_rank

### ---------------------------------------- ###

def parse_obo(obo_file):

    # This creates a tree structure for deep search
    
    # Read obo file and create hierarchy
    
    id_to_name = {}
    hpo_hierarchy = []
    for term in open(obo_file).read().split('\n\n'):
        
        if term.startswith('[Term]'):
            
            # Parse term info
            
            term = term.split('\n')
            
            term_ids = [t.replace('alt_id: ', '').replace('id: ', '') for t in term if t.startswith('id: ') or t.startswith('alt_id: ')]
            
            term_name = [t.replace('name: ', '') for t in term if t.startswith('name: ')]
            term_name = term_name[0] if len(term_name) else ''
            
            parents = [t[:t.index(' !')].replace('is_a: ', '') for t in term if t.startswith('is_a: ')]
            
            # Update id_to_name
            
            for t in term_ids:
                
                id_to_name[t] = term_name
            
            # Update hpo_hierarchy
            
            if not len(parents):
                
                for t in term_ids:
                    
                    hpo_hierarchy.append(('', t))
            
            else:
            
                for p in parents:
                    
                    for t in term_ids:
                        
                        hpo_hierarchy.append((p, t))

    hpo_hierarchy = pd.DataFrame(hpo_hierarchy, columns=['parent', 'child'])
    
    return id_to_name, hpo_hierarchy

### ---------------------------------------- ###

def hpo_overlap(hierarchy, targets, hits, max_dist):
    
    target_hits, overlap = {}, 0
    
    for h in hits:
        
        if h in targets:
            
            target_hits[h] = 0
        
        else:
            
            h_parents = get_parent_terms([h], hierarchy)
            
            h_children = get_children_terms([h], hierarchy)
            
            min_dist = ('', max_dist + 1)
            
            for t in targets:
                
                if t in h_parents:
                    
                    dist = get_parent_child_distance(t, h, hierarchy)
                    
                    if dist < min_dist[1]:
                        
                        min_dist = (t, dist)
                
                elif t in h_children:

                    dist = get_parent_child_distance(h, t, hierarchy)
                    
                    if dist < min_dist[1]:
                        
                        min_dist = (t, dist)
                    
                else:
                    
                    continue
            
            if min_dist[1] <= max_dist:
                
                if min_dist[0] not in target_hits.keys():
                    
                    target_hits[min_dist[0]] = min_dist[1]
                
                else:
                    
                    if min_dist[1] < target_hits[min_dist[0]]:
                        
                        target_hits[min_dist[0]] = min_dist[1]
            
    overlap = len(target_hits) / len(targets)
        
    return target_hits, overlap

### ---------------------------------------- ###

def get_children_terms(c, parent_child):
    
    children = set(parent_child.loc[parent_child.parent.isin(c), 'child'].to_list())
    
    children.update(c)

    if children == c:
        
        return c
    
    else:
        
        return get_children_terms(children, parent_child)

### ---------------------------------------- ###

def get_parent_terms(p, parent_child):
    
    parents = set(parent_child.loc[parent_child.child.isin(p), 'parent'].to_list())
    
    parents.update(p)

    if parents == p:
        
        return p
    
    else:
        
        return get_parent_terms(parents, parent_child)

### ---------------------------------------- ###

def get_parent_child_distance(p, c, parent_child):
    
    children = [p]
    
    n = 0
    
    while c not in children:
        
        n += 1
        
        children = set(parent_child.loc[parent_child.parent.isin(children), 'child'].to_list())
    
    return n

### ---------------------------------------- ###

import pandas as pd
import numpy as np

from sys import argv

### Parse args

hpo_id_to_name, hpo_hierarchy, target_hpo_terms, genes2hpo, genes_rank = parse_args()

### Only consider the top 1000 genes

N = 1000

genes_rank = genes_rank.iloc[:N,]

### Score genes based on HPO overlap

# N.B. Genes may be associated with an HPO term with a certain hierarchical distance from a targer
# one. Only terms with distance < 3 nodes will be considered

max_branch_length = 3

genes_rank.loc[:, 'hpo_hits'] = np.repeat('', genes_rank.shape[0])
genes_rank.loc[:, 'hpo_hits_distance'] = np.repeat('', genes_rank.shape[0])
genes_rank.loc[:, 'hpo_overlap'] = np.repeat(0, genes_rank.shape[0])

genes_rank = genes_rank.reset_index(drop=True)

for idx,gene_info in genes_rank.iterrows():
    
    if gene_info['gene_symbol'] not in genes2hpo['gene_symbol'].values:
        
        continue
    
    gene_hpo = genes2hpo.loc[genes2hpo['gene_symbol'] == gene_info['gene_symbol'], 'hpo_id'].values
    
    target_hits, overlap = hpo_overlap(hpo_hierarchy, target_hpo_terms, gene_hpo, max_branch_length)
    
    genes_rank.loc[idx, 'hpo_hits'] = ';'.join(target_hits.keys())
    
    genes_rank.loc[idx, 'hpo_hits_distance'] = ';'.join([str(v) for v in target_hits.values()])
    
    genes_rank.loc[idx, 'hpo_overlap'] = overlap

genes_rank.sort_values('hpo_overlap', ascending=False, inplace=True)

genes_rank.to_csv('gene_ranks_with_hpo.tsv.gz', sep='\t', index=False, header=True)
