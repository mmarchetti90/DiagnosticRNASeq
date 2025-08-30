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
    
    # Ranking type (one of 'integrated_data', 'aberrant_splicing', 'aberrant_expression', 'allelic_imbalance')
    
    genes_rank_type = argv[argv.index('--rank_type') + 1]
    
    return hpo_id_to_name, hpo_hierarchy, target_hpo_terms, genes2hpo, genes_rank, genes_rank_type

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

def kneedle(vector, sort_vector=True):
    
    """
    Kneedle to find threshold cutoff.
    """
    
    if sort_vector:
        
        vector = np.sort(vector)[::-1]
    
    # Find gradient and intercept
    x0, x1 = 0, len(vector)
    y0, y1 = max(vector), min(vector)
    gradient = (y1 - y0) / (x1 - x0)
    intercept = y0
    
    # Compute difference vector
    difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(vector)]
    
    # Find max of difference_vector and define cutoff
    cutoff_index = difference_vector.index(max(difference_vector))
    cutoff_value = vector[cutoff_index]
    
    return cutoff_index, cutoff_value

### ---------------------------------------- ###

def hpo_distances(hierarchy, targets, hits, max_distance=3):
    
    distances = []
    
    for t in targets:
        
        # Processing hits that are parents of t
        
        t_parents = get_parent_terms([t], hierarchy)
        
        hits_sub_1 = [h for h in hits if h in t_parents]
        
        for hs in hits_sub_1:
            
            dist = get_parent_child_distance(hs, t, hierarchy)
            
            if dist <= max_distance:
            
                distances.append([hs, t, dist, 'y'])
        
        # Processing hits that are children of t
        
        t_childrens = get_children_terms([t], hierarchy)
        
        hits_sub_2 = [h for h in hits if h in t_childrens and h != t]
        
        for hs in hits_sub_2:
            
            dist = get_parent_child_distance(t, hs, hierarchy)
            
            if dist <= max_distance:
            
                distances.append([hs, t, dist, 'y'])
        
        # Processing terms from other branches
        
        for hs in hits:
            
            if hs in (hits_sub_1 + hits_sub_2):
                
                continue
            
            _, dist = find_closest_common_ancestor(hs, t, hierarchy)
            
            if dist <= max_distance:
                
                distances.append([hs, t, dist, 'n'])
    
    distances = pd.DataFrame(distances, columns=['hit_hpo', 'target_hpo', 'distance', 'same_branch'])
    
    # Cleanup so that each hit HPO is associated with only the closest HPO term
    # N.B. target terms in the same branch as the hit are prioritized
    
    distances_cleaned = []
    
    for h in np.unique(distances['hit_hpo'].values):
        
        distances_sub = distances.loc[distances['hit_hpo'] == h,].copy()
    
        distances_sub = distances_sub.replace({'same_branch' : {'y' : '0', 'n' : '1'}})
    
        distances_sub = distances_sub.sort_values(by=['same_branch', 'distance'], ascending=True)
        
        distances_cleaned.append(distances_sub.iloc[0,].values)
    
    distances_cleaned = pd.DataFrame(distances_cleaned, columns=['hit_hpo', 'target_hpo', 'distance', 'same_branch'])
    
    distances_cleaned = distances_cleaned.replace({'same_branch' : {'0' : 'y', '1' : 'n'}})
    
    return distances_cleaned

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

def find_closest_common_ancestor(element_1, element_2, parent_child):
    
    # Get ancestors in common
    
    parents_1 = get_parent_terms([element_1], parent_child)
    
    parents_2 = get_parent_terms([element_2], parent_child)
    
    common_ancestors = [p1 for p1 in parents_1 if p1 in parents_2 and p1 != '']
    
    # Find shortest distance
    
    elements_distance = []
    
    for ca in common_ancestors:
        
        distance_1 = get_parent_child_distance(ca, element_1, parent_child)
        
        distance_2 = get_parent_child_distance(ca, element_2, parent_child)
        
        elements_distance.append((ca, distance_1 + distance_2 + 1))
    
    elements_distance.sort(key=lambda ed: ed[1])
    
    # Closest ancestor
    
    closest_ancestor = elements_distance[0]
    
    return closest_ancestor

### ---------------------------------------- ###

import pandas as pd
import numpy as np

from sys import argv

### Parse args

hpo_id_to_name, hpo_hierarchy, target_hpo_terms, genes2hpo, genes_rank, genes_rank_type = parse_args()

### Processing

if genes_rank.shape[0] == 0:
    
    genes_rank.loc[:, ['hpo_hits', 'hpo_hits_distance', 'hpo_overlap']] = []
    
    genes_rank.to_csv(f'gene_ranks_with_hpo-{genes_rank_type}.tsv.gz', sep='\t', index=False, header=True)

else:
    
    ### Cutoff final_pval using a kneedle (but keep N genes max)
    
    rank_type_pval_columns = {'integrated_data' : 'final_pval',
                              'aberrant_splicing' : 'clusters_min_pval',
                              'aberrant_expression' : 'padj',
                              'allelic_imbalance' : 'compound_p'}
    
    pval_column = rank_type_pval_columns[genes_rank_type]
    
    genes_rank = genes_rank.sort_values(by=pval_column)
    
    minus_log_pval = - np.log10(genes_rank[pval_column].values + 1e-6)
    
    _, thr = kneedle(minus_log_pval.copy(), True)
    
    genes_rank = genes_rank.loc[minus_log_pval >= thr,]
    
    N = 3000
    
    genes_rank = genes_rank.iloc[:N,]
    
    ### Init score matrix for HPO terms
    
    # N.B. Genes may be associated with an HPO term with a certain hierarchical distance from a targer
    # one. Only terms with distance < 3 nodes will be considered
    
    max_branch_length = 3

    hpo_hits = genes2hpo.loc[genes2hpo['gene_symbol'].isin(genes_rank['gene_symbol']), 'hpo_id'].values
    
    hpo_hits = np.unique(hpo_hits)

    hit_hpo_distances_to_target = hpo_distances(hpo_hierarchy, target_hpo_terms, hpo_hits, max_branch_length)
    
    ### Score genes based on HPO overlap
    
    genes_rank.loc[:, 'hpo_hits'] = np.repeat('', genes_rank.shape[0])
    genes_rank.loc[:, 'hpo_hits_distance'] = np.repeat('', genes_rank.shape[0])
    genes_rank.loc[:, 'hpo_hits_same_branch'] = np.repeat('', genes_rank.shape[0])
    genes_rank.loc[:, 'hpo_overlap'] = np.repeat(0., genes_rank.shape[0])
    
    genes_rank = genes_rank.reset_index(drop=True)
    
    for idx,gene_info in genes_rank.iterrows():
        
        if gene_info['gene_symbol'] not in genes2hpo['gene_symbol'].values:
            
            continue
        
        # Get gene's HPO
        
        gene_hpo = genes2hpo.loc[genes2hpo['gene_symbol'] == gene_info['gene_symbol'], 'hpo_id'].values
        
        # Find hits between the gene's HPO terms and the target HPO terms
        
        gene_hpo_hits = hit_hpo_distances_to_target.loc[hit_hpo_distances_to_target['hit_hpo'].isin(gene_hpo),]
        
        # Cleanup
        
        gene_hpo_hits = gene_hpo_hits.replace({'same_branch' : {'y' : '0', 'n' : '1'}})
        
        gene_hpo_hits = gene_hpo_hits.sort_values(by=['target_hpo', 'same_branch', 'distance'], ascending=True)
        
        gene_hpo_hits = gene_hpo_hits.drop_duplicates('target_hpo', keep='first')
        
        gene_hpo_hits = gene_hpo_hits.replace({'same_branch' : {'0' : 'y', '1' : 'n'}})
        
        # Annotated
        
        genes_rank.loc[idx, 'hpo_hits'] = ';'.join(gene_hpo_hits['target_hpo'].values)
        
        genes_rank.loc[idx, 'hpo_hits_distance'] = ';'.join([str(d) for d in gene_hpo_hits['distance'].values])
        
        genes_rank.loc[idx, 'hpo_hits_same_branch'] = ';'.join(gene_hpo_hits['same_branch'].values)
        
        genes_rank.loc[idx, 'hpo_overlap'] = gene_hpo_hits.shape[0] / len(target_hpo_terms)
    
    ### Sort
    
    genes_rank.loc[:, '-log10(pval)'] = - np.log10(genes_rank[pval_column].values + 1e-6)
    
    genes_rank.sort_values(['hpo_overlap', '-log10(pval)'], ascending=False, inplace=True)
    
    genes_rank = genes_rank.drop('-log10(pval)', axis=1)
    
    ### Save to file
    
    genes_rank.to_csv(f'gene_ranks_with_hpo-{genes_rank_type}.tsv.gz', sep='\t', index=False, header=True)
