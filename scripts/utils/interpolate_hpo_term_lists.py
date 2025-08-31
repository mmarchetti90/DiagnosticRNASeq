#!/usr/bin/env python3

"""
This script interpolates two sets of HPO terms and calculates their distance

Arguments
---------
list1
    Path to txt file with list of HPO terms
    1 HPO term per row
list2
    Path to txt file with list of HPO terms
    1 HPO term per row
obo
    Path to obo file for HPO annotations
"""

### ---------------------------------------- ###

def parse_args():
    
    obo_path = argv[argv.index('--obo') + 1]
    hpo_id_to_name, hpo_hierarchy = parse_obo(obo_path)

    list1_hpo_terms_path = argv[argv.index('--list1') + 1]
    list1_hpo_terms = [term for term in open(list1_hpo_terms_path, 'r').read().split('\n') if len(term)]

    list2_hpo_terms_path = argv[argv.index('--list2') + 1]
    list2_hpo_terms = [term for term in open(list2_hpo_terms_path, 'r').read().split('\n') if len(term)]
    
    return hpo_id_to_name, hpo_hierarchy, list1_hpo_terms, list2_hpo_terms

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

def hpo_overlap(hpo_id_to_name_dict, hierarchy, hits, targets):
    
    impossible_distance = 100000
    
    distances = []
    
    for h in hits:
        
        h_name = hpo_id_to_name_dict[h]
        
        if h in targets:
            
            distances.append((h, h_name, h, h_name, 0, 'y'))
            
            continue
        
        else:
            
            h_parents = get_parent_terms([h], hierarchy)
            
            h_children = get_children_terms([h], hierarchy)
            
            min_dist = ('', impossible_distance)
            
            for t in targets:
                
                if t in h_parents:
                    
                    dist = get_parent_child_distance(t, h, hierarchy)
                    
                    if dist < min_dist[1]:
                        
                        min_dist = (t, dist, 'y')
                
                elif t in h_children:

                    dist = get_parent_child_distance(h, t, hierarchy)
                    
                    if dist < min_dist[1]:
                        
                        min_dist = (t, dist, 'y')
                    
                else:
                    
                    _, dist = find_closest_common_ancestor(h, t, hierarchy)
                    
                    if dist < min_dist[1]:
                        
                        min_dist = (t, dist, 'n')
            
            if min_dist[1] == impossible_distance:
                
                min_dist = (h, h_name, 'NA', 'NA', 'NA', 'NA')
            
            else:
                
                min_dist = (h, h_name, min_dist[0], hpo_id_to_name_dict[min_dist[0]], min_dist[1], min_dist[2])
            
            distances.append(min_dist)
    
    distances = pd.DataFrame(distances, columns=['hit_hpo', 'hit_hpo_name', 'closest_target_hpo', 'closest_target_hpo_name', 'distance', 'same_branch'])
    
    return distances

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

from sys import argv

### Load data

hpo_id_to_name, hpo_hierarchy, list1_hpo_terms, list2_hpo_terms = parse_args()

### Compute distance of list1 to list2

distances = hpo_overlap(hpo_id_to_name, hpo_hierarchy, list1_hpo_terms, list2_hpo_terms)

distances.to_csv('list1_vs_list2_hpo_distance.tsv', sep='\t', index=False)

### Compute distance of list2 to list1

distances = hpo_overlap(hpo_id_to_name, hpo_hierarchy, list2_hpo_terms, list1_hpo_terms)

distances.to_csv('list2_vs_list1_hpo_distance.tsv', sep='\t', index=False)
