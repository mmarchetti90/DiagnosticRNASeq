#!/usr/bin/env python3

"""
This script clusters HPO terms based on their distance
Only terms in the same branch are clustered
Clusters are based on graph leaves (i.e. it used terms with no childrens among the other HPO terms)

N.B. Algorithm is stringent

Arguments
---------
hpo_list
    Path to txt file with list of HPO terms
    1 HPO term per row
obo
    Path to obo file for HPO annotations
max_hpo_distance
    Maximum distance between HPO terms to allow clustering
    (Default = 3)
"""

### ---------------------------------------- ###

def parse_args():
    
    obo_path = argv[argv.index('--obo') + 1]
    hpo_id_to_name, hpo_hierarchy = parse_obo(obo_path)

    hpo_terms_path = argv[argv.index('--hpo_list') + 1]
    hpo_terms = [term for term in open(hpo_terms_path, 'r').read().split('\n') if len(term)]
    hpo_terms = list(set(hpo_terms))
    
    if '--max_hpo_distance' in argv:
        
        max_hpo_distance = int(argv[argv.index('--max_hpo_distance') + 1])
        
    else:
        
        max_hpo_distance = 3
    
    return hpo_id_to_name, hpo_hierarchy, hpo_terms, max_hpo_distance

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

def hpo_distance(hpo_id_to_name_dict, hierarchy, terms):
    
    impossible_distance = -1
    
    distances = np.zeros((len(terms), len(terms)))
    
    for h1 in terms:
        
        idx1 = terms.index(h1)
        
        h1_parents = get_parent_terms([h1], hierarchy)
        
        h1_children = get_children_terms([h1], hierarchy)
        
        for h2 in terms[terms.index(h1) + 1:]:
            
            idx2 = terms.index(h2)
            
            if h2 in h1_parents:
                
                dist = get_parent_child_distance(h2, h1, hierarchy)
            
            elif h2 in h1_children:
                
                dist = get_parent_child_distance(h1, h2, hierarchy)
            
            else:
                
                dist = impossible_distance
            
            distances[idx1, idx2] = dist
            
            distances[idx2, idx1] = dist
            
    distances = pd.DataFrame(distances, index=terms, columns=terms)
    
    return distances

### ---------------------------------------- ###

def cluster_terms(terms, hierarchy, distances, max_distance=3):
    
    log = []
    
    # Find childrenless terms (i.e. leaves)
    # N.B. a node with childrens is still considered a leaf if none of the children is in hpo_terms
    
    leaves = []
    
    for h in terms:
        
        childrens = [c for c in get_children_terms([h], hierarchy) if c in terms]
        
        if len(childrens) == 1:
            
            leaves.append(h)
    
    original_leaves_num = len(leaves)
    
    log.append(f'Found {original_leaves_num} initial leaves')
    
    # Temporarily remove terms shared by leaves
    
    sharing = {}
    
    for l in leaves:
        
        parents = get_parent_terms([l], hierarchy)
        
        for p in parents:
            
            if p not in terms:
                
                continue
            
            if p in sharing.keys():
                
                sharing[p] += 1
            
            else:
                
                sharing[p] = 1
    
    to_be_removed = [s for s,sn in sharing.items() if sn > 1]
    
    terms = [h for h in terms if h not in to_be_removed]
    
    log.append(f'Temporarily removed {len(to_be_removed)} HPO terms shared between leaves ancestors')
    
    # Init clusters
    
    clusters = {h : (leaves.index(h) if h in leaves
                     else -1)
                for h in terms}
    
    clusters_parents = {cl : h for h,cl in clusters.items() if cl != -1}
    
    # Cluster
    
    toggle = True
    
    while toggle:
        
        old_clusters = clusters.copy()
        
        for l in leaves:
            
            # Extract cluster elements
            
            l_clust = clusters[l]
            
            cluster_elements = [h for h,cl in clusters.items() if cl == l_clust]
            
            cluster_parent = clusters_parents[l_clust]
            
            # Find closest non-cluster element which is a parent of cluster_parent
            
            closest_term = [(p, distances.loc[p, cluster_parent]) for p in get_parent_terms([cluster_parent], hierarchy) if p in terms and p not in cluster_elements]
            
            closest_term.sort(key=lambda ct: ct[1])
            
            if not len(closest_term):
                
                continue
            
            else:
                
                closest_term, closest_term_distance = closest_term[0]
                
                # If closest term is within max_distance from cluster_parent, then add to cluster
                # Else, create new leaf
                
                if closest_term_distance <= max_distance:
                
                    clusters[closest_term] = l_clust
                    
                    clusters_parents[l_clust] = closest_term
                    
                    break
                
                else:
                    
                    leaves.append(closest_term)
                    
                    new_cluster_n = max(clusters.values()) + 1
                    
                    clusters[closest_term] = new_cluster_n
                    
                    clusters_parents[new_cluster_n] = closest_term
                    
                    break
        
        if clusters == old_clusters:
            
            toggle = False
        
    log.append(f'Created {len(leaves) - original_leaves_num} new leaves')
        
    log.append(f'Found {max(clusters.values()) + 1} clusters')
    
    # Keep shared terms that don't have any children among other shared terms
    
    shared_terms = to_be_removed
    
    shared_terms = [s for s in shared_terms if np.isin(list(get_children_terms([s], hierarchy)), shared_terms).sum() == 1]
    
    # Simplify clusters using shared parent terms previously removed
    
    for s in shared_terms:
        
        # Find clusters sharing s
        
        s_children = get_children_terms([s], hierarchy)
    
        clusters_to_merge = [cl for cl,h in clusters_parents.items() if h in s_children]
        
        # Only keep those whose distance from s is <= max_distance
        
        clusters_to_merge = [cl for cl in clusters_to_merge if get_parent_child_distance(s, clusters_parents[cl], hierarchy) <= max_distance]
        
        if len(clusters_to_merge):
            
            new_cluster_n = min(clusters_to_merge)
            
            # Update clusters
            
            elements_to_update = [h for h,cl in clusters.items() if cl in clusters_to_merge] + [s]
            
            for e in elements_to_update:
                
                clusters[e] = new_cluster_n
            
            # Update clusters_parents
            
            for cl in clusters_to_merge:
                
                if cl == new_cluster_n:
                    
                    clusters_parents[new_cluster_n] = s
                
                else:
                    
                    _ = clusters_parents.pop(cl, None)
            
            # Log
            
            log.append(f'Merged clusters {clusters_to_merge} into cluster {new_cluster_n}')
    
    # Add unclustered terms as individual clusters
    
    for term,cl in clusters.items():
        
        if cl == -1:
            
            new_cl = max(clusters.values()) + 1
            
            clusters[term] = new_cl
            
            clusters_parents[new_cl] = term
    
    # Re-index clusters
    
    reindexing = {old : new for new,old in enumerate(np.sort(list(clusters_parents.keys())))}
    
    clusters = {h : reindexing[cl] for h,cl in clusters.items() if cl != -1}
    
    clusters_parents = {reindexing[cl] : h for cl,h in clusters_parents.items()}

    log.append(f'Reduced data to {max(clusters.values()) + 1} clusters')
    
    return clusters, clusters_parents, log

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

import numpy as np
import pandas as pd

from sys import argv

### Load data

hpo_id_to_name, hpo_hierarchy, hpo_terms, max_hpo_distance = parse_args()

### Compute pairwise distances

pairwise_distances = hpo_distance(hpo_id_to_name, hpo_hierarchy, hpo_terms)

### Cluster

hpo_clusters, hpo_clusters_parents, clustering_log = cluster_terms(hpo_terms, hpo_hierarchy, pairwise_distances, max_hpo_distance)

### Simplifying HPO terms

simplified_hpo_terms = list(hpo_clusters_parents.values())

### Saving data

with open('simplified_hpo_terms.txt', 'w') as hpo_out:
    
    hpo_out.write('\n'.join(simplified_hpo_terms))

with open('hpo_clusters.txt', 'w') as clusters_out:
    
    clusters_out.write('\n'.join(['hpo\thpo_name\tcluster'] + [f'{h}\t{hpo_id_to_name[h]}\t{cl}' for h,cl in hpo_clusters.items()]))

with open('clustering_log.txt', 'w') as log_out:
    
    log_out.write('\n'.join(clustering_log))
