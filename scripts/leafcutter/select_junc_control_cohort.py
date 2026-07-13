#!/usr/bin/env python3

### IMPORTS -------------------------------- ###

import numpy as np
import pandas as pd
import pickle as pk
import seaborn as sns

from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from sys import argv

### FUNCTIONS ------------------------------ ###

def parse_args() -> tuple[str, str]:
    
    # Sample files manifest
    sample_manifest = argv[argv.index('--samples') + 1]
    
    # Control files manifest
    control_manifest = argv[argv.index('--controls') + 1]
    
    # Min reads to support junctions
    min_reads = int(argv[argv.index('--min_reads') + 1]) if '--min_reads' in argv else 50
    
    # Max control cohort size
    max_ctrls = int(argv[argv.index('--max_ctrls') + 1]) if '--max_ctrls' in argv else 50
    
    return sample_manifest, control_manifest, min_reads, max_ctrls

### ---------------------------------------- ###

def load_junc_files(manifest_path: str, min_reads: int=50) -> tuple[list[str], pd.DataFrame]:
    
    sample_ids = []
    junctions = pd.DataFrame({col : [] for col in ['contig', 'start', 'end', 'strand']})
    merge_outer_toggle = True
    for s_path in open(manifest_path, 'r').readlines():
        # Fix path and define sample id
        s_path = s_path.strip()
        s_id = s_path.split('/')[-1].replace('.junc', '')
        sample_ids.append(s_id)
        # Load and filter data
        s_data = pd.read_csv(s_path, sep='\t', header=None)
        s_data.columns = ['contig', 'start', 'end', 'id', s_id, 'strand']
        s_data_tot_reads = s_data[s_id].sum()
        s_data = s_data.loc[s_data[s_id] >= min_reads, ['contig', 'start', 'end', 'strand', s_id]]
        s_data[s_id] = s_data[s_id] * (1e6 / s_data_tot_reads)
        # Merge to junction df
        if merge_outer_toggle:
            junctions = pd.merge(junctions, s_data, how='outer', on=['contig', 'start', 'end', 'strand'])
            merge_outer_toggle = False
        else:
            junctions = pd.merge(junctions, s_data, how='inner', on=['contig', 'start', 'end', 'strand'])
    
    return sample_ids, junctions

### ---------------------------------------- ###

def normalize_counts(cnts):
    
    #cnts = np.log1p(cnts * (1e6 / cnts.sum(axis=0)))
    cnts = np.log1p(cnts)
    
    return cnts

### ---------------------------------------- ###

def reduce_dimensions(cnts: pd.DataFrame) -> tuple[PCA, pd.DataFrame, np.array, int]:
    
    # Fit PCA
    pca_model = PCA()
    pca_data = pca_model.fit_transform(cnts.T)
    pca_data = pd.DataFrame(pca_data, index=cnts.columns, columns=[f'PC{i+1}' for i in range(pca_data.shape[1])])
    
    # Extract explained variance
    explained_variance = pca_model.explained_variance_
    explained_variance = explained_variance * (100 / explained_variance.sum())
    explained_variance = explained_variance[:20]
    
    # Extract optimal components
    optimal_components, _ = kneedle(explained_variance)
    optimal_components += 1
    
    # Plot explained variance
    plt.figure(figsize=(7, 4))
    plt.plot([str(i+1) for i in range(explained_variance.shape[0])], explained_variance / 100, 'b', marker='o', markersize=5, linewidth=1)
    plt.vlines(optimal_components + 0.5, 0, max(explained_variance) / 100, linestyle='dashed', color='red', linewidth=1)
    plt.text(optimal_components + 1, max(explained_variance) / 110, 'Optimal PCA components', color='red')
    plt.xlabel('PC')
    plt.ylabel('Explained Variance (%)')
    plt.tight_layout()
    plt.savefig('pca_explained_variance.png', bbox_inches='tight', dpi=600)
    plt.close()
    
    # Export data
    pca_data.to_csv('pca_transform.tsv', sep='\t', index=True, header=True)
    pk.dump(pca_model, open('pca.pkl', 'wb'))
    
    return pca_model, pca_data, explained_variance, optimal_components

### ---------------------------------------- ###

def kneedle(vector: np.array, sort_vector: bool=True) -> tuple[int, float]:
    
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

def plot_pca(pca_data: pd.DataFrame, explained_variance: np.array, info: pd.DataFrame, info_name: str='', out_plot: str='pca.png'):
    
    # Prepare data for plotting
    plot_data = pca_data.loc[:, ['PC1', 'PC2']].copy()
    plot_data.loc[:, info_name] = [info.loc[info['sample'] == sample, 'condition'].values[0] if sample in info['sample'].values else 'NA' for sample in pca_data.index]
    
    # Plot data
    plt.figure(figsize=(5, 5))
    ax = sns.scatterplot(data=plot_data,
                         x='PC1',
                         y='PC2',
                         hue=info_name,
                         palette='tab10',
                         s=50,
                         marker='o',
                         edgecolors='black',
                         linewidths=1)
    plt.xlabel(f'PC1 ({round(explained_variance[0], 2)}%)')
    plt.ylabel(f'PC2 ({round(explained_variance[1], 2)}%)')
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.savefig(out_plot, bbox_inches='tight', dpi=600)
    plt.close()

### MAIN ----------------------------------- ###

if __name__ == '__main__':
    
    # Parse CLI args
    sample_manifest, control_manifest, min_reads, max_ctrls = parse_args()
    
    # Load samples, filter for reads, then merge
    sample_ids, sample_junctions = load_junc_files(sample_manifest, min_reads)
    
    # Add controls
    control_ids, control_junctions = load_junc_files(control_manifest, min_reads)
    all_junctions = pd.merge(sample_junctions, control_junctions, how='inner', on=['contig', 'start', 'end', 'strand'])
    del sample_junctions, control_junctions
    
    if all_junctions.shape[0]:
    
        # PCA
        normalized_counts = normalize_counts(all_junctions[sample_ids + control_ids])
        _, pca_data, explained_variance, optimal_components = reduce_dimensions(normalized_counts)
        pca_data = pca_data.iloc[:, :optimal_components+1]
        
        # For each sample, find closest controls (weighted distance in PCA space)
        pca_data = pca_data.T
        selected_controls = set()
        for s_id in sample_ids:
            distances = np.sqrt(
                np.sum(
                    np.power(
                        pca_data[control_ids].sub(pca_data[s_id].values, axis=0).mul(explained_variance[:optimal_components+1] / explained_variance[:optimal_components+1].sum(), axis=0),
                        2
                    ),
                    axis=0
                )
            )
            distances = distances.sort_values()
            closest_samples = distances.index.values[:max_ctrls].tolist()
            selected_controls.update(closest_samples)
        
        # Plot PCA
        pca_data = pca_data.T
        sample_info = pd.DataFrame({
            'sample' : sample_ids + control_ids,
            'condition' : ['sample' for _ in sample_ids] + ['selected_control' if c_id in selected_controls else 'discarded_control' for c_id in control_ids]
        })
        plot_pca(pca_data, explained_variance, sample_info, 'sample_type', 'controls_selection.png')
        
        # Generate updated control manifest
        updated_control_manifest = '\n'.join([c_path.strip() for c_path,c_id in zip(open(control_manifest, 'r').readlines(), control_ids) if c_id in selected_controls])
    
    else:
        
        updated_control_manifest = ''
    
    # Output updated control manifest
    with open('junc_selected_control_manifest.txt', 'w') as manifest_out:
        manifest_out.write(updated_control_manifest)
