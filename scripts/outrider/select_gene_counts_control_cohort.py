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
    
    # Counts matrix (first column must be gene names)
    counts_matrix_path = argv[argv.index('--counts') + 1]
    
    # Control files IDs
    control_manifest = argv[argv.index('--controls') + 1]
    
    # Min reads to support junctions
    min_reads = int(argv[argv.index('--min_reads') + 1]) if '--min_reads' in argv else 50
    
    # Max control cohort size
    max_ctrls = int(argv[argv.index('--max_ctrls') + 1]) if '--max_ctrls' in argv else 100
    
    return counts_matrix_path, control_manifest, min_reads, max_ctrls

### ---------------------------------------- ###

def normalize_counts(cnts):
    
    cnts = np.log1p(cnts * (1e6 / cnts.sum(axis=0)))
    
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
    counts_matrix_path, control_manifest, min_reads, max_ctrls = parse_args()
    
    # Load counts, then filter
    counts_matrix_raw = pd.read_csv(counts_matrix_path, sep='\t', index_col=0)
    counts_matrix_raw = counts_matrix_raw.loc[counts_matrix_raw.sum(axis=1) != 0,]
    counts_matrix = counts_matrix_raw.loc[(counts_matrix_raw >= min_reads).sum(axis=1) >= 3,]
    
    # Load control IDs
    control_ids = open(control_manifest, 'r').read().strip().split('\n')
    sample_ids = [col for col in counts_matrix.columns if col not in control_ids]
    
    if len(control_ids):
    
        # PCA
        normalized_counts = normalize_counts(counts_matrix)
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
        
        # Output updated control manifest
        updated_control_manifest = '\n'.join([c_id for c_id in control_ids if c_id in selected_controls])
        with open('ctrl_ids.filtered.txt', 'w') as manifest_out:
            manifest_out.write(updated_control_manifest)
        
        # Output updated counts matrix
        updated_counts_matrix = counts_matrix_raw[sample_ids + list(selected_controls)]
        updated_counts_matrix.to_csv('MergedGeneCounts_Subset.tsv', sep='\t', index=True)
    
    else:
        
        # Output mock control manifest
        updated_control_manifest = ''
        with open('mock.filtered.txt', 'w') as manifest_out:
            manifest_out.write(updated_control_manifest)
        
        # Output same raw counts matrix as input
        counts_matrix_raw.to_csv('MergedGeneCounts_Subset.tsv', sep='\t', index=True)
