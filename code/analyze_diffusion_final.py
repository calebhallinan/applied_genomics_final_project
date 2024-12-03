# import packages

import scanpy as sc
from sklearn import metrics
import h5py
from scipy.sparse import csc_matrix
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# read in RDS file
from rds2py import read_rds
# import kmeans
from sklearn.cluster import KMeans
from pathlib import Path


############################################################################################################


### read in ground truth data ###

# metadata
metadata = pd.read_csv("/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/sc_meta.csv", index_col=0)
metadata

# gt matrix
gt_matrix = read_rds('/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/mtx_500_genes.RDS')
# make an adata object
adata_gt = sc.AnnData(gt_matrix.T)
adata_gt.obs = metadata.reset_index(drop=True)
adata_gt.var = pd.DataFrame(index=range(adata_gt.shape[1]))
# normalize and cluster
sc.pp.log1p(adata_gt)
sc.pp.neighbors(adata_gt)
sc.tl.leiden(adata_gt)
kmeans = KMeans(n_clusters=25, random_state=0).fit(adata_gt.X)
adata_gt.obs['kmeans'] = kmeans.labels_


############################################################################################################


### read in different conditions ###


# get all conditions
file_paths = [
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_1p1max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_1p01max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_1p05max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_1p2max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_1p3max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_1p4max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_1p5max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_1p6max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_1p7max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_1p8max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_1p9max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_2max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_3max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_4max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_5max_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_1pct_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_2pct_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_4pct_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_5pct_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_6pct_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_8pct_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_10pct_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_12pct_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_14pct_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_16pct_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_18pct_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_20pct_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_30pct_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_40pct_parallelized.RDS',
    '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_50pct_parallelized.RDS'
]



# Dictionary to store processed AnnData objects
adata_dict = {}

for file_path in file_paths:
    # Generate a dynamic variable name
    var_name = f"adata_{Path(file_path).stem.replace('ld_mtx_', '').replace('_parallelized', '')}"
    
    # Read data and create AnnData object
    data = read_rds(file_path)
    adata = sc.AnnData(data.T)
    adata.obs = metadata.reset_index(drop=True)
    adata.var = pd.DataFrame(index=range(data.shape[0]))
    
    # Preprocessing and clustering
    sc.pp.log1p(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)
    kmeans = KMeans(n_clusters=25, random_state=0).fit(adata.X)
    adata.obs['kmeans'] = kmeans.labels_
    
    # Save to dictionary
    adata_dict[var_name] = adata


############################################################################################################


### compare clustering results ###


# Initialize ARI DataFrame
ari_df = pd.DataFrame(columns=['ARI_leiden', 'ARI_kmeans'], index=['GT'] + list(adata_dict.keys()))

# Compute ARI for ground truth (self-comparison)
ari_df.loc['GT', 'ARI_leiden'] = metrics.adjusted_rand_score(adata_gt.obs['leiden'], adata_gt.obs['leiden'])
ari_df.loc['GT', 'ARI_kmeans'] = metrics.adjusted_rand_score(adata_gt.obs['kmeans'], adata_gt.obs['kmeans'])

# Compute ARI for each key in the dictionary
for condition, adata in adata_dict.items():
    ari_df.loc[condition, 'ARI_leiden'] = metrics.adjusted_rand_score(adata_gt.obs['leiden'], adata.obs['leiden'])
    ari_df.loc[condition, 'ARI_kmeans'] = metrics.adjusted_rand_score(adata_gt.obs['kmeans'], adata.obs['kmeans'])

# Print the ARI DataFrame
print(ari_df)

# change the index to be more readable
ari_df.index = [
    'Ground Truth',
    '1.1max',
    '1.01max',
    '1.05max',
    '1.2max',
    '1.3max',
    '1.4max',
    '1.5max',
    '1.6max',
    '1.7max',
    '1.8max',
    '1.9max',
    '2max',
    '3max',
    '4max',
    '5max',
    '1pct',
    '2pct',
    '4pct',
    '5pct',
    '6pct',
    '8pct',
    '10pct',
    '12pct',
    '14pct',
    '16pct',
    '18pct',
    '20pct',
    '30pct',
    '40pct',
    '50pct'
]

# reorder to make 1p01max and 1p05 after Ground Truth
ari_df = ari_df.reindex(['Ground Truth', '1.01max', '1.05max', '1.1max', '1.2max', '1.3max', '1.4max', '1.5max', '1.6max', 
                         '1.7max', '1.8max', '1.9max', '2max', '3max', '4max', '5max', '1pct', '2pct', '4pct', '5pct', 
                         '6pct', '8pct', '10pct', '12pct', '14pct', '16pct', '18pct', '20pct', '30pct', '40pct', '50pct'])

# separate the ari for max and pct
ari_df_max = ari_df.loc[ari_df.index.str.contains('max') | ari_df.index.str.contains('Truth')]
ari_df_pct = ari_df.loc[ari_df.index.str.contains('pct') | ari_df.index.str.contains('Truth')]


############################################################################################################


### plot ARI ###

# Plot ARI for leiden as a scatterplot with line connecting the points
plt.figure(figsize=(10, 6))
plt.scatter(ari_df_max.index, ari_df_max['ARI_leiden'], label='leiden')
plt.plot(ari_df_max.index, ari_df_max['ARI_leiden'], linestyle='-', marker='o')
plt.scatter(ari_df_max.index, ari_df_max['ARI_kmeans'], label='kmeans')
plt.plot(ari_df_max.index, ari_df_max['ARI_kmeans'], linestyle='-', marker='o')
plt.xticks(rotation=90)
plt.xlabel('Condition')
plt.ylabel('ARI')
plt.title('ARI for leiden and kmeans clustering')
plt.legend()


# Plot ARI for kmeans as a scatterplot with line connecting the points
plt.figure(figsize=(10, 6))
plt.scatter(ari_df_pct.index, ari_df_pct['ARI_leiden'], label='leiden')
plt.plot(ari_df_pct.index, ari_df_pct['ARI_leiden'], linestyle='-', marker='o')
plt.scatter(ari_df_pct.index, ari_df_pct['ARI_kmeans'], label='kmeans')
plt.plot(ari_df_pct.index, ari_df_pct['ARI_kmeans'], linestyle='-', marker='o')
plt.xticks(rotation=90)
plt.xlabel('Condition')
plt.ylabel('ARI')
plt.title('ARI for leiden and kmeans clustering')
plt.legend()



############################################################################################################


