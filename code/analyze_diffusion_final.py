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
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.mixture import GaussianMixture
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
sc.tl.pca(adata_gt, n_comps=30)
sc.pp.neighbors(adata_gt, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata_gt)
# test more resolutions
# sc.tl.leiden(adata_gt, key_added="leiden_res0_25", resolution=0.25)
# sc.tl.leiden(adata_gt, key_added="leiden_res0_5", resolution=0.5)
# sc.tl.leiden(adata_gt, key_added="leiden_res1", resolution=1.0)
# sc.tl.leiden(adata_gt, key_added="leiden_res2", resolution=2.0)
sc.tl.louvain(adata_gt)
kmeans = KMeans(n_clusters=25, random_state=0).fit(adata_gt.obsm['X_pca'])
adata_gt.obs['kmeans'] = kmeans.labels_
# perform agglomerative clustering
agglomerative = AgglomerativeClustering(n_clusters=25).fit(adata_gt.obsm['X_pca'])
adata_gt.obs['agglomerative'] = agglomerative.labels_
# perform GMM
gmm = GaussianMixture(n_components=25, random_state=0).fit(adata_gt.obsm['X_pca'])
adata_gt.obs['gmm'] = gmm.predict(adata_gt.obsm['X_pca'])


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
    sc.tl.pca(adata, n_comps=30)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(adata)
    # test more resolutions
    # sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
    # sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
    # sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)
    # sc.tl.leiden(adata, key_added="leiden_res2", resolution=2.0)
    sc.tl.louvain(adata)
    kmeans = KMeans(n_clusters=25, random_state=0).fit(adata.obsm['X_pca'])
    adata.obs['kmeans'] = kmeans.labels_
    # perform agglomerative clustering
    agglomerative = AgglomerativeClustering(n_clusters=25).fit(adata.obsm['X_pca'])
    adata.obs['agglomerative'] = agglomerative.labels_
    # perform GMM
    gmm = GaussianMixture(n_components=25, random_state=0).fit(adata.obsm['X_pca'])
    adata.obs['gmm'] = gmm.predict(adata.obsm['X_pca'])

    # Save to dictionary
    adata_dict[var_name] = adata


############################################################################################################


### compare clustering results ###


# Initialize ARI DataFrame
ari_df = pd.DataFrame(columns=['ARI_leiden', 'ARI_kmeans', 'ARI_louvain', 'ARI_agglomerative'], index=['GT'] + list(adata_dict.keys()))

# Compute ARI for ground truth (self-comparison)
ari_df.loc['GT', 'ARI_leiden'] = metrics.adjusted_rand_score(adata_gt.obs['leiden'], adata_gt.obs['leiden'])
ari_df.loc['GT', 'ARI_kmeans'] = metrics.adjusted_rand_score(adata_gt.obs['kmeans'], adata_gt.obs['kmeans'])
ari_df.loc['GT', 'ARI_louvain'] = metrics.adjusted_rand_score(adata_gt.obs['louvain'], adata_gt.obs['louvain'])
ari_df.loc['GT', 'ARI_agglomerative'] = metrics.adjusted_rand_score(adata_gt.obs['agglomerative'], adata_gt.obs['agglomerative'])
ari_df.loc['GT', 'ARI_gmm'] = metrics.adjusted_rand_score(adata_gt.obs['gmm'], adata_gt.obs['gmm'])

# add new leiden resolutions
# ari_df.loc['GT', 'ARI_leiden_res0_25'] = metrics.adjusted_rand_score(adata_gt.obs['leiden_res0_25'], adata_gt.obs['leiden_res0_25'])
# ari_df.loc['GT', 'ARI_leiden_res0_5'] = metrics.adjusted_rand_score(adata_gt.obs['leiden_res0_5'], adata_gt.obs['leiden_res0_5'])
# ari_df.loc['GT', 'ARI_leiden_res1'] = metrics.adjusted_rand_score(adata_gt.obs['leiden_res1'], adata_gt.obs['leiden_res1'])
# ari_df.loc['GT', 'ARI_leiden_res2'] = metrics.adjusted_rand_score(adata_gt.obs['leiden_res2'], adata_gt.obs['leiden_res2'])

# Compute ARI for each key in the dictionary
for condition, adata in adata_dict.items():
    ari_df.loc[condition, 'ARI_leiden'] = metrics.adjusted_rand_score(adata_gt.obs['leiden'], adata.obs['leiden'])
    ari_df.loc[condition, 'ARI_kmeans'] = metrics.adjusted_rand_score(adata_gt.obs['kmeans'], adata.obs['kmeans'])
    ari_df.loc[condition, 'ARI_louvain'] = metrics.adjusted_rand_score(adata_gt.obs['louvain'], adata.obs['louvain'])
    ari_df.loc[condition, 'ARI_agglomerative'] = metrics.adjusted_rand_score(adata_gt.obs['agglomerative'], adata.obs['agglomerative'])
    ari_df.loc[condition, 'ARI_gmm'] = metrics.adjusted_rand_score(adata_gt.obs['gmm'], adata.obs['gmm'])
    # # add new leiden resolutions
    # ari_df.loc[condition, 'ARI_leiden_res0_25'] = metrics.adjusted_rand_score(adata_gt.obs['leiden_res0_25'], adata.obs['leiden_res0_25'])
    # ari_df.loc[condition, 'ARI_leiden_res0_5'] = metrics.adjusted_rand_score(adata_gt.obs['leiden_res0_5'], adata.obs['leiden_res0_5'])
    # ari_df.loc[condition, 'ARI_leiden_res1'] = metrics.adjusted_rand_score(adata_gt.obs['leiden_res1'], adata.obs['leiden_res1'])
    # ari_df.loc[condition, 'ARI_leiden_res2'] = metrics.adjusted_rand_score(adata_gt.obs['leiden_res2'], adata.obs['leiden_res2'])

# Print the ARI DataFrame
print(ari_df)

adata_gt.obs['louvain'].value_counts()

adata_dict['adata_10pct'].obs['louvain'].value_counts()

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

ari_df_max.index
x_axis = [.9, 1.01, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3, 3.5]

### plot ARI ###

# Plot ARI for leiden as a scatterplot with line connecting the points
plt.figure(figsize=(14, 6))
plt.scatter(x_axis, ari_df_max['ARI_leiden'], label='Leiden', s=100)
plt.plot(x_axis, ari_df_max['ARI_leiden'], linestyle='-', marker='o')
# add new leiden resolutions
# plt.scatter(x_axis, ari_df_max['ARI_leiden_res0_25'], label='Leiden_res0_25', s=100)
# plt.plot(x_axis, ari_df_max['ARI_leiden_res0_25'], linestyle='-', marker='o')
# plt.scatter(x_axis, ari_df_max['ARI_leiden_res0_5'], label='Leiden_res0_5', s=100)
# plt.plot(x_axis, ari_df_max['ARI_leiden_res0_5'], linestyle='-', marker='o')
# plt.scatter(x_axis, ari_df_max['ARI_leiden_res1'], label='Leiden_res1', s=100)
# plt.plot(x_axis, ari_df_max['ARI_leiden_res1'], linestyle='-', marker='o')
# plt.scatter(x_axis, ari_df_max['ARI_leiden_res2'], label='Leiden_res2', s=100)
# plt.plot(x_axis, ari_df_max['ARI_leiden_res2'], linestyle='-', marker='o')
plt.scatter(x_axis, ari_df_max['ARI_kmeans'], label='K-means', s=100)
plt.plot(x_axis, ari_df_max['ARI_kmeans'], linestyle='-', marker='o')
plt.scatter(x_axis, ari_df_max['ARI_louvain'], label='Louvain', s=100)
plt.plot(x_axis, ari_df_max['ARI_louvain'], linestyle='-', marker='o')
plt.scatter(x_axis, ari_df_max['ARI_agglomerative'], label='Agglomerative', s=100)
plt.plot(x_axis, ari_df_max['ARI_agglomerative'], linestyle='-', marker='o')
# plt.scatter(x_axis, ari_df_max['ARI_gmm'], label='GMM', s=100)
# plt.plot(x_axis, ari_df_max['ARI_gmm'], linestyle='-', marker='o')
plt.xlabel('Condition', fontsize=14)
plt.ylabel('Adjusted Rand Index (ARI)', fontsize=14) 
plt.title('ARI for Max Distance a transcript can diffuse', fontsize=22)
# make the legened bigger
plt.legend(fontsize=12)
sns.despine()
# make size of text larger
plt.yticks(fontsize=14)
plt.ylim(.4,1.02)
# change the x ticks to be ari_df_max.index
# Set custom labels for the x-axis
plt.xticks(ticks=x_axis, labels=ari_df_max.index, rotation=90,fontsize=12)



x_axis1 = [0, 1, 2, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 23, 26, 29]

# Plot ARI for kmeans as a scatterplot with line connecting the points
plt.figure(figsize=(14, 6))
plt.scatter(x_axis1, ari_df_pct['ARI_leiden'], label='Leiden', s=100)
plt.plot(x_axis1, ari_df_pct['ARI_leiden'], linestyle='-', marker='o')
plt.scatter(x_axis1, ari_df_pct['ARI_kmeans'], label='K-means', s=100)
plt.plot(x_axis1, ari_df_pct['ARI_kmeans'], linestyle='-', marker='o')
plt.scatter(x_axis1, ari_df_pct['ARI_louvain'], label='Louvain', s=100)
plt.plot(x_axis1, ari_df_pct['ARI_louvain'], linestyle='-', marker='o')
plt.scatter(x_axis1, ari_df_pct['ARI_agglomerative'], label='Agglomerative', s=100)
plt.plot(x_axis1, ari_df_pct['ARI_agglomerative'], linestyle='-', marker='o')
plt.xlabel('Condition', fontsize=14)
plt.ylabel('Adjusted Rand Index (ARI)', fontsize=14)
plt.title('ARI for Total Percentage of transcripts that can diffuse', fontsize=22)
plt.legend(fontsize=12)
sns.despine()
# make size of text larger
plt.yticks(fontsize=12)
plt.ylim(.4,1.02)
plt.xticks(ticks=x_axis1, labels=ari_df_pct.index, rotation=90,fontsize=12)


############################################################################################################


import pandas as pd
import numpy as np
from scipy.optimize import linear_sum_assignment
import seaborn as sns
import matplotlib.pyplot as plt
import colorcet as cc

# Extract cluster labels
gt_clusters = adata_gt.obs['kmeans']
other_clusters = adata_dict["adata_1pct"].obs['kmeans']

# Create a contingency table
contingency = pd.crosstab(gt_clusters, other_clusters)

# Solve the assignment problem (maximize overlap, so minimize negative values)
row_ind, col_ind = linear_sum_assignment(-contingency.values)

# Create a mapping from other_clusters to gt_clusters
cluster_mapping = {contingency.columns[col]: contingency.index[row] for row, col in zip(row_ind, col_ind)}

# Map clusters in the second dataset to ground truth clusters
adata_dict["adata_1pct"].obs['kmeans_mapped'] = adata_dict["adata_1pct"].obs['kmeans'].map(cluster_mapping).fillna(-1)

# Create a shared color palette for all clusters
all_clusters = np.union1d(adata_gt.obs['kmeans'].unique(), adata_dict["adata_1pct"].obs['kmeans_mapped'].unique())
palette = sns.color_palette(cc.glasbey, n_colors=len(all_clusters))
hue_mapping = {cluster: color for cluster, color in zip(all_clusters, palette)}



# Plot the ground truth
plt.figure(figsize=(12, 8))
sns.scatterplot(
    x=adata_gt.obs['x'],
    y=adata_gt.obs['y'],
    hue=adata_gt.obs['kmeans'],
    size=adata_gt.obs["area"],
    palette=hue_mapping
)
plt.title("Ground Truth")
sns.despine()
# Do not plot legend
plt.legend([], [], frameon=False)

# Plot the different conditions with mapped clusters
plt.figure(figsize=(12, 8))
sns.scatterplot(
    x=adata_dict["adata_1pct"].obs['x'],
    y=adata_dict["adata_1pct"].obs['y'],
    hue=adata_dict["adata_1pct"].obs['kmeans_mapped'],
    size=adata_dict["adata_1pct"].obs["area"],
    palette=hue_mapping
)
plt.title("1pct Clusters")
sns.despine()
# Do not plot legend
plt.legend([], [], frameon=False)
plt.show()




# Extract cluster labels
gt_clusters = adata_gt.obs['kmeans']
other_clusters = adata_dict["adata_50pct"].obs['kmeans']

# Create a contingency table
contingency = pd.crosstab(gt_clusters, other_clusters)

# Solve the assignment problem (maximize overlap, so minimize negative values)
row_ind, col_ind = linear_sum_assignment(-contingency.values)

# Create a mapping from other_clusters to gt_clusters
cluster_mapping = {contingency.columns[col]: contingency.index[row] for row, col in zip(row_ind, col_ind)}

# Map clusters in the second dataset to ground truth clusters
adata_dict["adata_50pct"].obs['kmeans_mapped'] = adata_dict["adata_50pct"].obs['kmeans'].map(cluster_mapping).fillna(-1)

# Create a shared color palette for all clusters
all_clusters = np.union1d(adata_gt.obs['kmeans'].unique(), adata_dict["adata_50pct"].obs['kmeans_mapped'].unique())
palette = sns.color_palette(cc.glasbey, n_colors=len(all_clusters))
hue_mapping = {cluster: color for cluster, color in zip(all_clusters, palette)}



# Plot the different conditions with mapped clusters
plt.figure(figsize=(12, 8))
sns.scatterplot(
    x=adata_dict["adata_50pct"].obs['x'],
    y=adata_dict["adata_50pct"].obs['y'],
    hue=adata_dict["adata_50pct"].obs['kmeans_mapped'],
    size=adata_dict["adata_50pct"].obs["area"],
    palette=hue_mapping
)
plt.title("50pct Clusters")
sns.despine()
# Do not plot legend
plt.legend([], [], frameon=False)
plt.show()




############################################################################################################


# plot umap for ground truth
sc.tl.umap(adata_gt)
sc.tl.umap(adata_dict["adata_50pct"])

# plot umap for ground truth
sc.pl.umap(adata_gt, color='leiden', palette='tab20', title='Ground Truth with Leiden', show=False)
# Customize the legend position
plt.legend(
    loc='upper center', 
    bbox_to_anchor=(0.5, -0.05),  # Position the legend below the plot
    ncol=4,  # Number of columns for the legend
    title='Leiden Clusters'  # Optional: Add a title to the legend
)
plt.show()



# plot umap for 50pct
sc.pl.umap(adata_dict["adata_50pct"], color='leiden', palette='tab20', title='50pct with Leiden Clusters', show=False)
# Customize the legend position
plt.legend(
    loc='upper center', 
    bbox_to_anchor=(0.5, -0.05),  # Position the legend below the plot
    ncol=4,  # Number of columns for the legend
    title='Leiden Clusters'  # Optional: Add a title to the legend
)
plt.show()




