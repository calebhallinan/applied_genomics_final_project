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



############################################################################################################

# set env and read in data
datadir = '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/'
# metadata
metadata = pd.read_csv("/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/sc_srt_meta.csv", index_col=0)
metadata

# Open the HDF5 file
file_path = f'{datadir}/filtered_srt_mtx.h5'
with h5py.File(file_path, 'r') as h5file:
    # Access the datasets under 'mm10-1.2.0_premrna'
    data = h5file['mm10-1.2.0_premrna/data'][:]
    indices = h5file['mm10-1.2.0_premrna/indices'][:]
    indptr = h5file['mm10-1.2.0_premrna/indptr'][:]
    shape = h5file['mm10-1.2.0_premrna/shape'][:]

    # Create a sparse matrix using the loaded data
    sparse_matrix = csc_matrix((data, indices, indptr), shape=shape)

# Print the sparse matrix to verify it was loaded correctly
# print(sparse_matrix)
sparse_matrix.toarray().shape

# make an adata object
adata = sc.AnnData(sparse_matrix.T)
adata.obs = metadata.reset_index(drop=True)
adata.var = pd.DataFrame(index=range(adata.shape[1]))


############################################################################################################


# read in RDS file
from rds2py import read_rds

ld10_d3 = read_rds('/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/filtered_ld10_d3_mtx.RDS')
ld10 = read_rds('/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/filtered_ld10_mtx.RDS')
ld20 = read_rds('/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/filtered_ld20_mtx.RDS')

# make adata objects for each
adata_ld10_d3 = sc.AnnData(ld10_d3.T)
adata_ld10_d3.obs = metadata.reset_index(drop=True)
adata_ld10_d3.var = pd.DataFrame(index=range(adata_ld10_d3.shape[1]))

adata_ld10 = sc.AnnData(ld10.T)
adata_ld10.obs = metadata.reset_index(drop=True)
adata_ld10.var = pd.DataFrame(index=range(adata_ld10.shape[1]))

adata_ld20 = sc.AnnData(ld20.T)
adata_ld20.obs = metadata.reset_index(drop=True)
adata_ld20.var = pd.DataFrame(index=range(adata_ld20.shape[1]))



############################################################################################################

# log normalize the data
sc.pp.log1p(adata)
sc.pp.log1p(adata_ld10_d3)
sc.pp.log1p(adata_ld10)
sc.pp.log1p(adata_ld20)

# cluster the data

sc.pp.neighbors(adata)
sc.tl.leiden(adata)

sc.pp.neighbors(adata_ld10_d3)
sc.tl.leiden(adata_ld10_d3)

sc.pp.neighbors(adata_ld10)
sc.tl.leiden(adata_ld10)

sc.pp.neighbors(adata_ld20)
sc.tl.leiden(adata_ld20)

# kmeans clustering
from sklearn.cluster import KMeans

kmeans = KMeans(n_clusters=25, random_state=0).fit(adata.X)
adata.obs['kmeans'] = kmeans.labels_

kmeans = KMeans(n_clusters=25, random_state=0).fit(adata_ld10_d3.X)
adata_ld10_d3.obs['kmeans'] = kmeans.labels_

kmeans = KMeans(n_clusters=25, random_state=0).fit(adata_ld10.X)
adata_ld10.obs['kmeans'] = kmeans.labels_

kmeans = KMeans(n_clusters=25, random_state=0).fit(adata_ld20.X)
adata_ld20.obs['kmeans'] = kmeans.labels_


############################################################################################################


# compare the two using ARI

# make a table of ARI values
ari_df = pd.DataFrame(columns=['ARI'], index=['LD10', 'LD10_D3', 'LD20'])
ari_df.loc['LD10', 'ARI'] = metrics.adjusted_rand_score(adata.obs['leiden'], adata_ld10.obs['leiden'])
ari_df.loc['LD10_D3', 'ARI'] = metrics.adjusted_rand_score(adata.obs['leiden'], adata_ld10_d3.obs['leiden'])
ari_df.loc['LD20', 'ARI'] = metrics.adjusted_rand_score(adata.obs['leiden'], adata_ld20.obs['leiden'])
ari_df


# make a table of ARI values for kmeans
ari_df = pd.DataFrame(columns=['ARI'], index=['LD10', 'LD10_D3', 'LD20'])
ari_df.loc['LD10', 'ARI'] = metrics.adjusted_rand_score(adata.obs['kmeans'], adata_ld10.obs['kmeans'])
ari_df.loc['LD10_D3', 'ARI'] = metrics.adjusted_rand_score(adata.obs['kmeans'], adata_ld10_d3.obs['kmeans'])
ari_df.loc['LD20', 'ARI'] = metrics.adjusted_rand_score(adata.obs['kmeans'], adata_ld20.obs['kmeans'])
ari_df



import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment

# Create subplots for side-by-side comparison
fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# Plot for LD10
crosstab = pd.crosstab(adata_ld10.obs['leiden'], adata.obs['leiden'], normalize="columns")
cost_matrix = crosstab.max().max() - crosstab.values
row_ind, col_ind = linear_sum_assignment(cost_matrix)
optimal_crosstab = crosstab.iloc[row_ind, :].iloc[:, col_ind]
sns.heatmap(optimal_crosstab, cmap='coolwarm', ax=axes[0], cbar=False)
axes[0].set_xlabel("Ground Truth Clusters", fontsize=14)
axes[0].set_ylabel("Predicted Clusters", fontsize=14)
axes[0].set_title("Normalized Confusion Matrix - leiden (LD10)", fontsize=18)
axes[0].tick_params(axis='both', which='major', labelsize=10)

# Plot for LD10_D3
crosstab = pd.crosstab(adata_ld10_d3.obs['leiden'], adata.obs['leiden'], normalize="columns")
cost_matrix = crosstab.max().max() - crosstab.values
row_ind, col_ind = linear_sum_assignment(cost_matrix)
optimal_crosstab = crosstab.iloc[row_ind, :].iloc[:, col_ind]
sns.heatmap(optimal_crosstab, cmap='coolwarm', ax=axes[1], cbar=False)
axes[1].set_xlabel("Ground Truth Clusters", fontsize=14)
axes[1].set_ylabel("Predicted Clusters", fontsize=14)
axes[1].set_title("Normalized Confusion Matrix - leiden (LD10_D3)", fontsize=18)
axes[1].tick_params(axis='both', which='major', labelsize=10)

# Plot for LD20
crosstab = pd.crosstab(adata_ld20.obs['leiden'], adata.obs['leiden'], normalize="columns")
cost_matrix = crosstab.max().max() - crosstab.values
row_ind, col_ind = linear_sum_assignment(cost_matrix)
optimal_crosstab = crosstab.iloc[row_ind, :].iloc[:, col_ind]
sns.heatmap(optimal_crosstab, cmap='coolwarm', ax=axes[2], cbar=False)
axes[2].set_xlabel("Ground Truth Clusters", fontsize=14)
axes[2].set_ylabel("Predicted Clusters", fontsize=14)
axes[2].set_title("Normalized Confusion Matrix - leiden (LD20)", fontsize=18)
axes[2].tick_params(axis='both', which='major', labelsize=10)

plt.tight_layout()
plt.show()




# Create subplots for side-by-side comparison
fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# Plot for LD10
crosstab = pd.crosstab(adata_ld10.obs['kmeans'], adata.obs['kmeans'], normalize="columns")
cost_matrix = crosstab.max().max() - crosstab.values
row_ind, col_ind = linear_sum_assignment(cost_matrix)
optimal_crosstab = crosstab.iloc[row_ind, :].iloc[:, col_ind]
sns.heatmap(optimal_crosstab, cmap='coolwarm', ax=axes[0], cbar=False)
axes[0].set_xlabel("Ground Truth Clusters", fontsize=14)
axes[0].set_ylabel("Predicted Clusters", fontsize=14)
axes[0].set_title("Normalized Confusion Matrix - kmeans (LD10)", fontsize=18)
axes[0].tick_params(axis='both', which='major', labelsize=10)

# Plot for LD10_D3
crosstab = pd.crosstab(adata_ld10_d3.obs['kmeans'], adata.obs['kmeans'], normalize="columns")
cost_matrix = crosstab.max().max() - crosstab.values
row_ind, col_ind = linear_sum_assignment(cost_matrix)
optimal_crosstab = crosstab.iloc[row_ind, :].iloc[:, col_ind]
sns.heatmap(optimal_crosstab, cmap='coolwarm', ax=axes[1], cbar=False)
axes[1].set_xlabel("Ground Truth Clusters", fontsize=14)
axes[1].set_ylabel("Predicted Clusters", fontsize=14)
axes[1].set_title("Normalized Confusion Matrix - kmeans (LD10_D3)", fontsize=18)
axes[1].tick_params(axis='both', which='major', labelsize=10)

# Plot for LD20
crosstab = pd.crosstab(adata_ld20.obs['kmeans'], adata.obs['kmeans'], normalize="columns")
cost_matrix = crosstab.max().max() - crosstab.values
row_ind, col_ind = linear_sum_assignment(cost_matrix)
optimal_crosstab = crosstab.iloc[row_ind, :].iloc[:, col_ind]
sns.heatmap(optimal_crosstab, cmap='coolwarm', ax=axes[2], cbar=False)
axes[2].set_xlabel("Ground Truth Clusters", fontsize=14)
axes[2].set_ylabel("Predicted Clusters", fontsize=14)
axes[2].set_title("Normalized Confusion Matrix - kmeans (LD20)", fontsize=18)
axes[2].tick_params(axis='both', which='major', labelsize=10)

plt.tight_layout()
plt.show()
