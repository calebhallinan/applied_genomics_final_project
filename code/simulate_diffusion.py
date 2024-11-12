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


### read in full transcriptome data ###

# Open the HDF5 file
file_path = f'{datadir}/mtx.h5'
with h5py.File(file_path, 'r') as h5file:
    # Access the datasets under 'mm10-1.2.0_premrna'
    data = h5file['mm10-1.2.0_premrna/data'][:]
    indices = h5file['mm10-1.2.0_premrna/indices'][:]
    indptr = h5file['mm10-1.2.0_premrna/indptr'][:]
    shape = h5file['mm10-1.2.0_premrna/shape'][:]

    # Create a sparse matrix using the loaded data
    sparse_matrix_full = csc_matrix((data, indices, indptr), shape=shape)

# Print the sparse matrix to verify it was loaded correctly
# print(sparse_matrix)

# make an adata object
adata_full = sc.AnnData(sparse_matrix_full.T)
# adata_gt.obs = metadata
adata_full.var = pd.DataFrame(index=range(adata_full.shape[1]))


############################################################################################################

# focus on just 191 genes rn

# normalize data
# sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50)

# perform lieiden clustering
sc.tl.leiden(adata, key_added="leiden", resolution=0.8)
# adata.obs['kmeans'] = kmeans.labels_.astype(str)

# plot results
sns.scatterplot(x='x', y='y', hue='leiden', data=adata.obs, size=0.1)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


############################################################################################################

### reread in data so it is not normalized ###

# set env and read in data
datadir = '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/'
# metadata_diff
metadata_diff = pd.read_csv("/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/sc_srt_meta.csv", index_col=0)
metadata_diff

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

# make an adata_diff object
adata_diff = sc.AnnData(sparse_matrix.T)
adata_diff.obs = metadata_diff.reset_index(drop=True)
adata_diff.var = pd.DataFrame(index=range(adata_diff.shape[1]))


# simulate diffusion by finding the nearest neighbors of each cell
# and subtracting cell i's count data and randomly adding it to one of its neighbors

# get nearest neighbors
from sklearn.neighbors import NearestNeighbors

# Assuming 'x' and 'y' are coordinates in adata_diff.obs
coords = adata_diff.obs[['x', 'y']].values

knn = NearestNeighbors(n_neighbors=5)
knn.fit(coords)
distances, indices = knn.kneighbors(coords)



ct = 0
# simulate diffusion
for i in range(adata_diff.X.shape[0]):
    # get neighbors
    neighbors = indices[i]
    # for loop 
    for r_neighbor in neighbors:
        # get 100 random genes to simulate diffusion
        for gene in range(100):
            # get random gene
            rand_gene = np.random.choice(adata_diff.X.shape[1])
            # check if gene is expressed in cell i
            if adata_diff.X[i, rand_gene] == 0:
                continue
            elif adata_diff.X[r_neighbor, rand_gene] > 2:
                # subtract gene count from cell i
                adata_diff.X[i, rand_gene] -= 1
                # add gene count to random neighbor
                adata_diff.X[r_neighbor, rand_gene] += 1
                ct +=1
            else:
                continue

print(ct)


# normalize data
sc.pp.log1p(adata_diff)

# perform lieiden clustering
sc.tl.pca(adata_diff, n_comps=50, svd_solver='arpack')
sc.pp.neighbors(adata_diff, n_neighbors=20, n_pcs=50) 
sc.tl.leiden(adata_diff, key_added="leiden", resolution=0.8)


# plot results
sns.scatterplot(x='x', y='y', hue='leiden', data=adata_diff.obs, size=0.1)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


# compare the two using ARI
ari = metrics.adjusted_rand_score(adata.obs['leiden'], adata_diff.obs['leiden'])
print(ari)


# Create the crosstab matrix
crosstab = pd.crosstab(adata.obs['leiden'], adata_diff.obs['leiden'])

# Convert the crosstab to a cost matrix by subtracting from the max value
cost_matrix = crosstab.max().max() - crosstab.values

# Use the Hungarian algorithm to find the optimal assignment
row_ind, col_ind = linear_sum_assignment(cost_matrix)

# Create a new crosstab with the optimal assignment
optimal_crosstab = crosstab.iloc[row_ind, col_ind]

optimal_crosstab

# generate a heatmap of clsuter assignments
sns.clustermap(optimal_crosstab, cmap='coolwarm', annot=True, cbar=False, figsize=(10,10), row_cluster=True, col_cluster=True)


sns.heatmap(optimal_crosstab, cmap='coolwarm')


############################################################################################################

### reread in data so it is not normalized ###

# set env and read in data
datadir = '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/'
# metadata_diff
metadata_diff = pd.read_csv("/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/sc_srt_meta.csv", index_col=0)
metadata_diff

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

# make an adata_diff object
adata_diff = sc.AnnData(sparse_matrix.T)
adata_diff.obs = metadata_diff.reset_index(drop=True)
adata_diff.var = pd.DataFrame(index=range(adata_diff.shape[1]))


# simulate diffusion by finding the nearest neighbors of each cell
# and subtracting cell i's count data and randomly adding it to one of its neighbors

# get nearest neighbors
from sklearn.neighbors import NearestNeighbors
from scipy.optimize import linear_sum_assignment

# Assuming 'x' and 'y' are coordinates in adata_diff.obs
coords = adata_diff.obs[['x', 'y']].values

knn = NearestNeighbors(n_neighbors=5)
knn.fit(coords)
distances, indices = knn.kneighbors(coords)



ct = 0
# simulate diffusion
for i in range(adata_diff.X.shape[0]):
    # get neighbors
    neighbors = indices[i]
    # for loop 
    for r_neighbor in neighbors:
        # get 100 random genes to simulate diffusion
        for gene in range(100):
            # get random gene
            rand_gene = np.random.choice(adata_diff.X.shape[1])
            # check if gene is expressed in cell i
            if adata_diff.X[i, rand_gene] == 0:
                continue
            elif adata_diff.X[r_neighbor, rand_gene] > 1:
                # subtract gene count from cell i
                adata_diff.X[i, rand_gene] -= 2
                # add gene count to random neighbor
                adata_diff.X[r_neighbor, rand_gene] += 2
                ct +=1
            else:
                continue

print(ct)


# normalize data
sc.pp.log1p(adata_diff)

# perform lieiden clustering
sc.tl.pca(adata_diff, n_comps=50, svd_solver='arpack')
sc.pp.neighbors(adata_diff, n_neighbors=20, n_pcs=50) 
sc.tl.leiden(adata_diff, key_added="leiden", resolution=0.8)


# plot results
sns.scatterplot(x='x', y='y', hue='leiden', data=adata_diff.obs, size=0.1)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


# compare the two using ARI
# remove NA
adata_diff.obs = adata_diff.obs.dropna()
ari = metrics.adjusted_rand_score(adata.obs['leiden'], adata_diff.obs['leiden'])
print(ari)



