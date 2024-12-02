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


############################################################################################################


### read in different conditions ###

# read in matrices for all conditions
ld1p1max = read_rds('/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_1p1max_parallelized.RDS')
ld1p5max = read_rds('/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_1p5max_parallelized.RDS')
ld2max = read_rds('/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_2max_parallelized.RDS')
ld5max = read_rds('/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_5max_parallelized.RDS')
ld5pct = read_rds('/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_5pct_parallelized.RDS')
ld10pct = read_rds('/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_10pct_parallelized.RDS')
ld20pct = read_rds('/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_20pct_parallelized.RDS')
ld50pct = read_rds('/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/sim_data/prox_assign/ld_mtx_50pct_parallelized.RDS')

# make an adata objects
adata_ld1p1max = sc.AnnData(ld1p1max.T)
adata_ld1p1max.obs = metadata.reset_index(drop=True)
adata_ld1p1max.var = pd.DataFrame(index=range(adata_ld1p1max.shape[1]))

adata_ld1p5max = sc.AnnData(ld1p5max.T)
adata_ld1p5max.obs = metadata.reset_index(drop=True)
adata_ld1p5max.var = pd.DataFrame(index=range(adata_ld1p5max.shape[1]))

adata_ld2max = sc.AnnData(ld2max.T)
adata_ld2max.obs = metadata.reset_index(drop=True)
adata_ld2max.var = pd.DataFrame(index=range(adata_ld2max.shape[1]))

adata_ld5max = sc.AnnData(ld5max.T)
adata_ld5max.obs = metadata.reset_index(drop=True)
adata_ld5max.var = pd.DataFrame(index=range(adata_ld5max.shape[1]))

adata_ld5pct = sc.AnnData(ld5pct.T)
adata_ld5pct.obs = metadata.reset_index(drop=True)
adata_ld5pct.var = pd.DataFrame(index=range(adata_ld5pct.shape[1]))

adata_ld10pct = sc.AnnData(ld10pct.T)
adata_ld10pct.obs = metadata.reset_index(drop=True)
adata_ld10pct.var = pd.DataFrame(index=range(adata_ld10pct.shape[1]))

adata_ld20pct = sc.AnnData(ld20pct.T)
adata_ld20pct.obs = metadata.reset_index(drop=True)
adata_ld20pct.var = pd.DataFrame(index=range(adata_ld20pct.shape[1]))

adata_ld50pct = sc.AnnData(ld50pct.T)
adata_ld50pct.obs = metadata.reset_index(drop=True)
adata_ld50pct.var = pd.DataFrame(index=range(adata_ld50pct.shape[1]))


############################################################################################################


np.sum(adata_ld1p1max.X.toarray() == adata_ld5max.X.toarray()) / (adata_ld1p1max.X.toarray().shape[0] * adata_ld1p1max.X.toarray().shape[1])


### normalize and cluster all conditions ###


# for loop to normalize and cluster all conditions
for adata in [adata_gt, adata_ld1p1max, adata_ld1p5max, adata_ld2max, adata_ld5max, adata_ld5pct, adata_ld10pct, adata_ld20pct, adata_ld50pct]:
    sc.pp.log1p(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)
    kmeans = KMeans(n_clusters=25, random_state=0).fit(adata.X)
    adata.obs['kmeans'] = kmeans.labels_
    print("one done")


############################################################################################################


### compare clustering results ###


# compare clustering results
ari_df = pd.DataFrame(columns=['ARI_leiden', "ARI_kmeans"], index=['GT', 'LD1p1max', 'LD1p5max', 'LD2max', 'LD5max', 'LD5pct', 'LD10pct', 'LD20pct', 'LD50pct'])
ari_df.loc['GT', 'ARI_leiden'] = metrics.adjusted_rand_score(adata_gt.obs['leiden'], adata_gt.obs['leiden'])
ari_df.loc['GT', 'ARI_kmeans'] = metrics.adjusted_rand_score(adata_gt.obs['kmeans'], adata_gt.obs['kmeans'])
ari_df.loc['LD1p1max', 'ARI_leiden'] = metrics.adjusted_rand_score(adata_gt.obs['leiden'], adata_ld1p1max.obs['leiden'])
ari_df.loc['LD1p1max', 'ARI_kmeans'] = metrics.adjusted_rand_score(adata_gt.obs['kmeans'], adata_ld1p1max.obs['kmeans'])
ari_df.loc['LD1p5max', 'ARI_leiden'] = metrics.adjusted_rand_score(adata_gt.obs['leiden'], adata_ld1p5max.obs['leiden'])
ari_df.loc['LD1p5max', 'ARI_kmeans'] = metrics.adjusted_rand_score(adata_gt.obs['kmeans'], adata_ld1p5max.obs['kmeans'])
ari_df.loc['LD2max', 'ARI_leiden'] = metrics.adjusted_rand_score(adata_gt.obs['leiden'], adata_ld2max.obs['leiden'])
ari_df.loc['LD2max', 'ARI_kmeans'] = metrics.adjusted_rand_score(adata_gt.obs['kmeans'], adata_ld2max.obs['kmeans'])
ari_df.loc['LD5max', 'ARI_leiden'] = metrics.adjusted_rand_score(adata_gt.obs['leiden'], adata_ld5max.obs['leiden'])
ari_df.loc['LD5max', 'ARI_kmeans'] = metrics.adjusted_rand_score(adata_gt.obs['kmeans'], adata_ld5max.obs['kmeans'])
ari_df.loc['LD5pct', 'ARI_leiden'] = metrics.adjusted_rand_score(adata_gt.obs['leiden'], adata_ld5pct.obs['leiden'])
ari_df.loc['LD5pct', 'ARI_kmeans'] = metrics.adjusted_rand_score(adata_gt.obs['kmeans'], adata_ld5pct.obs['kmeans'])
ari_df.loc['LD10pct', 'ARI_leiden'] = metrics.adjusted_rand_score(adata_gt.obs['leiden'], adata_ld10pct.obs['leiden'])
ari_df.loc['LD10pct', 'ARI_kmeans'] = metrics.adjusted_rand_score(adata_gt.obs['kmeans'], adata_ld10pct.obs['kmeans'])
ari_df.loc['LD20pct', 'ARI_leiden'] = metrics.adjusted_rand_score(adata_gt.obs['leiden'], adata_ld20pct.obs['leiden'])
ari_df.loc['LD20pct', 'ARI_kmeans'] = metrics.adjusted_rand_score(adata_gt.obs['kmeans'], adata_ld20pct.obs['kmeans'])
ari_df.loc['LD50pct', 'ARI_leiden'] = metrics.adjusted_rand_score(adata_gt.obs['leiden'], adata_ld50pct.obs['leiden'])
ari_df.loc['LD50pct', 'ARI_kmeans'] = metrics.adjusted_rand_score(adata_gt.obs['kmeans'], adata_ld50pct.obs['kmeans'])
print(ari_df)



import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment

# plot each confusion matrix
for adata in [adata_gt, adata_ld1p1max, adata_ld1p5max, adata_ld2max, adata_ld5max, adata_ld5pct, adata_ld10pct, adata_ld20pct, adata_ld50pct]:
    crosstab = pd.crosstab(adata_gt.obs['leiden'], adata.obs['leiden'], normalize="columns")
    cost_matrix = crosstab.max().max() - crosstab.values
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    optimal_crosstab = crosstab.iloc[row_ind, :].iloc[:, col_ind]
    sns.heatmap(optimal_crosstab, cmap='coolwarm')
    plt.xlabel("Ground Truth Clusters")
    plt.ylabel("Predicted Clusters")
    plt.title("Normalized Confusion Matrix - leiden")
    plt.show()



# Create subplots for side-by-side comparison
fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# Plot for LD10
crosstab = pd.crosstab(adata_ld1p1max.obs['leiden'], adata_gt.obs['leiden'], normalize="columns")
cost_matrix = crosstab.max().max() - crosstab.values
row_ind, col_ind = linear_sum_assignment(cost_matrix)
optimal_crosstab = crosstab.iloc[row_ind, :].iloc[:, col_ind]
sns.heatmap(optimal_crosstab, cmap='coolwarm', ax=axes[0], cbar=False)
axes[0].set_xlabel("Ground Truth Clusters", fontsize=14)
axes[0].set_ylabel("Predicted Clusters", fontsize=14)
axes[0].set_title("Normalized Confusion Matrix - leiden (LD10)", fontsize=18)
axes[0].tick_params(axis='both', which='major', labelsize=10)

