import scanpy as sc
from sklearn import metrics

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

datadir = '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/data/'
memory=[1,2,3,4,5,6,7,8,9,10]
during_time=[1,2,3,4,5,6,7,8,9,10]
ari_list=[1,2,3,4,5,6,7,8,9,10]

new_data = sc.read_h5ad(f'{datadir}/151673.h5ad')


#Run 10 repeated experiments
for i in range(10):
    print(i)

    import tracemalloc
    import time
    from sklearn.cluster import KMeans
 
    tracemalloc.start()
    start_time=time.time()
 
    adata = sc.read_h5ad(f'{datadir}/151673.h5ad')

    
    n_clusters=7    #according own dataset
    sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50) # 20

    # perform kmeans clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(adata.obsm["X_pca"])
    adata.obs['kmeans'] = kmeans.labels_.astype(str)

    obs_df = adata.obs.dropna()
    ari = metrics.adjusted_rand_score(obs_df['sce.layer_guess'], obs_df['kmeans'])
    ari_list[i]=ari

    end_time=time.time()
    during=end_time-start_time
  
    size, peak = tracemalloc.get_traced_memory()

    tracemalloc.stop()
    
    memory[i]=peak /1024/1024
    during_time[i]=during
    
   
    print('memory blocks peak:{:>10.4f} MB'.format(memory[i]))
    print('time: {:.4f} s'.format(during_time[i]))
    print('ARI:{}'.format(ari_list[i]))

    new_data.obs['pred_{}'.format(i+1)]=adata.obs['kmeans']


# plot the results
sc.pl.spatial(new_data, color=['pred_10'])
sc.pl.spatial(new_data, color=['sce.layer_guess'])
