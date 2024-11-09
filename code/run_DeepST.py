

import os 
from deepst.DeepST import run
import matplotlib.pyplot as plt
from pathlib import Path
import scanpy as sc


# NOTE: using a different path here, will need the filtered_feature_bc_matrix.h5 for this
data_path = '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/code/DeepST/data/DLPFC/'
datadir = '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/data/'
data_name = '151673' #### project name
save_path = "../Results" #### save path
n_domains = 7 ###### the number of spatial domains.

deepen = run(save_path = save_path,
	task = "Identify_Domain", #### DeepST includes two tasks, one is "Identify_Domain" and the other is "Integration"
	pre_epochs = 50, ####  choose the number of training
	epochs = 100, #### choose the number of training
	use_gpu = True)
###### Read in 10x Visium data, or user can read in themselves.
adata = deepen._get_adata(platform="Visium", data_path=data_path, data_name=data_name)
# get adata from other folder for the gt
adata_4gt = sc.read_h5ad(f'{datadir}/151673.h5ad')
# set gt
adata.obs['sce.layer_guess'] = adata_4gt.obs['sce.layer_guess']


###### Segment the Morphological Image
adata = deepen._get_image_crop(adata, data_name=data_name) 

###### Data augmentation. spatial_type includes three kinds of "KDTree", "BallTree" and "LinearRegress", among which "LinearRegress"
###### is only applicable to 10x visium and the remaining omics selects the other two.
###### "use_morphological" defines whether to use morphological images.
adata = deepen._get_augment(adata, spatial_type="LinearRegress", use_morphological=True)

###### Build graphs. "distType" includes "KDTree", "BallTree", "kneighbors_graph", "Radius", etc., see adj.py
graph_dict = deepen._get_graph(adata.obsm["spatial"], distType = "BallTree")

###### Enhanced data preprocessing
data = deepen._data_process(adata, pca_n_comps = 200)

###### Training models
deepst_embed = deepen._fit(
		data = data,
		graph_dict = graph_dict,)
###### DeepST outputs
adata.obsm["DeepST_embed"] = deepst_embed

###### Define the number of space domains, and the model can also be customized. If it is a model custom priori = False.
adata = deepen._get_cluster_data(adata, n_domains=n_domains, priori = True)

###### Spatial localization map of the spatial domain
sc.pl.spatial(adata, color='DeepST_refine_domain', frameon = False, spot_size=150)

# print ARI
from sklearn import metrics
# make none a category in the gt bc there are some NA
mask = adata.obs['sce.layer_guess'][adata.obs['sce.layer_guess'].isna()].index
adata.obs['sce.layer_guess'] = adata.obs['sce.layer_guess'].cat.add_categories("None")
adata.obs['sce.layer_guess'].loc[mask] = "None"

ari = metrics.adjusted_rand_score(adata.obs['sce.layer_guess'], adata.obs['DeepST_refine_domain'])
print(ari)

sc.pl.spatial(adata, color=['sce.layer_guess'])
sc.pl.spatial(adata, color=['DeepST_refine_domain'])
