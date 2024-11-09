import scanpy as sc
import torch
import matplotlib.pyplot as plt
import imageio
from numpy import newaxis
import os

from proust.Train import *
from proust.cluster import *
from proust.prep import *


# reference: https://github.com/JianingYao/proust/blob/master/tutorial_IFDLPFC.ipynb


##########################################################################################

### set env and read in data ###

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
os.environ['R_HOME'] = '/Library/Frameworks/R.framework/Resources' #replace with your R_HOME path
# set seed and number of clusters
seed = 1998
n_clusters = '7' #predefined number of clusters

# read in data
datadir = '/Users/calebhallinan/Desktop/jhu/classes/applied_genomics_final_proj/data/'
adata = sc.read_h5ad(f'{datadir}/151673.h5ad')
adata.var_names_make_unique()

# get image
image = adata.uns['spatial']["151673"]['images']['hires']
image.shape
# image = image[1:6, :, :]
# set pixels
x_pixel = adata.obsm['spatial'][:, 1].astype(int)
y_pixel = adata.obsm['spatial'][:, 0].astype(int)


##########################################################################################


# Extract img features
Img_learn(adata, image, device=device)
print("Finish extracting image features!")
# Reconstruct protein and gene features
model = proust(adata, device=device, random_seed=seed)
adata = model.train()