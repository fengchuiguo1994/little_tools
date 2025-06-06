SpatialGlueTest
```
import os
import torch
import pandas as pd
import scanpy as sc

import SpatialGlue
from SpatialGlue.preprocess import fix_seed
from SpatialGlue.preprocess import clr_normalize_each_cell, pca
from SpatialGlue.preprocess import construct_neighbor_graph

device = torch.device('cuda:2' if torch.cuda.is_available() else 'cpu')
os.environ['R_HOME'] = '/share/apps/software/R/3.6.1/bin/R'

file_fold = '/data/home/ruanlab/huangxingyu/Haoxi20230215/spatial/SpatialGlueTest/'
adata_omics1 = sc.read_h5ad(file_fold + 'adata_RNA.h5ad')
adata_omics2 = sc.read_h5ad(file_fold + 'adata_ADT.h5ad')

adata_omics1.var_names_make_unique()
adata_omics2.var_names_make_unique()

data_type = '10x'
random_seed = 2022
fix_seed(random_seed)


# RNA
sc.pp.filter_genes(adata_omics1, min_cells=10)
sc.pp.highly_variable_genes(adata_omics1, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata_omics1, target_sum=1e4)
sc.pp.log1p(adata_omics1)
sc.pp.scale(adata_omics1)
adata_omics1_high =  adata_omics1[:, adata_omics1.var['highly_variable']]
adata_omics1.obsm['feat'] = pca(adata_omics1_high, n_comps=adata_omics2.n_vars-1)

# Protein
adata_omics2 = clr_normalize_each_cell(adata_omics2)
sc.pp.scale(adata_omics2)
adata_omics2.obsm['feat'] = pca(adata_omics2, n_comps=adata_omics2.n_vars-1)


data = construct_neighbor_graph(adata_omics1, adata_omics2, datatype=data_type)
```