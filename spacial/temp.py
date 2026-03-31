# -*- coding: utf-8 -*-
"""
Spacial
"""

import scanpy as sc
import pandas as pd
import numpy as np
import h5py
import anndata as ad
import matplotlib.pyplot as plt
import warnings
import os
import squidpy as sq
import spatialdata as sd
import spatialdata_io as sio
import spatialdata_plot


# sdata = sd.read_zarr("/home/marine-louarn/Documents/Xenium_Calderaro/Xenium_Sample_Data_Shared-20260306T132506Z-1-001/Xenium_Sample_Data_Shared/transfer_8030963_files_3dd3e3a6_P736175_04/cell_feature_matrix.zarr")
adata = sc.read_10x_h5("/home/marine-louarn/Documents/Xenium_Calderaro/Xenium_Sample_Data_Shared-20260306T132506Z-1-001/Xenium_Sample_Data_Shared/transfer_8030963_files_3dd3e3a6_P736175_04/cell_feature_matrix.h5")
df = pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Xenium_Sample_Data_Shared-20260306T132506Z-1-001/Xenium_Sample_Data_Shared/transfer_8030963_files_3dd3e3a6_P736175_04/cells.csv")

df.set_index(adata.obs_names, inplace=True)
adata.obs = df.copy()
adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()

adata.layers["counts"] = adata.X.copy()

adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True)
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
) #no count_mt > logic with probes not looking at mt genes
sc.pp.filter_cells(adata, min_genes=10)
sc.pp.filter_genes(adata, min_cells=3)

sc.pp.normalize_total(adata, inplace=True) 
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.pl.umap(adata, color=["leiden"], wspace=0.4)

sq.pl.spatial_scatter(adata, library_id="spatial", shape=None, color=["leiden"], wspace=0.4)

sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
sq.gr.nhood_enrichment(adata, cluster_key="leiden")
sq.pl.nhood_enrichment(adata, cluster_key="leiden", title="Neighborhood enrichment adata")
