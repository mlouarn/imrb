import spatialdata as sd
import sopa
from spatialdata_io import xenium
from pathlib import Path
import shutil
import spatialdata_plot
import scanpy as sc
import harmonypy as hm
import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import anndata
from collections import defaultdict
import functions_spatial as func

adata_s2=sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample2_issues.h5ad")

endoth= list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample2/Sample2_seqgeq_test.txt_Endoth.csv", skiprows=5).columns)
epith= list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample2/Sample2_seqgeq_test.txt_Epithelial.csv", skiprows=5).columns)
fibro= list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample2/Sample2_seqgeq_test.txt_Fibrobl.csv", skiprows=5).columns)
myelo= list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample2/Sample2_seqgeq_test.txt_MYELOID_For_Export.csv", skiprows=5).columns)

adata_s2.obs['CellType']= 'unassigned'
adata_s2.obs.loc[adata_s2.obs['cell_id'].isin(endoth), 'CellType'] = 'Endothelial'
adata_s2.obs.loc[adata_s2.obs['cell_id'].isin(epith), 'CellType'] = 'Epithelial'
adata_s2.obs.loc[adata_s2.obs['cell_id'].isin(fibro), 'CellType'] = 'Fibroblast'
adata_s2.obs.loc[adata_s2.obs['cell_id'].isin(myelo), 'CellType'] = 'Myeloid'
sc.tl.rank_genes_groups(adata_s2, groupby="CellType", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata_s2, groupby="CellType", standard_scale="var", n_genes=5)
markers = sc.get.rank_genes_groups_df(adata_s2,group=None)

markers_R=pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample2/markers_handAnnotation.csv")
markers_R_p_val=markers_R[markers_R['p_val_adj']<0.05]
marker_other = markers_R_p_val[~markers_R_p_val['cluster'].isin(['myeloid','Unknown'])]['gene']
marker_myelo = markers_R_p_val[markers_R_p_val['cluster'].isin(['myeloid'])]['gene']
marker_to_exclude = marker_myelo[~marker_myelo.isin(list(set(marker_other)&set(marker_myelo)))]
genes_to_keep = adata_s2.var_names.copy()
genes_to_keep = [x for x in genes_to_keep if x not in list(marker_to_exclude)]

adata_myelo = adata_s2[adata_s2.obs['CellType']=='Myeloid'].copy()
adata_myelo.obs['leiden_ICPC_full']=adata_myelo.obs['leiden_ICPC'].copy()
adata_myelo.obsm['umapICPC_full']=adata_myelo.obsm['umapICPC'].copy()
adata_myelo.X= adata_myelo.layers["counts"]

adata_myelo_red=adata_myelo[:,genes_to_keep].copy()
sc.pp.normalize_total(adata_myelo_red)
sc.pp.log1p(adata_myelo_red)
adata_myelo_red.layers["lognorm"] = adata_myelo_red.X.copy()
sc.pp.highly_variable_genes(adata_myelo_red, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata_myelo_red, zero_center=False)
adata_myelo_red.layers["scaled"] = adata_myelo_red.X.copy()

sc.pp.pca(adata_myelo_red, n_comps=50)
sc.pp.neighbors(adata_myelo_red, metric="cosine")
sc.tl.leiden(adata_myelo_red, flavor="igraph", n_iterations=-1, resolution=1)
sc.tl.umap(adata_myelo_red, min_dist=0.1)
#sc.pl.umap(adata_myelo, color="leiden")


func.ica(adata_myelo_red,n_components=50)
ics = pd.DataFrame(adata_myelo_red.varm['ICs'])
ics.index = adata_myelo_red.var.index
signature_v2 = pd.read_csv("/home/marine-louarn/ref/20251001  JDD_BreastK.csv")

for top in range(2,5):
    ic_tokeeps = func.ic_tokeep(adata_myelo_red,ics, signature_v2[signature_v2['Family']=='Myeloid'],"Cell_Subset",topX=top)
    for max_pc in range(10,50,10):
        adata_myelo_red.obsm["X_ic_pca" ] = np.concatenate((adata_myelo_red.obsm["X_ica" ][:,ic_tokeeps], adata_myelo_red.obsm["X_pca" ][:,list(range(0,max_pc))]),axis=1)
        adata_myelo_red.varm["IC_PCs"] = np.concatenate((adata_myelo_red.varm["ICs"][:,ic_tokeeps], adata_myelo_red.varm["PCs" ][:,list(range(0,max_pc))]),axis=1)

        sc.pp.neighbors(adata_myelo_red, metric="cosine",use_rep="X_ic_pca",key_added="neigh_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.leiden(adata_myelo_red, flavor="igraph", n_iterations=-1, resolution=1.5,neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc),key_added="leiden_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.umap(adata_myelo_red,key_added = "UMAP_IC"+str(top)+"_PC"+str(max_pc),neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc), min_dist=0.1)


sc.pl.embedding(adata_myelo_red, color="leiden_IC2_PC10",basis='UMAP_IC2_PC10')
list_tokeep= signature_v2['gene'].tolist()

func.to_seqgeq_myelo(adata_myelo_red,header_sg,"/home/marine-louarn/Documents/Xenium_Calderaro/20260618_Sample2_onlyMyelo_DEG_seqgeq_full.txt",adata_myelo_red.var_names.tolist())

adata_myelo_2= sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260616_Sample2myelo_ICPC.h5ad")
adata_myelo_2.obsm['UMAP_wo_DEG'] = adata_myelo_red.obsm['UMAP_IC2_PC10'].copy()
adata_myelo_2.obs['leiden_wo_DEG'] = adata_myelo_red.obs['leiden_IC2_PC10'].copy()
func.to_seqgeq_myelo(adata_myelo_2,header_sg,"/home/marine-louarn/Documents/Xenium_Calderaro/20260618_Sample2_onlyMyelo_DEG_seqgeq_full.txt",adata_myelo_2.var_names.tolist())
adata_myelo_2.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260618_Sample2myelo_ICPC.h5ad")
