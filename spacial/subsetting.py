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
import pyucell as uc
import functions_spatial as func

adata_s2=sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample2_issues.h5ad")
adata_s1=sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260603_Sample1_all_annot.h5ad")


#subset myelo
myelo_cell_s1 = list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Xemnium_Sample1_MNP_ToEXPORT.csv", skiprows=5).columns)
myelo_cell_s2 = list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample2/Sample2_seqgeq_test.txt_MYELOID_For_Export.csv", skiprows=5).columns)


adata_myelo = adata_s2[adata_s2.obs['cell_id'].isin(myelo_cell_s2)].copy()
adata_myelo.obs['leiden_ICPC_full']=adata_myelo.obs['leiden_ICPC'].copy()
adata_myelo.obsm['umapICPC_full']=adata_myelo.obsm['umapICPC'].copy()
adata_myelo.X= adata_myelo.layers["counts"]
sc.pp.normalize_total(adata_myelo)
sc.pp.log1p(adata_myelo)
adata_myelo.layers["lognorm"] = adata_myelo.X.copy()
sc.pp.highly_variable_genes(adata_myelo, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata_myelo, zero_center=False)
adata_myelo.layers["scaled"] = adata_myelo.X.copy()

sc.pp.pca(adata_myelo, n_comps=50)
sc.pp.neighbors(adata_myelo, metric="cosine")
sc.tl.leiden(adata_myelo, flavor="igraph", n_iterations=-1, resolution=1)
sc.tl.umap(adata_myelo, min_dist=0.1)
#sc.pl.umap(adata_myelo, color="leiden")


func.ica(adata_myelo,n_components=50)
ics = pd.DataFrame(adata_myelo.varm['ICs'])
ics.index = adata_myelo.var.index
signature_v2 = pd.read_csv("/home/marine-louarn/ref/20251001  JDD_BreastK.csv")

for top in range(2,5):
    ic_tokeeps = func.ic_tokeep(adata_myelo,ics, signature_v2[signature_v2['Family']=='Myeloid'],"Cell_Subset",topX=top)
    for max_pc in range(10,50,10):
        adata_myelo.obsm["X_ic_pca" ] = np.concatenate((adata_myelo.obsm["X_ica" ][:,ic_tokeeps], adata_myelo.obsm["X_pca" ][:,list(range(0,max_pc))]),axis=1)
        adata_myelo.varm["IC_PCs"] = np.concatenate((adata_myelo.varm["ICs"][:,ic_tokeeps], adata_myelo.varm["PCs" ][:,list(range(0,max_pc))]),axis=1)

        sc.pp.neighbors(adata_myelo, metric="cosine",use_rep="X_ic_pca",key_added="neigh_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.leiden(adata_myelo, flavor="igraph", n_iterations=-1, resolution=1.5,neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc),key_added="leiden_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.umap(adata_myelo,key_added = "UMAP_IC"+str(top)+"_PC"+str(max_pc),neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc), min_dist=0.1)


sc.pl.embedding(adata_myelo, color="leiden_IC2_PC10",basis='UMAP_IC2_PC10')
list_tokeep= signature_v2['gene'].tolist()
func.to_seqgeq_myelo(adata_myelo,header_sg,"/home/marine-louarn/Documents/Xenium_Calderaro/20260616_Sample2_onlyMyelo_seqgeq.txt",list_tokeep)
func.to_seqgeq_myelo(adata_myelo,header_sg,"/home/marine-louarn/Documents/Xenium_Calderaro/20260616_Sample2_onlyMyelo_seqgeq_full.txt",adata_myelo.var_names.tolist())

adata_myelo.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260616_Sample2myelo_ICPC.h5ad")
adata_myelo_2= sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260616_Sample2myelo_ICPC.h5ad")
func.h5ad_to_tsv(adata_myelo_2,"/home/marine-louarn/Documents/Xenium_Calderaro/Xenium_exp1/Sample2_myelo","UMAP_IC2_PC10")

adata_myelo.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260528_Sample1myelo_ICPC.h5ad")
adata_myelo_1= sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260529_Sample1_myelo_ICPC.h5ad")
sc.pl.embedding(adata_myelo_1, color="leiden_IC2_PC10",basis='UMAP_IC2_PC10')
func.h5ad_to_tsv(adata_myelo_1,"/home/marine-louarn/Documents/Xenium_Calderaro/Xenium_exp1/Sample1_myelo","UMAP_IC2_PC10")



#add annot
Bcells= list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1_Myelo/Bcells.csv", skiprows=5).columns)
cd1c_DC= list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1_Myelo/CD1C+_DC.csv", skiprows=5).columns)
folr2= list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1_Myelo/FOLR2_Mac.csv", skiprows=5).columns)
il4i1= list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1_Myelo/IL4I1_Mac.csv", skiprows=5).columns)
kc= list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1_Myelo/KC.csv", skiprows=5).columns)
dc1= list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1_Myelo/Maybe_DC1.csv", skiprows=5).columns)
mono= list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1_Myelo/mono.csv", skiprows=5).columns)
pmn= list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1_Myelo/PMN.csv", skiprows=5).columns)
trem= list(pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1_Myelo/TREM2_Mac.csv", skiprows=5).columns)

adata_myelo.obs['MNP_cellType']= 'unassigned'
adata_myelo.obs.loc[adata_myelo.obs['cell_id'].isin(Bcells), 'MNP_cellType'] = 'Bcells'
adata_myelo.obs.loc[adata_myelo.obs['cell_id'].isin(cd1c_DC), 'MNP_cellType'] = 'CD1C+_DC'
adata_myelo.obs.loc[adata_myelo.obs['cell_id'].isin(folr2), 'MNP_cellType'] = 'FOLR2_Mac'
adata_myelo.obs.loc[adata_myelo.obs['cell_id'].isin(il4i1), 'MNP_cellType'] = 'IL4I1_Mac'
adata_myelo.obs.loc[adata_myelo.obs['cell_id'].isin(kc), 'MNP_cellType'] = 'KC'
adata_myelo.obs.loc[adata_myelo.obs['cell_id'].isin(dc1), 'MNP_cellType'] = 'DC1'
adata_myelo.obs.loc[adata_myelo.obs['cell_id'].isin(mono), 'MNP_cellType'] = 'Mono'
adata_myelo.obs.loc[adata_myelo.obs['cell_id'].isin(pmn), 'MNP_cellType'] = 'PMN'
adata_myelo.obs.loc[adata_myelo.obs['cell_id'].isin(trem), 'MNP_cellType'] = 'TREM2_Mac'
sc.pl.embedding(adata_myelo, color="MNP_cellType",basis='umapICPC')

sdata1 = sd.read_zarr("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1.zarr")
sdata1.tables["table"]=adata_myelo
sdata1.tables["table"].obs["region"] = "cell_boundaries"
sdata1.set_table_annotates_spatialelement("table", region="cell_boundaries")
leinden = sdata1.tables['table'].obs['MNP_cellType']
leinden.index= sdata1.tables['table'].obs.cell_id
sdata1['cell_boundaries']['MNP_cellType'] = sdata1.tables['table'].obs['MNP_cellType']
sdata1.pl.render_shapes("cell_boundaries", color="MNP_cellType").pl.show()
plt.show()

adata_myelo_1bis = adata_myelo_1[adata_myelo_1.obs['MNP_cellType']!='Bcells'].copy()
adata_myelo_1bis.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260629_Sample1_myelo_noB_ICPC.h5ad")

#Tcells
uc.compute_ucell_scores(adata_s2, signatures={'Tcell':['CD27','CD3E','CD3G','GZMA','GNLY','FCGR3A','NKG7','TRAC','CD7','CD3D']}, chunk_size=500)
sc.pl.embedding(adata_s2, color=["Tcell_UCell","leiden_ICPC"], cmap="viridis", size=20,basis='umapICPC')

adata_T = adata_s1[adata_s1.obs['leiden_ICPC']=='11'].copy()
adata_T.obs['leiden_ICPC_full']=adata_T.obs['leiden_ICPC'].copy()
adata_T.obsm['umapICPC_full']=adata_T.obsm['umapICPC'].copy()
adata_T.X= adata_T.layers["counts"]
sc.pp.normalize_total(adata_T)
sc.pp.log1p(adata_T)
adata_T.layers["lognorm"] = adata_T.X.copy()
sc.pp.highly_variable_genes(adata_T, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata_T, zero_center=False)
adata_T.layers["scaled"] = adata_T.X.copy()

sc.pp.pca(adata_T, n_comps=50)
sc.pp.neighbors(adata_T, metric="cosine")
sc.tl.leiden(adata_T, flavor="igraph", n_iterations=-1, resolution=1)
sc.tl.umap(adata_T, min_dist=0.1)


func.ica(adata_T,n_components=50)
ics = pd.DataFrame(adata_T.varm['ICs'])
ics.index = adata_T.var.index
signature_v2 = pd.read_csv("/home/marine-louarn/ref/20251001  JDD_BreastK.csv")

for top in range(2,5):
    ic_tokeeps = func.ic_tokeep(adata_T,ics, signature_v2[signature_v2['Family']=='Myeloid'],"Cell_Subset",topX=top)
    for max_pc in range(10,50,10):
        adata_T.obsm["X_ic_pca" ] = np.concatenate((adata_T.obsm["X_ica" ][:,ic_tokeeps], adata_T.obsm["X_pca" ][:,list(range(0,max_pc))]),axis=1)
        adata_T.varm["IC_PCs"] = np.concatenate((adata_T.varm["ICs"][:,ic_tokeeps], adata_T.varm["PCs" ][:,list(range(0,max_pc))]),axis=1)

        sc.pp.neighbors(adata_T, metric="cosine",use_rep="X_ic_pca",key_added="neigh_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.leiden(adata_T, flavor="igraph", n_iterations=-1, resolution=1.5,neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc),key_added="leiden_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.umap(adata_T,key_added = "UMAP_IC"+str(top)+"_PC"+str(max_pc),neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc), min_dist=0.1)


sc.pl.embedding(adata_T, color="leiden_IC2_PC10",basis='UMAP_IC2_PC10')

adata_T.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260721_Sample1-T_ICPC.h5ad")
list_tokeep= signature_v2['gene'].tolist()
func.to_seqgeq_myelo(adata_T,"/home/marine-louarn/Documents/Header_SeqGeq.txt","/home/marine-louarn/Documents/Xenium_Calderaro/20260721_Sample1_onlyT_seqgeq.txt",list_tokeep)
func.to_seqgeq_myelo(adata_T,"/home/marine-louarn/Documents/Header_SeqGeq.txt","/home/marine-louarn/Documents/Xenium_Calderaro/20260721_Sample1_onlyT_seqgeq_full.txt",adata_T.var_names.tolist())

#markers
sc.tl.rank_genes_groups(adata_myelo, 'leiden_ICPC', method='wilcoxon', key_added = "wilcoxon")
markers_myelo = sc.get.rank_genes_groups_df(adata_myelo,group=None,key='wilcoxon')
deg_mnp = pd.read_csv("/home/marine-louarn/ref/DEG_MNP_Fig1E.csv")

#pseudobulk
myelo_pseudobulk = sc.get.aggregate(adata_myelo, by=["leiden_ICPC"], func="sum", layer="counts")
ref_mnp = sc.read_h5ad("/home/marine-louarn/ref/2021_MNP_Verse.h5ad")
mnp_pseudobulk = sc.get.aggregate(ref_mnp, by=["MegaCluster"], func="sum", layer="counts")

genes_to_keep = list(set(deg_mnp['Gene']) & set(adata_myelo.var_names))
myelo_pseudobulk.X = myelo_pseudobulk.layers['sum']
myelo_pseudobulk_mat = myelo_pseudobulk.to_df().T
myelo_pseudobulk_mat = myelo_pseudobulk_mat[myelo_pseudobulk_mat.index.isin(genes_to_keep)]

mnp_pseudobulk.X = mnp_pseudobulk.layers['sum']
mnp_pseudobulk_mat = mnp_pseudobulk.to_df().T
mnp_pseudobulk_mat = mnp_pseudobulk_mat[mnp_pseudobulk_mat.index.isin(genes_to_keep)]

mat_both = pd.concat([mnp_pseudobulk_mat, myelo_pseudobulk_mat], axis=1)
cor = mat_both.corr()
cor = cor[mnp_pseudobulk_mat.columns]
cor = cor[cor.index.isin(myelo_pseudobulk_mat.columns)]
sns.heatmap(cor, annot=True)
plt.show()