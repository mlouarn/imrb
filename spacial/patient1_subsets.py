import spatialdata as sd
import sopa
from spatialdata_io import xenium
from pathlib import Path
import shutil
import spatialdata_plot
import scanpy as sc
import squidpy as sq
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

adata_p1= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20260918_Patient1_cellposeProseg.h5ad")
header_sg="/home/marine-louarn/Documents/Header_SeqGeq.txt"
signature_v2 = pd.read_csv("/home/marine-louarn/ref/20260928_JDD_BreastK_DC1update.csv")
top=2
max_pc=10

'''PC only:
tumour: #1 #3 #4 #5 #6 #14 #15 #16 

PC/IC:
endothelial: #11
Fibroblasts: #2
T cells: #8 #12
Myeloid: #7
B cells+plasma cells: #9 #10
'''
adata_p1.obs['Cell_Type']='0'
adata_p1.obs.loc[adata_p1.obs['leiden_ICPC'].isin(['1','3','4','5','6','14','15','16']), 'Cell_Type'] = 'tumour'
adata_p1.obs.loc[adata_p1.obs['leiden_ICPC'].isin(['0']), 'Cell_Type'] = 'epithelial'
adata_p1.obs.loc[adata_p1.obs['leiden_ICPC'].isin(['11']), 'Cell_Type'] = 'endothelial'
adata_p1.obs.loc[adata_p1.obs['leiden_ICPC'].isin(['2']), 'Cell_Type'] = 'Fibroblasts'
adata_p1.obs.loc[adata_p1.obs['leiden_ICPC'].isin(['8','12']), 'Cell_Type'] = 'T cells'
adata_p1.obs.loc[adata_p1.obs['leiden_ICPC'].isin(['7']), 'Cell_Type'] = 'Myeloid'
adata_p1.obs.loc[adata_p1.obs['leiden_ICPC'].isin(['9','10']), 'Cell_Type'] = 'B cells and plasma cells'
adata_p1.obs.loc[adata_p1.obs['leiden_ICPC'].isin(['13','17']), 'Cell_Type'] = 'unassigned'

sc.pl.embedding(adata_p1, color="Cell_Type",basis='umapICPC')

#tumour
adata_tumor = adata_p1[adata_p1.obs['leiden_ICPC'].isin(['1','3','4','5','6','14','15','16'])].copy()
adata_tumor.obs['leiden_ICPC_full']=adata_tumor.obs['leiden_ICPC'].copy()
adata_tumor.obsm['umapICPC_full']=adata_tumor.obsm['umapICPC'].copy()
adata_tumor.X= adata_tumor.layers["counts"]
sc.pp.normalize_total(adata_tumor)
sc.pp.log1p(adata_tumor)
adata_tumor.layers["lognorm"] = adata_tumor.X.copy()
sc.pp.highly_variable_genes(adata_tumor, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata_tumor, zero_center=False)
adata_tumor.layers["scaled"] = adata_tumor.X.copy()

sc.pp.pca(adata_tumor, n_comps=50)
sc.pp.neighbors(adata_tumor, metric="cosine")
sc.tl.leiden(adata_tumor, flavor="igraph", n_iterations=-1, resolution=1)
sc.tl.umap(adata_tumor, min_dist=0.1)

adata_tumor.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_tumour.h5ad")
func.to_seqgeq_myelo(adata_tumor,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_tumour.txt",adata_tumor.var_names)

#B & plasma
adata_B = adata_p1[adata_p1.obs['leiden_ICPC'].isin(['9','10'])].copy()
adata_B.obs['leiden_ICPC_full']=adata_B.obs['leiden_ICPC'].copy()
adata_B.obsm['umapICPC_full']=adata_B.obsm['umapICPC'].copy()
adata_B.X= adata_B.layers["counts"]
sc.pp.normalize_total(adata_B)
sc.pp.log1p(adata_B)
adata_B.layers["lognorm"] = adata_B.X.copy()
sc.pp.highly_variable_genes(adata_B, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata_B, zero_center=False)
adata_B.layers["scaled"] = adata_B.X.copy()

sc.pp.pca(adata_B, n_comps=50)
sc.pp.neighbors(adata_B, metric="cosine")
sc.tl.leiden(adata_B, flavor="igraph", n_iterations=-1, resolution=1)
sc.tl.umap(adata_B, min_dist=0.1)

adata_B.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_Bplasma.h5ad")
func.to_seqgeq_myelo(adata_B,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_Bplasma.txt",adata_B.var_names)


#endothelial
adata_endoth = adata_p1[adata_p1.obs['leiden_ICPC'].isin(['11'])].copy()
adata_endoth.obs['leiden_ICPC_full']=adata_endoth.obs['leiden_ICPC'].copy()
adata_endoth.obsm['umapICPC_full']=adata_endoth.obsm['umapICPC'].copy()
adata_endoth.X= adata_endoth.layers["counts"]
sc.pp.normalize_total(adata_endoth)
sc.pp.log1p(adata_endoth)
adata_endoth.layers["lognorm"] = adata_endoth.X.copy()
sc.pp.highly_variable_genes(adata_endoth, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata_endoth, zero_center=False)
adata_endoth.layers["scaled"] = adata_endoth.X.copy()

sc.pp.pca(adata_endoth, n_comps=50)
sc.pp.neighbors(adata_endoth, metric="cosine")
sc.tl.leiden(adata_endoth, flavor="igraph", n_iterations=-1, resolution=1)
sc.tl.umap(adata_endoth, min_dist=0.1)

func.ica(adata_endoth,n_components=50)
ics = pd.DataFrame(adata_endoth.varm['ICs'])
ics.index = adata_endoth.var.index
for top in range(2,6):
    ic_tokeeps = func.ic_tokeep(adata_endoth,ics, signature_v2[signature_v2['Family']=='Endothelial'],"Major_Cell_Populations",topX=top)
    for max_pc in range(10,60,10):
        adata_endoth.obsm["X_ic_pca" ] = np.concatenate((adata_endoth.obsm["X_ica" ][:,ic_tokeeps], adata_endoth.obsm["X_pca" ][:,list(range(0,max_pc))]),axis=1)
        adata_endoth.varm["IC_PCs"] = np.concatenate((adata_endoth.varm["ICs"][:,ic_tokeeps], adata_endoth.varm["PCs" ][:,list(range(0,max_pc))]),axis=1)
        sc.pp.neighbors(adata_endoth, metric="cosine",use_rep="X_ic_pca",key_added="neigh_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.leiden(adata_endoth, flavor="igraph", n_iterations=-1, resolution=1.5,neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc),key_added="leiden_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.umap(adata_endoth,key_added = "UMAP_IC"+str(top)+"_PC"+str(max_pc),neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc), min_dist=0.1)

adata_endoth.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_endothelial.h5ad")
func.to_seqgeq_myelo(adata_endoth,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_endothelial.txt",adata_endoth.var_names)

filtered_signature = signature_v2[signature_v2['Family']=='Endothelial']
filtered_signature = filtered_signature[filtered_signature['gene'].isin(adata_p1.var_names.tolist())]
filtered_signature = filtered_signature.dropna(how='any')

marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Major_Cell_Populations']].append(row['gene'])

marker_endoth=dict(marker_gene_dictionary)
uc.compute_ucell_scores(adata_endoth, signatures=marker_endoth, chunk_size=500)
list_sign_endoth = ["{}{}".format(i,"_UCell") for i in list(marker_endoth.keys())]
sc.pl.embedding(adata_endoth, color=list_sign_endoth,basis='UMAP_IC2_PC10')


#T
adata_T = adata_p1[adata_p1.obs['leiden_ICPC'].isin(['8','12'])].copy()
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

for top in range(2,6):
    ic_tokeeps = func.ic_tokeep(adata_T,ics, signature_v2[signature_v2['Major_Cell_Populations']=='T cell'],"Cell_Subset",topX=top)
    for max_pc in range(10,60,10):
        adata_T.obsm["X_ic_pca" ] = np.concatenate((adata_T.obsm["X_ica" ][:,ic_tokeeps], adata_T.obsm["X_pca" ][:,list(range(0,max_pc))]),axis=1)
        adata_T.varm["IC_PCs"] = np.concatenate((adata_T.varm["ICs"][:,ic_tokeeps], adata_T.varm["PCs" ][:,list(range(0,max_pc))]),axis=1)
        sc.pp.neighbors(adata_T, metric="cosine",use_rep="X_ic_pca",key_added="neigh_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.leiden(adata_T, flavor="igraph", n_iterations=-1, resolution=1.5,neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc),key_added="leiden_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.umap(adata_T,key_added = "UMAP_IC"+str(top)+"_PC"+str(max_pc),neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc), min_dist=0.1)

adata_T.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_T.h5ad")
func.to_seqgeq_myelo(adata_T,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_T.txt",adata_T.var_names)

filtered_signature = signature_v2[signature_v2['Major_Cell_Populations']=='T cell']
filtered_signature = filtered_signature[filtered_signature['gene'].isin(adata_p1.var_names.tolist())]
filtered_signature = filtered_signature.dropna(how='any')

marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Cell_Subset']].append(row['gene'])

marker_T=dict(marker_gene_dictionary)
uc.compute_ucell_scores(adata_T, signatures=marker_T, chunk_size=500)
list_sign_T = ["{}{}".format(i,"_UCell") for i in list(marker_T.keys())]
sc.pl.embedding(adata_T, color=list_sign_T,basis='UMAP_IC2_PC10')

#Fibroblast
adata_Fibro = adata_p1[adata_p1.obs['leiden_ICPC'].isin(['2'])].copy()
adata_Fibro.obs['leiden_ICPC_full']=adata_Fibro.obs['leiden_ICPC'].copy()
adata_Fibro.obsm['umapICPC_full']=adata_Fibro.obsm['umapICPC'].copy()
adata_Fibro.X= adata_Fibro.layers["counts"]
sc.pp.normalize_total(adata_Fibro)
sc.pp.log1p(adata_Fibro)
adata_Fibro.layers["lognorm"] = adata_Fibro.X.copy()
sc.pp.highly_variable_genes(adata_Fibro, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata_Fibro, zero_center=False)
adata_Fibro.layers["scaled"] = adata_Fibro.X.copy()

sc.pp.pca(adata_Fibro, n_comps=50)
sc.pp.neighbors(adata_Fibro, metric="cosine")
sc.tl.leiden(adata_Fibro, flavor="igraph", n_iterations=-1, resolution=1)
sc.tl.umap(adata_Fibro, min_dist=0.1)

func.ica(adata_Fibro,n_components=50)
ics = pd.DataFrame(adata_Fibro.varm['ICs'])
ics.index = adata_Fibro.var.index
for top in range(2,6):
    ic_tokeeps = func.ic_tokeep(adata_Fibro,ics, signature_v2[signature_v2['Family']=='Fibroblast'],"Major_Cell_Populations",topX=top)
    for max_pc in range(10,60,10):
        adata_Fibro.obsm["X_ic_pca" ] = np.concatenate((adata_Fibro.obsm["X_ica" ][:,ic_tokeeps], adata_Fibro.obsm["X_pca" ][:,list(range(0,max_pc))]),axis=1)
        adata_Fibro.varm["IC_PCs"] = np.concatenate((adata_Fibro.varm["ICs"][:,ic_tokeeps], adata_Fibro.varm["PCs" ][:,list(range(0,max_pc))]),axis=1)
        sc.pp.neighbors(adata_Fibro, metric="cosine",use_rep="X_ic_pca",key_added="neigh_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.leiden(adata_Fibro, flavor="igraph", n_iterations=-1, resolution=1.5,neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc),key_added="leiden_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.umap(adata_Fibro,key_added = "UMAP_IC"+str(top)+"_PC"+str(max_pc),neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc), min_dist=0.1)

adata_Fibro.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_fibroblast.h5ad")
func.to_seqgeq_myelo(adata_Fibro,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_fibroblast.txt",adata_Fibro.var_names)

filtered_signature = signature_v2[signature_v2['Family']=='Fibroblast']
filtered_signature = filtered_signature[filtered_signature['gene'].isin(adata_p1.var_names.tolist())]
filtered_signature = filtered_signature.dropna(how='any')

marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Cell_Subset']].append(row['gene'])

marker_fibro=dict(marker_gene_dictionary)
uc.compute_ucell_scores(adata_Fibro, signatures=marker_fibro, chunk_size=500)
list_sign_fibro = ["{}{}".format(i,"_UCell") for i in list(marker_fibro.keys())]
sc.pl.embedding(adata_Fibro, color=list_sign_fibro,basis='UMAP_IC2_PC10')

#Myeloid
adata_Myelo = adata_p1[adata_p1.obs['leiden_ICPC'].isin(['7'])].copy()
adata_Myelo.obs['leiden_ICPC_full']=adata_Myelo.obs['leiden_ICPC'].copy()
adata_Myelo.obsm['umapICPC_full']=adata_Myelo.obsm['umapICPC'].copy()
adata_Myelo.X= adata_Myelo.layers["counts"]
sc.pp.normalize_total(adata_Myelo)
sc.pp.log1p(adata_Myelo)
adata_Myelo.layers["lognorm"] = adata_Myelo.X.copy()
sc.pp.highly_variable_genes(adata_Myelo, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata_Myelo, zero_center=False)
adata_Myelo.layers["scaled"] = adata_Myelo.X.copy()

sc.pp.pca(adata_Myelo, n_comps=50)
sc.pp.neighbors(adata_Myelo, metric="cosine")
sc.tl.leiden(adata_Myelo, flavor="igraph", n_iterations=-1, resolution=1)
sc.tl.umap(adata_Myelo, min_dist=0.1)

func.ica(adata_Myelo,n_components=50)
ics = pd.DataFrame(adata_Myelo.varm['ICs'])
ics.index = adata_Myelo.var.index
for top in range(2,6):
    ic_tokeeps = func.ic_tokeep(adata_Myelo,ics, signature_v2[signature_v2['Family']=='Myeloid'],"Cell_Subset",topX=top)
    for max_pc in range(10,60,10):
        adata_Myelo.obsm["X_ic_pca" ] = np.concatenate((adata_Myelo.obsm["X_ica" ][:,ic_tokeeps], adata_Myelo.obsm["X_pca" ][:,list(range(0,max_pc))]),axis=1)
        adata_Myelo.varm["IC_PCs"] = np.concatenate((adata_Myelo.varm["ICs"][:,ic_tokeeps], adata_Myelo.varm["PCs" ][:,list(range(0,max_pc))]),axis=1)
        sc.pp.neighbors(adata_Myelo, metric="cosine",use_rep="X_ic_pca",key_added="neigh_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.leiden(adata_Myelo, flavor="igraph", n_iterations=-1, resolution=1.5,neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc),key_added="leiden_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.umap(adata_Myelo,key_added = "UMAP_IC"+str(top)+"_PC"+str(max_pc),neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc), min_dist=0.1)

adata_Myelo.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_myeloid_DC1update.h5ad")
func.to_seqgeq_myelo(adata_Myelo,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_myeloid_DC1update.txt",adata_Myelo.var_names)

filtered_signature = signature_v2[signature_v2['Family']=='Myeloid']
filtered_signature = filtered_signature[filtered_signature['gene'].isin(adata_p1.var_names.tolist())]
filtered_signature = filtered_signature.dropna(how='any')

marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Cell_Subset']].append(row['gene'])

marker_myelo=dict(marker_gene_dictionary)
uc.compute_ucell_scores(adata_Myelo, signatures=marker_myelo, chunk_size=500)
list_sign_myelo = ["{}{}".format(i,"_UCell") for i in list(marker_myelo.keys())]
#sc.pl.embedding(adata_Myelo, color=list_sign_myelo,basis='spatial')
sc.pl.embedding(adata_Myelo, color=list_sign_myelo,basis='UMAP_IC2_PC10')

#readall + Markers
adata_Myelo= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_myeloid.h5ad")
adata_Fibro= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_fibroblast.h5ad")
adata_T= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_T.h5ad")
adata_endoth= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_endothelial.h5ad")
adata_B= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_Bplasma.h5ad")
adata_tumor= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_tumour.h5ad")

adata=[adata_Myelo,adata_Fibro, adata_T, adata_endoth,adata_B, adata_tumor]
name=['Myelo','Fibro', 'T', 'endoth','B','Tumour']
for i in range(0,len(adata)):
    sc.tl.rank_genes_groups(adata[i], 'leiden', method='wilcoxon', key_added = "wilcoxon")
    markers = sc.get.rank_genes_groups_df(adata[i],group=None,key='wilcoxon')
    #sc.pl.rank_genes_groups_heatmap(adata[i], n_genes=5, key="wilcoxon", groupby="leiden", show_gene_labels=True)
    #sc.pl.rank_genes_groups_dotplot(adata[i],n_genes=5, key="wilcoxon", groupby="leiden")
    markers.to_csv('/home/marine-louarn/Documents/Lung_Spatial/Patient1/20261007_marker_'+name[i]+'_PC.csv')

adata=[adata_Myelo,adata_Fibro, adata_T, adata_endoth]
name=['Myelo','Fibro', 'T', 'endoth']
for i in range(0,len(adata)):
    sc.tl.rank_genes_groups(adata[i], 'leiden_IC2_PC10', method='wilcoxon', key_added = "wilcoxon")
    markers = sc.get.rank_genes_groups_df(adata[i],group=None,key='wilcoxon')
    sc.pl.rank_genes_groups_heatmap(adata[i], n_genes=5, key="wilcoxon", groupby="leiden_IC2_PC10", show_gene_labels=True)
    sc.pl.rank_genes_groups_dotplot(adata[i],n_genes=5, key="wilcoxon", groupby="leiden_IC2_PC10")
    markers.to_csv('/home/marine-louarn/Documents/Lung_Spatial/Patient1/20261007_marker_'+name[i]+'_IC2_PC10.csv')


#annotations
adata_Myelo= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_myeloid.h5ad")
dc1= list(pd.read_csv("/home/marine-louarn/Documents/Lung_Spatial/Patient1/Myeloid/Export_leidenICPC#5_DC1.csv", skiprows=5).columns)
folr2= list(pd.read_csv("/home/marine-louarn/Documents/Lung_Spatial/Patient1/Myeloid/Export_leiden_new_ICPC_#15inPC#8_Mac.FOLR2.csv", skiprows=5).columns)
marco= list(pd.read_csv("/home/marine-louarn/Documents/Lung_Spatial/Patient1/Myeloid/Export_leiden_new_ICPC_#15inPC#8_Mac.MARCO.csv", skiprows=5).columns)

adata_Myelo.obs['Cell_population']= 'MNP.unassigned'
adata_Myelo.obs.loc[adata_Myelo.obs['leiden'].isin(['3']), 'Cell_population'] = 'Mac.Fibro.Conta'
adata_Myelo.obs.loc[adata_Myelo.obs['leiden'].isin(['14']), 'Cell_population'] = 'MMNP.PLVAP.vascular'
adata_Myelo.obs.loc[adata_Myelo.obs['leiden'].isin(['13']), 'Cell_population'] = 'MNP.Bcell.conta'
adata_Myelo.obs.loc[adata_Myelo.obs['leiden'].isin(['12']), 'Cell_population'] = 'Mac.CSF3R'
adata_Myelo.obs.loc[adata_Myelo.obs['leiden'].isin(['11']), 'Cell_population'] = 'MNP.Epithelial.conta'
adata_Myelo.obs.loc[adata_Myelo.obs['leiden'].isin(['10']), 'Cell_population'] = 'Mac.IL4I1'
adata_Myelo.obs.loc[adata_Myelo.obs['leiden'].isin(['7']), 'Cell_population'] = 'Prolif.MNP'
adata_Myelo.obs.loc[adata_Myelo.obs['leiden'].isin(['5']), 'Cell_population'] = 'cMo'
adata_Myelo.obs.loc[adata_Myelo.obs['leiden'].isin(['4']), 'Cell_population'] = 'IL1B_Mono'
adata_Myelo.obs.loc[adata_Myelo.obs['leiden'].isin(['1']), 'Cell_population'] = 'Mac.TREM2'
adata_Myelo.obs.loc[adata_Myelo.obs['leiden'].isin(['0']), 'Cell_population'] = 'Mo.Infl.LYZ'
adata_Myelo.obs.loc[adata_Myelo.obs['leiden_IC2_PC10'].isin(['12']), 'Cell_population'] = 'T.MNP.conta'
adata_Myelo.obs.loc[adata_Myelo.obs['leiden_IC2_PC10'].isin(['18']), 'Cell_population'] = 'mregDC'
adata_Myelo.obs.loc[adata_Myelo.obs['leiden_IC2_PC10'].isin(['9']), 'Cell_population'] = 'DC2.3'
adata_Myelo.obs.loc[adata_Myelo.obs['cell_id'].isin(dc1), 'Cell_population'] = 'DC1'
adata_Myelo.obs.loc[adata_Myelo.obs['cell_id'].isin(folr2), 'Cell_population'] = 'Mac.FOLR2'
adata_Myelo.obs.loc[adata_Myelo.obs['cell_id'].isin(marco), 'Cell_population'] = 'Mac.MARCO'

adata_T= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_T.h5ad")
adata_T.obs['Cell_population']= 'T.unassigned'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['0']), 'Cell_population'] = 'Treg.CCR8'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['1']), 'Cell_population'] = 'T_stroma'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['2']), 'Cell_population'] = 'T.EGR1'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['3']), 'Cell_population'] = 'T.conta.macro'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['4']), 'Cell_population'] = 'T_B'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['5']), 'Cell_population'] = 'T_DC'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['7']), 'Cell_population'] = 'T_endoth'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['8']), 'Cell_population'] = 'T.conta.tumour'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['9']), 'Cell_population'] = 'T.CD8'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['10']), 'Cell_population'] = 'T_epith'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['11']), 'Cell_population'] = 'T.CCL1'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['12']), 'Cell_population'] = 'T.prolif'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['13']), 'Cell_population'] = 'T.conta.Myelo'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['14']), 'Cell_population'] = 'T_stroma'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['15']), 'Cell_population'] = 'T_mesenchym'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['16']), 'Cell_population'] = 'T_vasculaire'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['17']), 'Cell_population'] = 'T.IL22'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['18']), 'Cell_population'] = 'T.helper'
adata_T.obs.loc[adata_T.obs['leiden'].isin(['19']), 'Cell_population'] = 'T.conta.B'

adata_B= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_Bplasma.h5ad")
adata_B.obs['Cell_population']= 'B.unassigned'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['0']), 'Cell_population'] = 'B_stroma'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['1']), 'Cell_population'] = 'Plasma'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['2']), 'Cell_population'] = 'B.conta.T'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['3']), 'Cell_population'] = 'B.CD83'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['4']), 'Cell_population'] = 'Plasma.conta.T'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['5']), 'Cell_population'] = 'B_endo'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['6']), 'Cell_population'] = 'B.CD83'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['7']), 'Cell_population'] = 'Plasma'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['8']), 'Cell_population'] = 'B_macro'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['9']), 'Cell_population'] = 'B.prolf'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['10']), 'Cell_population'] = 'B.CD83'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['11']), 'Cell_population'] = 'B_stroma'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['12']), 'Cell_population'] = 'B.conta.Myelo'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['13']), 'Cell_population'] = 'B_perivascular'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['15']), 'Cell_population'] = 'B.conta.epith'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['16']), 'Cell_population'] = 'Plasma'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['17']), 'Cell_population'] = 'B_Plasma'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['18']), 'Cell_population'] = 'B_Plasma'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['19']), 'Cell_population'] = 'Plasma.CXCL10'
adata_B.obs.loc[adata_B.obs['leiden'].isin(['20']), 'Cell_population'] = 'B.CD83'

adata_Fibro= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_fibroblast.h5ad")
adata_Fibro.obs['Cell_population']= 'Fibro.unassigned'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['0']), 'Cell_population'] = 'Fibro.conta.B'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['1']), 'Cell_population'] = 'Fibro.conta.T'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['2']), 'Cell_population'] = 'Fibro.iCAF'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['3']), 'Cell_population'] = 'Fibro.EGR1'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['4']), 'Cell_population'] = 'Fibro.CAF'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['5']), 'Cell_population'] = 'Fibro_stroma'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['6']), 'Cell_population'] = 'Fibro.conta.Epith'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['7']), 'Cell_population'] = 'Fibro.COL1A1'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['8']), 'Cell_population'] = 'Fibro.iCAF-like'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['9']), 'Cell_population'] = 'Fibro.myCAF'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['10']), 'Cell_population'] = 'Fibro.remodeling.CAF'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['11']), 'Cell_population'] = 'Fibro.CAF'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['12']), 'Cell_population'] = 'Fibro.HLA-II'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['13']), 'Cell_population'] = 'Fibro.CAF'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['14']), 'Cell_population'] = 'Fibro.prolif'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['15']), 'Cell_population'] = 'Fibro.FBLN1'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['16']), 'Cell_population'] = 'Fibro.CAF.IL6'
adata_Fibro.obs.loc[adata_Fibro.obs['leiden_IC2_PC10'].isin(['17']), 'Cell_population'] = 'Fibro.myCAF'

adata_endoth= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_endothelial.h5ad")
adata_endoth.obs['Cell_population']= 'Endoth.unassigned'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['0']), 'Cell_population'] = 'Endoth.PROX1'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['1']), 'Cell_population'] = 'Endoth.conta.fibro'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['2']), 'Cell_population'] = 'Endoth.HLA'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['3']), 'Cell_population'] = 'Endoth.conta.B'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['4']), 'Cell_population'] = 'Endoth.ANGPT2'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['5']), 'Cell_population'] = 'Endoth.PLVAP'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['7']), 'Cell_population'] = 'Endoth.vascular'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['8']), 'Cell_population'] = 'Endoth.IFN'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['9']), 'Cell_population'] = 'Endoth.conta.macro'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['10']), 'Cell_population'] = 'Endoth.conta.fribro'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['11']), 'Cell_population'] = 'Endoth.conta.T'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['12']), 'Cell_population'] = 'Endoth.ACKR3'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['13']), 'Cell_population'] = 'Endoth.perivascular'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['14']), 'Cell_population'] = 'Endoth.ACKR3'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['15']), 'Cell_population'] = 'Endoth.conta.epith'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['16']), 'Cell_population'] = 'Endoth.perivascular'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['17']), 'Cell_population'] = 'Endoth.conta.myelo'
adata_endoth.obs.loc[adata_endoth.obs['leiden_IC2_PC10'].isin(['18']), 'Cell_population'] = 'Endoth.blood'

adata_tumor= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20260928_Patient1_tumour.h5ad")
adata_tumor.obs['Cell_population']= 'Endoth.unassigned'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['0']), 'Cell_population'] = 'Tumour.conta.Myelo'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['1']), 'Cell_population'] = 'Tumour.INF'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['2']), 'Cell_population'] = 'Tumour.MUC1'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['3']), 'Cell_population'] = 'Tumour.prolif'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['4']), 'Cell_population'] = 'Tumour.INF'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['5']), 'Cell_population'] = 'Tumour.plasticity'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['6']), 'Cell_population'] = 'Tumour.basal'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['7']), 'Cell_population'] = 'Tumour.prolif'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['8']), 'Cell_population'] = 'Tumour_endoth'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['9']), 'Cell_population'] = 'Tumour.EPCAM'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['10']), 'Cell_population'] = 'Tumour.EMT'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['11']), 'Cell_population'] = 'Tumour.pEMT'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['12']), 'Cell_population'] = 'Tumour.INF'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['13']), 'Cell_population'] = 'Tumour_T'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['14']), 'Cell_population'] = 'Tumour_B'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['15']), 'Cell_population'] = 'Tumour.CCL4'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['16']), 'Cell_population'] = 'Tumour.hypoxic'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['17']), 'Cell_population'] = 'Tumour.stress'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['18']), 'Cell_population'] = 'Tumour_CAF'
adata_tumor.obs.loc[adata_tumor.obs['leiden'].isin(['19']), 'Cell_population'] = 'Tumour.basal'

adata_Myelo.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20261008_Patient1_myeloid.h5ad")
adata_Fibro.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20261008_Patient1_fibroblast.h5ad")
adata_T.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20261008_Patient1_T.h5ad")
adata_endoth.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20261008_Patient1_endothelial.h5ad")
adata_B.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20261008_Patient1_Bplasma.h5ad")
adata_tumor.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20261008_Patient1_tumour.h5ad")

#remap to full object

adata_p1= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20260918_Patient1_cellposeProseg.h5ad")
adata_p1.obs.loc[adata_p1.obs['leiden_ICPC'].isin(['0']), 'Cell_Population'] = 'epithelial'
adata_list=[adata_Myelo,adata_Fibro, adata_T, adata_endoth,adata_B,adata_tumor]
for adata in adata_list:
    adata_p1.obs.loc[adata_p1.obs['cell_id'].isin(adata.obs['cell_id']), 'Cell_Population'] = adata.obs['Cell_population'].astype(str)
    adata_p1.obs.loc[adata_p1.obs['cell_id'].isin(adata.obs['cell_id']), 'Subset_leiden'] = adata.obs['leiden'].astype(str)
    adata_p1.obs.loc[adata_p1.obs['cell_id'].isin(adata.obs['cell_id']), 'Subset_umap_1'] = adata.obsm['X_umap'][:,0]
    adata_p1.obs.loc[adata_p1.obs['cell_id'].isin(adata.obs['cell_id']), 'Subset_umap_2'] = adata.obsm['X_umap'][:,1]
    if 'leiden_IC2_PC10' in adata.obs.columns :
        adata_p1.obs.loc[adata_p1.obs['cell_id'].isin(adata.obs['cell_id']), 'Subset_leiden_IC_PC'] = adata.obs['leiden_IC2_PC10'].astype(str)
        adata_p1.obs.loc[adata_p1.obs['cell_id'].isin(adata.obs['cell_id']), 'Subset_umap_ICPC_1'] = adata.obsm['UMAP_IC2_PC10'][:,0]
        adata_p1.obs.loc[adata_p1.obs['cell_id'].isin(adata.obs['cell_id']), 'Subset_umap_ICPC_2'] = adata.obsm['UMAP_IC2_PC10'][:,1]

adata_p1.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20261008_Patient1_withAllSubtypes.h5ad")
