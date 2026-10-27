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

adata_p2= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20260922_Patient2_cellposeProseg.h5ad")
header_sg="/home/marine-louarn/Documents/Header_SeqGeq.txt"
signature_v2 = pd.read_csv("/home/marine-louarn/ref/20260928_JDD_BreastK_DC1update.csv")
top=2
max_pc=10

'''PC only:
tumour: 38 37 34 33 31 30 29 27 25 24 23 22 21 16 15 13 12 11 4 & part of 19
B cells: #2
Plasma : 3
Mast : 0
pDC : 19

PC/IC:
endothelial: #17 
Fibroblasts: #40 39 32 18 9 7 5 
T cells: 26 20 14 8 6 
Myeloid: #1 36

'''
adata_p2.obs['Cell_Type']='unasigned'
adata_p2.obs.loc[adata_p2.obs['leiden_ICPC'].isin(['38','37','34','33','31','30','29','27','25','24','23','22','21','16','15','13','12','11','4','19']), 'Cell_Type'] = 'tumour'
adata_p2.obs.loc[adata_p2.obs['leiden_ICPC'].isin(['17']), 'Cell_Type'] = 'endothelial'
adata_p2.obs.loc[adata_p2.obs['leiden_ICPC'].isin(['40','39','32','18','9','7','5']), 'Cell_Type'] = 'Fibroblasts'
adata_p2.obs.loc[adata_p2.obs['leiden_ICPC'].isin(['26','20','14','8','6']), 'Cell_Type'] = 'T cells'
adata_p2.obs.loc[adata_p2.obs['leiden_ICPC'].isin(['1','36']), 'Cell_Type'] = 'Myeloid'
adata_p2.obs.loc[adata_p2.obs['leiden_ICPC'].isin(['2']), 'Cell_Type'] = 'B cells'
adata_p2.obs.loc[adata_p2.obs['leiden_ICPC'].isin(['3']), 'Cell_Type'] = 'Plasma cells'
pDC= list(pd.read_csv("Patient2/Lung_Merscope_P2_leiden_ICPC_19_U_pDC.csv", skiprows=5).columns)
adata_p2.obs.loc[adata_p2.obs['cell_id'].isin(pDC), 'Cell_Type'] = 'pDC'
adata_p2.obs.loc[adata_p2.obs['leiden_ICPC'].isin(['0']), 'Cell_Type'] = 'Mast Cells'
adata_p2.obs.loc[adata_p2.obs['leiden_ICPC'].isin(['35','28','10']), 'Cell_Type'] = 'unsure'
adata_p2.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20260922_Patient2_cellposeProseg.h5ad")

sc.pl.embedding(adata_p2, color="Cell_Type",basis='umapICPC')

#unsure
adata_unsure = adata_p2[adata_p2.obs['leiden_ICPC'].isin(['35','28','10'])].copy()
adata_unsure.obs['leiden_ICPC_full']=adata_unsure.obs['leiden_ICPC'].copy()
adata_unsure.obsm['umapICPC_full']=adata_unsure.obsm['umapICPC'].copy()
adata_unsure.X= adata_unsure.layers["counts"]
sc.pp.normalize_total(adata_unsure)
sc.pp.log1p(adata_unsure)
adata_unsure.layers["lognorm"] = adata_unsure.X.copy()
sc.pp.highly_variable_genes(adata_unsure, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata_unsure, zero_center=False)
adata_unsure.layers["scaled"] = adata_unsure.X.copy()

sc.pp.pca(adata_unsure, n_comps=50)
sc.pp.neighbors(adata_unsure, metric="cosine")
sc.tl.leiden(adata_unsure, flavor="igraph", n_iterations=-1, resolution=1)
sc.tl.umap(adata_unsure, min_dist=0.1)

adata_unsure.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_unsure.h5ad")
func.to_seqgeq_myelo(adata_unsure,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_unsure.txt",adata_unsure.var_names)


#tumour
adata_tumor = adata_p2[adata_p2.obs['leiden_ICPC'].isin(['38','37','34','33','31','30','29','27','25','24','23','22','21','16','15','13','12','11','4','19'])].copy()
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

adata_tumor.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_tumour.h5ad")
func.to_seqgeq_myelo(adata_tumor,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_tumour.txt",adata_tumor.var_names)

#B
adata_B = adata_p2[adata_p2.obs['leiden_ICPC'].isin(['3'])].copy()
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

adata_B.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_B.h5ad")
func.to_seqgeq_myelo(adata_B,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_B.txt",adata_B.var_names)

#Plasma
adata_plasma = adata_p2[adata_p2.obs['leiden_ICPC'].isin(['2'])].copy()
adata_plasma.obs['leiden_ICPC_full']=adata_plasma.obs['leiden_ICPC'].copy()
adata_plasma.obsm['umapICPC_full']=adata_plasma.obsm['umapICPC'].copy()
adata_plasma.X= adata_plasma.layers["counts"]
sc.pp.normalize_total(adata_plasma)
sc.pp.log1p(adata_plasma)
adata_plasma.layers["lognorm"] = adata_plasma.X.copy()
sc.pp.highly_variable_genes(adata_plasma, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata_plasma, zero_center=False)
adata_plasma.layers["scaled"] = adata_plasma.X.copy()

sc.pp.pca(adata_plasma, n_comps=50)
sc.pp.neighbors(adata_plasma, metric="cosine")
sc.tl.leiden(adata_plasma, flavor="igraph", n_iterations=-1, resolution=1)
sc.tl.umap(adata_plasma, min_dist=0.1)

adata_plasma.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_plasma.h5ad")
func.to_seqgeq_myelo(adata_plasma,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_plasma.txt",adata_plasma.var_names)

#endothelial
adata_endoth = adata_p2[adata_p2.obs['leiden_ICPC'].isin(['17'])].copy()
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

adata_endoth.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_endothelial.h5ad")
func.to_seqgeq_myelo(adata_endoth,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_endothelial.txt",adata_endoth.var_names)

filtered_signature = signature_v2[signature_v2['Family']=='Endothelial']
filtered_signature = filtered_signature[filtered_signature['gene'].isin(adata_p2.var_names.tolist())]
filtered_signature = filtered_signature.dropna(how='any')

marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Major_Cell_Populations']].append(row['gene'])

marker_endoth=dict(marker_gene_dictionary)
uc.compute_ucell_scores(adata_endoth, signatures=marker_endoth, chunk_size=500)
list_sign_endoth = ["{}{}".format(i,"_UCell") for i in list(marker_endoth.keys())]
sc.pl.embedding(adata_endoth, color=list_sign_endoth,basis='UMAP_IC2_PC10')


#T
adata_T = adata_p2[adata_p2.obs['leiden_ICPC'].isin(['26','20','14','8','6'])].copy()
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

adata_T.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_T.h5ad")
func.to_seqgeq_myelo(adata_T,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_T.txt",adata_T.var_names)

filtered_signature = signature_v2[signature_v2['Major_Cell_Populations']=='T cell']
filtered_signature = filtered_signature[filtered_signature['gene'].isin(adata_p2.var_names.tolist())]
filtered_signature = filtered_signature.dropna(how='any')

marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Cell_Subset']].append(row['gene'])

marker_T=dict(marker_gene_dictionary)
uc.compute_ucell_scores(adata_T, signatures=marker_T, chunk_size=500)
list_sign_T = ["{}{}".format(i,"_UCell") for i in list(marker_T.keys())]
sc.pl.embedding(adata_T, color=list_sign_T,basis='UMAP_IC2_PC10')

#Fibroblast
adata_Fibro = adata_p2[adata_p2.obs['leiden_ICPC'].isin(['40','39','32','18','9','7','5'])].copy()
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

adata_Fibro.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_fibroblast.h5ad")
func.to_seqgeq_myelo(adata_Fibro,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_fibroblast.txt",adata_Fibro.var_names)

filtered_signature = signature_v2[signature_v2['Family']=='Fibroblast']
filtered_signature = filtered_signature[filtered_signature['gene'].isin(adata_p2.var_names.tolist())]
filtered_signature = filtered_signature.dropna(how='any')

marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Cell_Subset']].append(row['gene'])

marker_fibro=dict(marker_gene_dictionary)
uc.compute_ucell_scores(adata_Fibro, signatures=marker_fibro, chunk_size=500)
list_sign_fibro = ["{}{}".format(i,"_UCell") for i in list(marker_fibro.keys())]
sc.pl.embedding(adata_Fibro, color=list_sign_fibro,basis='UMAP_IC2_PC10')

#Myeloid All
adata_Myelo = adata_p2[adata_p2.obs['leiden_ICPC'].isin(['1','19','36'])].copy()
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

adata_Myelo.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_myeloid.h5ad")
func.to_seqgeq_myelo(adata_Myelo,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_myeloid.txt",adata_Myelo.var_names)

filtered_signature = signature_v2[signature_v2['Family']=='Myeloid']
filtered_signature = filtered_signature[filtered_signature['gene'].isin(adata_p2.var_names.tolist())]
filtered_signature = filtered_signature.dropna(how='any')

marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Cell_Subset']].append(row['gene'])

marker_myelo=dict(marker_gene_dictionary)
uc.compute_ucell_scores(adata_Myelo, signatures=marker_myelo, chunk_size=500)
list_sign_myelo = ["{}{}".format(i,"_UCell") for i in list(marker_myelo.keys())]
#sc.pl.embedding(adata_Myelo, color=list_sign_myelo,basis='spatial')
sc.pl.embedding(adata_Myelo, color=list_sign_myelo,basis='UMAP_IC2_PC10')

#Myeloid selected
mnp10= list(pd.read_csv("Patient2/Annotation_Leiden_Clusters_Mono.Mac_Merzcope_Lung_Patient2/Extraire et rerun UMAP leiden/A_REANALYSERindépendemment_Subset_leiden_MNP_#10_MoMac.mix.IL4I1.TREM2.csv", skiprows=5).columns)
mnp7= list(pd.read_csv("Patient2/Annotation_Leiden_Clusters_Mono.Mac_Merzcope_Lung_Patient2/Extraire et rerun UMAP leiden/export_Subset_leiden_MNP_#7_MoMac.C1QC.LYZ.csv", skiprows=5).columns)
mnp1= list(pd.read_csv("Patient2/Annotation_Leiden_Clusters_Mono.Mac_Merzcope_Lung_Patient2/Extraire et rerun UMAP leiden/Subset_leiden_MNP_#1_DC2.3.csv", skiprows=5).columns)
mnp2= list(pd.read_csv("Patient2/Annotation_Leiden_Clusters_Mono.Mac_Merzcope_Lung_Patient2/Extraire et rerun UMAP leiden/Subset_leiden_MNP_#2_Mac.FOLR2.csv", skiprows=5).columns)
mnp8= list(pd.read_csv("Patient2/Annotation_Leiden_Clusters_Mono.Mac_Merzcope_Lung_Patient2/Extraire et rerun UMAP leiden/Subset_leiden_MNP_#8_Mac.SPP1.csv", skiprows=5).columns)
mnp12= list(pd.read_csv("Patient2/Annotation_Leiden_Clusters_Mono.Mac_Merzcope_Lung_Patient2/Extraire et rerun UMAP leiden/Subset_leiden_MNP_#12_Mac.IFNB1.csv", skiprows=5).columns)

adata_mnp = adata_p2[adata_p2.obs['cell_id'].isin(mnp10+mnp7+mnp1+mnp2+mnp8+mnp12)]
adata_mnp.obs['leiden_ICPC_full']=adata_mnp.obs['leiden_ICPC'].copy()
adata_mnp.obsm['umapICPC_full']=adata_mnp.obsm['umapICPC'].copy()
adata_mnp.X= adata_mnp.layers["counts"]
sc.pp.normalize_total(adata_mnp)
sc.pp.log1p(adata_mnp)
adata_mnp.layers["lognorm"] = adata_mnp.X.copy()
sc.pp.highly_variable_genes(adata_mnp, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata_mnp, zero_center=False)
adata_mnp.layers["scaled"] = adata_mnp.X.copy()

sc.pp.pca(adata_mnp, n_comps=50)
sc.pp.neighbors(adata_mnp, metric="cosine")
sc.tl.leiden(adata_mnp, flavor="igraph", n_iterations=-1, resolution=1)
sc.tl.umap(adata_mnp, min_dist=0.1)

func.ica(adata_mnp,n_components=50)
ics = pd.DataFrame(adata_mnp.varm['ICs'])
ics.index = adata_mnp.var.index
for top in range(2,6):
    ic_tokeeps = func.ic_tokeep(adata_mnp,ics, signature_v2[signature_v2['Family']=='Myeloid'],"Cell_Subset",topX=top)
    for max_pc in range(10,60,10):
        adata_mnp.obsm["X_ic_pca" ] = np.concatenate((adata_mnp.obsm["X_ica" ][:,ic_tokeeps], adata_mnp.obsm["X_pca" ][:,list(range(0,max_pc))]),axis=1)
        adata_mnp.varm["IC_PCs"] = np.concatenate((adata_mnp.varm["ICs"][:,ic_tokeeps], adata_mnp.varm["PCs" ][:,list(range(0,max_pc))]),axis=1)
        sc.pp.neighbors(adata_mnp, metric="cosine",use_rep="X_ic_pca",key_added="neigh_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.leiden(adata_mnp, flavor="igraph", n_iterations=-1, resolution=1.5,neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc),key_added="leiden_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.umap(adata_mnp,key_added = "UMAP_IC"+str(top)+"_PC"+str(max_pc),neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc), min_dist=0.1)

adata_mnp.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261022_Patient2_myeloid_extractedCells.h5ad")
func.to_seqgeq_myelo(adata_mnp,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261022_Patient2_myeloid_extractedCells.txt",adata_mnp.var_names)

#T CD8
cd8_3= list(pd.read_csv("Patient2/CD8/U_CD8A+_leiden_clusters_PC_3.csv", skiprows=5).columns)
cd8_10= list(pd.read_csv("Patient2/CD8/U_CD8A+_leiden_clusters_PC_10.csv", skiprows=5).columns)
cd8_11= list(pd.read_csv("Patient2/CD8/U_CD8A+_leiden_clusters_PC_11.csv", skiprows=5).columns)
cd8_13= list(pd.read_csv("Patient2/CD8/U_CD8A+_leiden_clusters_PC_13.csv", skiprows=5).columns)
cd8_14= list(pd.read_csv("Patient2/CD8/U_CD8A+_leiden_clusters_PC_14.csv", skiprows=5).columns)
cd8_15= list(pd.read_csv("Patient2/CD8/U_CD8A+_leiden_clusters_PC_15.csv", skiprows=5).columns)

adata_CD8 = adata_T[adata_T.obs['cell_id'].isin(cd8_3+cd8_10+cd8_11+cd8_13+cd8_14+cd8_15)]
adata_CD8.obs['leiden_T']=adata_CD8.obs['leiden'].copy()
adata_CD8.obsm['umapPC_T']=adata_CD8.obsm['X_umap'].copy()
adata_CD8.X= adata_CD8.layers["counts"]
sc.pp.normalize_total(adata_CD8)
sc.pp.log1p(adata_CD8)
adata_CD8.layers["lognorm"] = adata_CD8.X.copy()
sc.pp.highly_variable_genes(adata_CD8, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata_CD8, zero_center=False)
adata_CD8.layers["scaled"] = adata_CD8.X.copy()

sc.pp.pca(adata_CD8, n_comps=50)
sc.pp.neighbors(adata_CD8, metric="cosine")
sc.tl.leiden(adata_CD8, flavor="igraph", n_iterations=-1, resolution=1)
sc.tl.umap(adata_CD8, min_dist=0.1)

func.ica(adata_CD8,n_components=50)
ics = pd.DataFrame(adata_CD8.varm['ICs'])
ics.index = adata_CD8.var.index
for top in range(2,6):
    ic_tokeeps = func.ic_tokeep(adata_CD8,ics, signature_v2[signature_v2['Family']=='Myeloid'],"Cell_Subset",topX=top)
    for max_pc in range(10,60,10):
        adata_CD8.obsm["X_ic_pca" ] = np.concatenate((adata_CD8.obsm["X_ica" ][:,ic_tokeeps], adata_CD8.obsm["X_pca" ][:,list(range(0,max_pc))]),axis=1)
        adata_CD8.varm["IC_PCs"] = np.concatenate((adata_CD8.varm["ICs"][:,ic_tokeeps], adata_CD8.varm["PCs" ][:,list(range(0,max_pc))]),axis=1)
        sc.pp.neighbors(adata_CD8, metric="cosine",use_rep="X_ic_pca",key_added="neigh_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.leiden(adata_CD8, flavor="igraph", n_iterations=-1, resolution=1.5,neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc),key_added="leiden_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.umap(adata_CD8,key_added = "UMAP_IC"+str(top)+"_PC"+str(max_pc),neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc), min_dist=0.1)

adata_CD8.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261026_Patient2_CD8.h5ad")
func.to_seqgeq_myelo(adata_CD8,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261026_Patient2_CD8.txt",adata_CD8.var_names)

with open('/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261026_Patient2_CD8.txt', "a") as f:
    f.write("leiden_T" + "\t" + "\t".join(adata_CD8.obs['leiden_T']) + "\n")
    f.write("umapPC_T_1" + "\t" + "\t".join(adata_CD8.obsm['umapPC_T'][:,0].astype(str)) + "\n")
    f.write("umapPC_T_2" + "\t" + "\t".join(adata_CD8.obsm['umapPC_T'][:,1].astype(str)) + "\n")

with open('/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261026_Patient2_CD8.txt', "a") as f:
    f.write("Voisin_Plasma" + "\t" + "\t".join(adata_CD8.obs['Voisin_Plasma']) + "\n")
    f.write("Voisin_B" + "\t" + "\t".join(adata_CD8.obs['Voisin_B']) + "\n")
    f.write("Voisin_B_and_Plasma" + "\t" + "\t".join(adata_CD8.obs['Voisin_B_and_Plasma']) + "\n")
    f.write("Voisin_0.None_1.B_2.Plasma_3.both" + "\t" + "\t".join(adata_CD8.obs['Voisin_nb']) + "\n")

#readall + Markers
adata_Myelo= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_myeloid.h5ad")
adata_Fibro= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_fibroblast.h5ad")
adata_T= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_T.h5ad")
adata_endoth= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_endothelial.h5ad")
adata_B= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_B.h5ad")
adata_tumor= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_tumour.h5ad")
adata_unsure= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_unsure.h5ad")
adata_plasma= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_plasma.h5ad")

adata=[adata_Myelo,adata_Fibro, adata_T, adata_endoth,adata_B, adata_tumor,adata_unsure,adata_plasma]
name=['Myelo','Fibro', 'T', 'endoth','B','Tumour','Unsure','Plasma']
for i in range(0,len(adata)):
    sc.tl.rank_genes_groups(adata[i], 'leiden', method='wilcoxon', key_added = "wilcoxon")
    markers = sc.get.rank_genes_groups_df(adata[i],group=None,key='wilcoxon')
    #sc.pl.rank_genes_groups_heatmap(adata[i], n_genes=5, key="wilcoxon", groupby="leiden", show_gene_labels=True)
    #sc.pl.rank_genes_groups_dotplot(adata[i],n_genes=5, key="wilcoxon", groupby="leiden")
    markers.to_csv('/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261007_marker_'+name[i]+'_PC.csv')

adata=[adata_Myelo,adata_Fibro, adata_T, adata_endoth]
name=['Myelo','Fibro', 'T', 'endoth']
for i in range(0,len(adata)):
    sc.tl.rank_genes_groups(adata[i], 'leiden_IC2_PC10', method='wilcoxon', key_added = "wilcoxon")
    markers = sc.get.rank_genes_groups_df(adata[i],group=None,key='wilcoxon')
    sc.pl.rank_genes_groups_heatmap(adata[i], n_genes=5, key="wilcoxon", groupby="leiden_IC2_PC10", show_gene_labels=True)
    sc.pl.rank_genes_groups_dotplot(adata[i],n_genes=5, key="wilcoxon", groupby="leiden_IC2_PC10")
    markers.to_csv('/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261007_marker_'+name[i]+'_IC2_PC10.csv')

#Cell pop
cd8_3= list(pd.read_csv("Patient2/CD8/U_CD8A+_leiden_clusters_PC_3.csv", skiprows=5).columns)
cd8_10= list(pd.read_csv("Patient2/CD8/U_CD8A+_leiden_clusters_PC_10.csv", skiprows=5).columns)
cd8_11= list(pd.read_csv("Patient2/CD8/U_CD8A+_leiden_clusters_PC_11.csv", skiprows=5).columns)
cd8_13= list(pd.read_csv("Patient2/CD8/U_CD8A+_leiden_clusters_PC_13.csv", skiprows=5).columns)
cd8_14= list(pd.read_csv("Patient2/CD8/U_CD8A+_leiden_clusters_PC_14.csv", skiprows=5).columns)
cd8_15= list(pd.read_csv("Patient2/CD8/U_CD8A+_leiden_clusters_PC_15.csv", skiprows=5).columns)

adata_p2.obs['Cell_Population'] = 'unknown'
adata_p2.obs.loc[adata_p2.obs['cell_id'].isin(cd8_3+cd8_10+cd8_11+cd8_13+cd8_14+cd8_15), 'Cell_Population'] = 'T.CD8'



#remap to full object
adata_Myelo= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_myeloid.h5ad")
adata_Fibro= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_fibroblast.h5ad")
adata_T= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_T.h5ad")
adata_endoth= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_endothelial.h5ad")
adata_B= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_B.h5ad")
adata_tumor= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_tumour.h5ad")
adata_unsure= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_unsure.h5ad")
adata_plasma= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261013_Patient2_plasma.h5ad")


adata_list=[adata_Myelo,adata_Fibro, adata_T, adata_endoth,adata_B,adata_tumor,adata_unsure,adata_plasma]
for adata in adata_list:
    #adata_p2.obs.loc[adata_p2.obs['cell_id'].isin(adata.obs['cell_id']), 'Cell_Population'] = adata.obs['Cell_population'].astype(str)
    adata_p2.obs.loc[adata_p2.obs['cell_id'].isin(adata.obs['cell_id']), 'Subset_leiden'] = adata.obs['leiden'].astype(str)
    adata_p2.obs.loc[adata_p2.obs['cell_id'].isin(adata.obs['cell_id']), 'Subset_umap_1'] = adata.obsm['X_umap'][:,0]
    adata_p2.obs.loc[adata_p2.obs['cell_id'].isin(adata.obs['cell_id']), 'Subset_umap_2'] = adata.obsm['X_umap'][:,1]
    if 'leiden_IC2_PC10' in adata.obs.columns :
        adata_p2.obs.loc[adata_p2.obs['cell_id'].isin(adata.obs['cell_id']), 'Subset_leiden_IC_PC'] = adata.obs['leiden_IC2_PC10'].astype(str)
        adata_p2.obs.loc[adata_p2.obs['cell_id'].isin(adata.obs['cell_id']), 'Subset_umap_ICPC_1'] = adata.obsm['UMAP_IC2_PC10'][:,0]
        adata_p2.obs.loc[adata_p2.obs['cell_id'].isin(adata.obs['cell_id']), 'Subset_umap_ICPC_2'] = adata.obsm['UMAP_IC2_PC10'][:,1]
