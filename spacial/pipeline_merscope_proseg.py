#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:23:48 2026

@author: marine-louarn
"""

import spatialdata as sd
import sopa 
from spatialdata_io import xenium
from pathlib import Path
import shutil
import spatialdata_plot
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np
import pandas as pd
#import tangram as tg
import os
import novae
from collections import defaultdict
import importlib
import squidpy as sq
import pyucell as uc
import dask.dataframe as dd
import anndata

import functions_spatial as func #importlib.reload(func)

signature = pd.read_csv("/home/marine-louarn/ref/20260317_Human_Signature_Genes_for_scRNAseq_data_Extraction.csv")

"""Spacial
"""
#sdata = sopa.io.merscope("/home/marine-louarn/Documents/Lung_Spatial/HumanLungCancerPatient2")
#sdata.write("/home/marine-louarn/Documents/Lung_Spatial/Patient2_proseg.zarr")
sdata_p1 = sd.read_zarr("/home/marine-louarn/Documents/Lung_Spatial/Patient1_proseg.zarr") 
sdata_p2 = sd.read_zarr("/home/marine-louarn/Documents/Lung_Spatial/Patient2_proseg.zarr") 



sopa.segmentation.tissue(sdata_p2)
sopa.make_image_patches(sdata_p2, patch_width=1200, patch_overlap=50)
sopa.segmentation.cellpose(sdata_p2, channels=["DAPI"], diameter=35, gpu=True)
 
sopa.settings.parallelization_backend = "dask"
sopa.settings.dask_client_kwargs["n_workers"] = 2
sopa.settings.dask_client_kwargs["memory_limit"] = '16GiB'
sopa.settings.dask_client_kwargs["threads_per_worker"] = 4

sopa.make_transcript_patches(sdata_p2, patch_width=None, prior_shapes_key="cellpose_boundaries") #cellpose_boundaries
sopa.segmentation.proseg(sdata_p2, command_line_suffix= '--cell-compactness 0.01')
sopa.aggregate(sdata_p2)

#xenimum explorer
points_key = list(sdata_p2.points.keys())[0]
points = sdata_p2.points[points_key]
attrs = dict(points.attrs)
points_fixed = points.compute()
points_fixed = dd.from_pandas(points_fixed, npartitions=1)
points_fixed.attrs.update(attrs)
sdata_p2.points[points_key] = points_fixed

sopa.io.explorer.write('/home/marine-louarn/Documents/Lung_Spatial/Patient2_proseg', sdata_p2)
#sopa.io.explorer.write_cell_categories('/home/marine-louarn/Documents/Xenium_Calderaro/Xenium_explorer_s1', adata_s1)
#sd.sanitize_table(adata)

"""Transcript
"""
adata = sdata_p2["table"]
adata.layers["counts"] = adata.X.copy()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True)
'''sc.pl.violin(
    adata,
    ["n_genes_by_counts","total_counts", "pct_counts_mt"],
    jitter=0.4, 
    multi_panel=True,
)'''

thres = np.quantile(adata.obs["total_counts"], 0.99)
sc.pp.filter_cells(adata, min_counts=20)
sc.pp.filter_cells(adata, max_counts=thres)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=100)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()
sc.pp.scale(adata, zero_center=False)
adata.layers["scaled"] = adata.X.copy()

#PC
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, metric="cosine")
sc.tl.leiden(adata, flavor="igraph", n_iterations=-1, resolution=1,random_state=0)
sc.tl.umap(adata, min_dist=0.1)

#ICPC
signature = pd.read_csv("/home/marine-louarn/Documents/test/20260204_Signatures_Cell_populations_HUMAN.csv")

#sc.pp.pca(adata, n_comps=50)
func.ica(adata,n_components=50,max_iter=1000)
ics = pd.DataFrame(adata.varm['ICs'])
ics.index = adata.var.index


ic_tokeeps = func.ic_tokeep(adata,ics, signature,"Major_Cell_Populations")
adata.obsm["X_ic_pca" ] = np.concatenate((adata.obsm["X_ica" ][:,ic_tokeeps], adata.obsm["X_pca" ]),axis=1)
adata.varm["IC_PCs"] = np.concatenate((adata.varm["ICs"][:,ic_tokeeps], adata.varm["PCs" ]),axis=1)

sc.pp.neighbors(adata, metric="cosine",use_rep="X_ic_pca",key_added='neigh_ICPC')
sc.tl.leiden(adata, flavor="igraph", n_iterations=-1, resolution=1,neighbors_key='neigh_ICPC',key_added='leiden_ICPC')
sc.tl.umap(adata,key_added = 'umapICPC',neighbors_key='neigh_ICPC', min_dist=0.1)
#sc.pl.embedding(adata, color="leiden_ICPC",basis='umapICPC')

adata.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20260922_Patient2_cellposeProseg.h5ad")

header_sg="/home/marine-louarn/Documents/Header_SeqGeq.txt"
signature = pd.read_csv("/home/marine-louarn/Documents/test/20260204_Signatures_Cell_populations_HUMAN.csv")
genelist= signature['gene'].tolist()
func.to_seqgeq(adata,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/20260922_Patient2_cellposeProseg_full.txt",adata.var_names)

#analysis
adata_p1= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20260918_Patient1_cellposeProseg.h5ad")
adata_p2= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20260922_Patient2_cellposeProseg.h5ad")
#signature
signature = pd.read_csv("/home/marine-louarn/Documents/test/20260204_Signatures_Cell_populations_HUMAN.csv")

filtered_signature = signature[signature['gene'].isin(adata_p1.var_names.tolist())]
marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Family']].append(row['gene'])

marker_gene_dictionary=dict(marker_gene_dictionary)
uc.compute_ucell_scores(adata_p1, signatures=marker_gene_dictionary, chunk_size=500)
uc.compute_ucell_scores(adata_p2, signatures=marker_gene_dictionary, chunk_size=500)
sc.pl.embedding(adata_p1, color=["Myeloid_UCell","leiden_ICPC"], cmap="viridis", size=20,basis='umapICPC')
sc.pl.embedding(adata_p1, color="leiden_ICPC",basis='spatial')

func.to_seqgeq(adata_p1,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/20260922_Patient1_cellposeProseg_full.txt",adata_p1.var_names.tolist())
func.to_seqgeq(adata_p2,header_sg,"/home/marine-louarn/Documents/Lung_Spatial/20260922_Patient2_cellposeProseg_full.txt",adata_p2.var_names.tolist())

#novae
novae.spatial_neighbors(adata_p1, radius=80)
novae.plot.connectivities(adata_p1)
model = novae.Novae.from_pretrained("MICS-Lab/novae-human-0")
model.compute_representations(adata_p1, zero_shot=True)
model.fine_tune(adata_p1,accelerator="gpu", num_workers=4)
model.compute_representations(adata_p1,accelerator="gpu", num_workers=4)
model.assign_domains(adata_p1, resolution=1)
novae.plot.domains(adata_p1)
sq.pl.spatial_scatter(adata_p1, color="novae_domains_res1", shape=None, palette = 'Set3')
novae.plot.spatially_variable_genes(adata_p1, top_k=3, vmax="p95", cell_size=20)
novae.plot.pathway_scores(adata_p1, pathways="/home/marine-louarn/ref/h.all.v2026.1.Hs.json", figsize=(10, 7))
adata_p1.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20260923_Patient1_proseg_novae.h5ad")

sc.pl.embedding(adata_p1, color="novae_domains_res1",basis='spatial')
#adata_p1=sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20260923_Patient1_proseg_novae.h5ad")
#adata_p2=sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20260923_Patient2_proseg_novae.h5ad")

adata_novae = sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20261007_Patient1_novae.h5ad")
df= adata_novae.obs[['novae_domains_res1','Cell_Type']]
df["val"] = 1
newdf = pd.pivot_table(df, values="val", index="novae_domains_res1", columns="Cell_Type", aggfunc="sum")

sns.heatmap(newdf,annot=True)
plt.show()

df2=newdf.div(newdf.sum(axis=1), axis=0) #normalisation by row
sns.heatmap(df2,annot=True)
plt.show()

df3=newdf.T.div(newdf.T.sum(axis=1), axis=0) #normalisation by row
sns.heatmap(df3.T,annot=True)
plt.show()
adata_p1.obs['novae_domains_res1_nb']='0'
adata_p1.obs.loc[adata_p1.obs['novae_domains_res1'].isin(['L0']), 'novae_domains_res1_nb'] = '1'
adata_p1.obs.loc[adata_p1.obs['novae_domains_res1'].isin(['L1']), 'novae_domains_res1_nb'] = '2'
adata_p1.obs.loc[adata_p1.obs['novae_domains_res1'].isin(['L2']), 'novae_domains_res1_nb'] = '3'
adata_p1.obs.loc[adata_p1.obs['novae_domains_res1'].isin(['L3']), 'novae_domains_res1_nb'] = '4'

adata_p1.obs['novae_domains_7_nb']='0'
adata_p1.obs.loc[adata_p1.obs['novae_domains_7'].isin(['D1001']), 'novae_domains_7_nb'] = '1'
adata_p1.obs.loc[adata_p1.obs['novae_domains_7'].isin(['D1004']), 'novae_domains_7_nb'] = '2'
adata_p1.obs.loc[adata_p1.obs['novae_domains_7'].isin(['D1013']), 'novae_domains_7_nb'] = '3'
adata_p1.obs.loc[adata_p1.obs['novae_domains_7'].isin(['D1016']), 'novae_domains_7_nb'] = '4'
adata_p1.obs.loc[adata_p1.obs['novae_domains_7'].isin(['D980']), 'novae_domains_7_nb'] = '5'
adata_p1.obs.loc[adata_p1.obs['novae_domains_7'].isin(['D993']), 'novae_domains_7_nb'] = '6'
adata_p1.obs.loc[adata_p1.obs['novae_domains_7'].isin(['D994']), 'novae_domains_7_nb'] = '7'

adata_p1.obs['novae_domains_15_nb']='0'
adata_p1.obs.loc[adata_p1.obs['novae_domains_15'].isin(['D1004']), 'novae_domains_15_nb'] = '1'
adata_p1.obs.loc[adata_p1.obs['novae_domains_15'].isin(['D1000']), 'novae_domains_15_nb'] = '2'
adata_p1.obs.loc[adata_p1.obs['novae_domains_15'].isin(['D1001']), 'novae_domains_15_nb'] = '3'
adata_p1.obs.loc[adata_p1.obs['novae_domains_15'].isin(['D993']), 'novae_domains_15_nb'] = '4'
adata_p1.obs.loc[adata_p1.obs['novae_domains_15'].isin(['D1007']), 'novae_domains_15_nb'] = '5'
adata_p1.obs.loc[adata_p1.obs['novae_domains_15'].isin(['D1006']), 'novae_domains_15_nb'] = '6'
adata_p1.obs.loc[adata_p1.obs['novae_domains_15'].isin(['D980']), 'novae_domains_15_nb'] = '7'
adata_p1.obs.loc[adata_p1.obs['novae_domains_15'].isin(['D988']), 'novae_domains_15_nb'] = '8'
adata_p1.obs.loc[adata_p1.obs['novae_domains_15'].isin(['D1008']), 'novae_domains_15_nb'] = '9'
adata_p1.obs.loc[adata_p1.obs['novae_domains_15'].isin(['D994']), 'novae_domains_15_nb'] = '10'
adata_p1.obs.loc[adata_p1.obs['novae_domains_15'].isin(['D998']), 'novae_domains_15_nb'] = '11'
adata_p1.obs.loc[adata_p1.obs['novae_domains_15'].isin(['D990']), 'novae_domains_15_nb'] = '12'
adata_p1.obs.loc[adata_p1.obs['novae_domains_15'].isin(['D983']), 'novae_domains_15_nb'] = '13'
adata_p1.obs.loc[adata_p1.obs['novae_domains_15'].isin(['D976']), 'novae_domains_15_nb'] = '14'

adata_p1.obs['novae_domains_20_nb']='0'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D1000']), 'novae_domains_20_nb'] = '1'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D981']), 'novae_domains_20_nb'] = '2'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D1001']), 'novae_domains_20_nb'] = '3'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D993']), 'novae_domains_20_nb'] = '4'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D896']), 'novae_domains_20_nb'] = '5'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D1002']), 'novae_domains_20_nb'] = '6'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D980']), 'novae_domains_20_nb'] = '7'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D999']), 'novae_domains_20_nb'] = '8'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D992']), 'novae_domains_20_nb'] = '9'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D979']), 'novae_domains_20_nb'] = '10'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D988']), 'novae_domains_20_nb'] = '11'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D996']), 'novae_domains_20_nb'] = '12'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D994']), 'novae_domains_20_nb'] = '13'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D998']), 'novae_domains_20_nb'] = '14'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D990']), 'novae_domains_20_nb'] = '15'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D983']), 'novae_domains_20_nb'] = '16'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D976']), 'novae_domains_20_nb'] = '17'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D971']), 'novae_domains_20_nb'] = '18'
adata_p1.obs.loc[adata_p1.obs['novae_domains_20'].isin(['D910']), 'novae_domains_20_nb'] = '19'

with open("/home/marine-louarn/Documents/Lung_Spatial/20261008_Patient1_cellposeProseg_novae.txt", "w") as f:
    novae_domains_res1_nb = adata_p1.obs['novae_domains_res1_nb']
    novae_domains_res1_nb_str = [str(x) for x in novae_domains_res1_nb]
    f.write("novae_domains_res1" + "\t" + "\t".join(novae_domains_res1_nb_str) + "\n")
    novae_domains_7_nb = adata_p1.obs['novae_domains_7_nb']
    novae_domains_7_nb_str = [str(x) for x in novae_domains_7_nb]
    f.write("novae_domains_7" + "\t" + "\t".join(novae_domains_7_nb_str) + "\n")
    novae_domains_15_nb = adata_p1.obs['novae_domains_15_nb']
    novae_domains_15_nb_str = [str(x) for x in novae_domains_15_nb]
    f.write("novae_domains_15" + "\t" + "\t".join(novae_domains_15_nb_str) + "\n")
    novae_domains_20_nb = adata_p1.obs['novae_domains_20_nb']
    novae_domains_20_nb_str = [str(x) for x in novae_domains_20_nb]
    f.write("novae_domains_20" + "\t" + "\t".join(novae_domains_20_nb_str) + "\n")

import matplotlib.colors as mcolors
cols = mcolors.ListedColormap([ '#1f77b4ff','#ff7f0eff','#2ca02cff','#d62728ff','#9467bdff','#8c564bff','#e377c2ff','#7f7f7fff'])



df = adata_p1.obs[['Cell_Type','novae_domains_15_nb']]
df_new = df.pivot_table(columns="Cell_Type", index="novae_domains_15_nb", aggfunc=len, fill_value=0)
ax = df_new.plot.bar(stacked=True,cmap=cols)
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Population", framealpha=0);
plt.show()


#domains markers
sc.tl.rank_genes_groups(adata_p2, groupby="novae_domains_res1", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata_p2, groupby="novae_domains_res1", standard_scale="var", n_genes=5)

#Alexandre signature
sign_tum = pd.read_csv("/home/marine-louarn/Documents/Lung_Spatial/Alexandre_Paper_FOLR2/signature_tum.csv")
sign_norm = pd.read_csv("/home/marine-louarn/Documents/Lung_Spatial/Alexandre_Paper_FOLR2/signature.csv")

filtered_norm = sign_norm[sign_norm['gene'].isin(adata_p1.var_names.tolist())]
marker_norm= defaultdict(list)
for idx, row in filtered_norm.iterrows():
    name = row['cluster'] + "_normal"
    marker_norm[name].append(row['gene'])

marker_norm=dict(marker_norm)
uc.compute_ucell_scores(adata_p1, signatures=marker_norm, chunk_size=500)

filtered_tum = sign_tum[sign_tum['gene'].isin(adata_p1.var_names.tolist())]
marker_tum = defaultdict(list)
for idx, row in filtered_tum.iterrows():
    name = row['cluster'] + "_tumoral"
    marker_tum[name].append(row['gene'])

marker_tum=dict(marker_tum)
uc.compute_ucell_scores(adata_p1, signatures=marker_tum, chunk_size=500)

list_sign_tum = ["{}{}".format(i,"_UCell") for i in list(marker_tum.keys())]
sc.pl.embedding(adata_p1, color=list_sign_tum,basis='spatial')
list_sign_norm = ["{}{}".format(i,"_UCell") for i in list(marker_norm.keys())]
sc.pl.embedding(adata_p1, color=list_sign_norm,basis='spatial')

#Markers
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', key_added = "wilcoxon")
markers = sc.get.rank_genes_groups_df(adata,group=None,key='wilcoxon')
markers.to_csv('fullmarkers_leiden_res05.csv')
sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key="wilcoxon", groupby="leiden", show_gene_labels=True)
sc.pl.rank_genes_groups_dotplot(adata,n_genes=5, key="wilcoxon", groupby="leiden")

#dendro
adata_p1=sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20260805_Patient1_ICPC.h5ad")
adata_p2=sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20260805_Patient2_ICPC.h5ad")

#need to filter bf doing dict
signature_v2 = pd.read_csv("/home/marine-louarn/ref/20260928_JDD_BreastK_DC1update.csv")
filtered_signature= signature_v2[signature_v2['gene'].isin(adata_p1.var_names.tolist())]
marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Family']].append(row['gene'])

marker_gene_dictionary=dict(marker_gene_dictionary)
sc.pl.dotplot(adata_p1, marker_gene_dictionary, 'leiden', dendrogram=True)
uc.compute_ucell_scores(adata_p1, signatures=marker_gene_dictionary, chunk_size=500)

a=['Endothelial','Epithelial','Fibroblast','Lympho','Myeloid']
a = [x+"_UCell" for x in a]
df = pd.DataFrame(columns=a,index=adata_p1.obs['leiden_ICPC'].unique())
for j in adata_p1.obs['leiden_ICPC'].unique():
    i_index= adata_p1.obs[adata_p1.obs['leiden_ICPC']==j]
    mat=i_index[a].copy()
    mean=mat.mean(axis=0)
    df.loc[j]=pd.to_numeric(mean)

df = df[a].apply(pd.to_numeric)
df=df.T
sns.heatmap(df,cmap="crest")
plt.show()

#neighborhood
sq.gr.spatial_neighbors(adata_p1, coord_type="generic", delaunay=True)
sq.gr.centrality_scores(adata_p1, cluster_key="leiden_ICPC")
sq.pl.centrality_scores(adata_p1, cluster_key="leiden_ICPC", figsize=(16, 5))
plt.show()

sq.gr.co_occurrence(adata_p1,cluster_key="leiden_ICPC",)
sq.pl.co_occurrence(adata_p1,cluster_key="leiden_ICPC", clusters="7",figsize=(10, 10))
plt.show()
sq.pl.spatial_scatter(adata_p1, color="leiden_ICPC", shape=None, size=2,) #identic to sc.pl.embedding(adata_p1, color="leiden_ICPC",basis='spatial')
plt.show()
sq.gr.interaction_matrix(adata_p1, cluster_key="leiden_ICPC")
sq.pl.interaction_matrix(adata_p1, cluster_key="leiden_ICPC")
plt.show()

sq.gr.nhood_enrichment(adata_p1, cluster_key="leiden_ICPC")
sq.pl.nhood_enrichment(adata_p1, cluster_key="leiden_ICPC")
plt.show()
mode = "L"
sq.gr.ripley(adata_p1, cluster_key="leiden_ICPC", mode=mode)
sq.pl.ripley(adata_p1, cluster_key="leiden_ICPC", mode=mode)
plt.show()



#keep only with X neighboor
close_to = 'tumour'
target_indices = np.where(adata_p1.obs['Cell_Type'] == close_to)[0]
adj_matrix = adata_p1.obsp["spatial_connectivities"]
neighbor_counts = np.array(adj_matrix[:, target_indices].sum(axis=1))
np.histogram(neighbor_counts)
plt.hist(neighbor_counts)
plt.title("Tumoural neighbour")
plt.show()

adata_filtered = adata_p1[neighbor_counts > 5].copy()
sq.pl.spatial_scatter(adata_filtered, color="Cell_Type", shape=None, size=2,)
plt.show()

adata_p1.obs['Tumour_proximity']=neighbor_counts
adata_tmp = adata_p1[adata_p1.obs['Cell_Type'] == 'Myeloid'].copy()
adata_tmp.obs['Tumour_k5'] = ['Y' if int(i) >4 else 'N' for i in adata_tmp.obs['Tumour_proximity']]

df=adata_tmp.obs[['Cell_Population','Tumour_k5']]
df_new = df.pivot_table(columns="Cell_Population", index="Tumour_k5", aggfunc=len, fill_value=0)
ax = df_new.plot.bar(stacked=True,cmap='tab20')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Population", framealpha=0);
plt.show()

sq.pl.spatial_scatter(adata_tmp, color=["Cell_Population","Tumour_k5"], shape=None, size=2,palette='tab20')
plt.show()

import matplotlib.colors as mcolors
cols = mcolors.ListedColormap([ '#1f77b4ff','#ff7f0eff'])
sq.pl.spatial_scatter(adata_tmp, library_id="spatial",shape=None,color='Tumour_k5',palette=cols)
cols = mcolors.ListedColormap([ '#1f77b4ff','#ff7f0eff','#2ca02cff','#d62728ff','#9467bdff','#8c564bff','#e377c2ff','#7f7f7fff'])

df_new = df.pivot_table(columns="Cell_Population", index="Tumour_k5", aggfunc=len, fill_value=0).apply(lambda x: x/sum(x), axis=1)
ax = df_new.plot.bar(stacked=True,cmap='tab20')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Population", framealpha=0);
plt.show()

target_indices = np.where(adata_p1.obs['Cell_Type'] == 'B cells and plasma cells')[0]
adj_matrix = adata_p1.obsp["spatial_connectivities"]
neighbor_counts = np.array(adj_matrix[:, target_indices].sum(axis=1))
np.histogram(neighbor_counts)
plt.hist(neighbor_counts)
plt.title("B neighbour")
plt.show()

adata_filtered = adata_p1[neighbor_counts > 5].copy()
sq.pl.spatial_scatter(adata_filtered, color="Cell_Type", shape=None, size=2,)
plt.show()

adata_p1.obs['B_proximity']=neighbor_counts

adata_T = adata_p1[adata_p1.obs['leiden_ICPC'].isin(['8','12'])].copy()
adata_T.obs['B_neigh'] = ['Y' if int(i) >1 else 'N' for i in adata_T.obs['B_proximity']]
marker_B_neigh = sc.tl.rank_genes_groups(adata_T, groupby='B_neigh', key_added='group_DE_results')
sc.pl.rank_genes_groups_heatmap(adata_T, groupby='B_neigh',key='group_DE_results')
sc.pl.rank_genes_groups_dotplot(adata_T, groupby='B_neigh',key='group_DE_results')

'''
Integration
'''
#ingest 
adata_p1= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20260918_Patient1_cellposeProseg.h5ad")
adata_p2= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20260922_Patient2_cellposeProseg.h5ad")

overlap = list(set(adata_p1.var_names.tolist())&set(adata_p2.var_names.tolist()))
adata_ingest_1=adata_p1[:,overlap].copy()
adata_ingest_1.obsm['X_umap']=adata_ingest_1.obsm['umapICPC'].copy()
adata_ingest_1.uns['umap']=adata_ingest_1.uns['umapICPC'].copy()
adata_ingest_2=adata_p2[:,overlap].copy()

ingest = sc.tl.ingest(adata_ingest_2, adata_ingest_1, obs="leiden_ICPC",inplace=False)
ingest.uns["leiden_ICPC_colors"] = adata_ingest_1.uns["leiden_ICPC_colors"]  # fix colors
sc.pl.umap(ingest, color=["leiden_ICPC"], wspace=0.5)
sc.pl.umap(adata_ingest_1, color=["leiden_ICPC"], wspace=0.5)
sc.pl.embedding(adata_ingest_1, color=["leiden_ICPC"],basis='umapICPC')


s1s2_concat = anndata.concat([adata_ingest_1, ingest], label="Patient", keys=["Patient1", "Patient2"])
s1s2_concat.obs_names_make_unique()
s1s2_concat.obs["leiden_ICPC"] = (
    s1s2_concat.obs["leiden_ICPC"].astype("category").cat.reorder_categories(adata_ingest_1.obs["leiden_ICPC"].cat.categories)
)
sc.pl.umap(s1s2_concat, color=["Patient", "leiden_ICPC"])
sc.pl.embedding(s1s2_concat, color=["Patient", "leiden_ICPC"],basis='umapICPC')

s1s2_concat.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20261006_ingest.h5ad")

#bbknn: no good batch correction
adata_concat = anndata.concat([adata_p1, adata_p2], label="Patient", keys=["Patient1", "Patient2"])
adata_concat.obs_names_make_unique()
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key="Patient") 
sc.tl.umap(adata_concat,key_added = 'bbknn')
sc.pl.embedding(adata_concat, color=["Patient"],basis='bbknn')
adata_concat.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20261006_bbknn.h5ad")


#harmony
def integrate_harmony(adata, batch_effect_key):
    sc.pp.pca(adata, n_comps=50)
    harmony_output = hm.run_harmony(adata.obsm['X_pca'], 
                                    adata.obs, batch_effect_key)
    adata.obsm['X_pcaHar'] = harmony_output.Z_corr
    sc.pp.neighbors(adata, n_neighbors=30, 
            n_pcs=10, use_rep='X_pcaHar')
    sc.tl.umap(adata, min_dist=0.05,key_added = 'harmony')
    sc.tl.leiden(adata, flavor="igraph", n_iterations=-1,key_added='leiden_harmony')
    return(adata)

adata_concat = anndata.concat([adata_p1, adata_p2], label="Patient", keys=["Patient1", "Patient2"])
adata_concat.obs_names_make_unique()
harmony_int = integrate_harmony(adata_concat,'Patient')

sc.pl.embedding(harmony_int, color=["leiden_harmony","Patient"],basis="harmony")
harmony_int.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20261006_harmony.h5ad")

ingest= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20261006_ingest.h5ad")
harmony= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20261006_harmony.h5ad")
bbknn= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/20261006_bbknn.h5ad")