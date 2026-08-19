#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Random Tests
Created on 2026-08-12

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
import novae
import biomart
from collections import defaultdict
import pyucell as uc
import statistics 
import functions_spatial as func
import plotly.graph_objects as go

#test for different conda for novae computation (res1 : novae, spacial: spacialPy)
test_novae = sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260811_Sample1_novaeV2.h5ad")
df= test_novae.obs[['novae_domains_res1','novae_domains_res1_spacial']]
df["val"] = 1
newdf = pd.pivot_table(df, values="val", index="novae_domains_res1", columns="novae_domains_res1_spacial", aggfunc="sum")

sns.heatmap(newdf,annot=True)
plt.show()

df2=newdf.div(newdf.sum(axis=1), axis=0) #normalisation by row
sns.heatmap(df2,annot=True)
plt.show()

df3=newdf.T.div(newdf.T.sum(axis=1), axis=0) #normalisation by row
sns.heatmap(df3.T,annot=True)
plt.show()

#sankey chart
df2 = df.groupby(['novae_domains_res1','novae_domains_res1_spacial']).size().reset_index()
df1 = df.groupby(['novae_domains_res1','novae_domains_res1_spacial']).size().reset_index()
df1.columns = ['novae_domains_res1', 'novae_domains_res1_spacial', 'value']
#df1['novae_domains_res1']=int(''.join(i for i in df1['novae_domains_res1'].astype(str) if i.isdigit()))
df1['novae_domains_res1'] = df1['novae_domains_res1'].astype(str).map(lambda x: x.lstrip('L')).astype(int)-1
df1['novae_domains_res1_spacial'] = df1['novae_domains_res1_spacial'].astype(str).map(lambda x: x.lstrip('L')).astype(int)+1+max(df1['novae_domains_res1'])
df1['value']=pd.to_numeric(df1['value'])
links_dict = df1.to_dict(orient='list')
unique_source_target = list(pd.unique(df1[['novae_domains_res1', 'novae_domains_res1_spacial']].values.ravel('K')))

fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = ['L1','L2','L3','L4','L5','L6','L7','L8','L9','L10','L11','L0_v2','L1_v2','L2_v2','L3_v2','L4_v2','L5_v2','L6_v2','L7_v2','L9_v2'],
      color = "blue"
      
    ),
    link = dict(
      source = links_dict["novae_domains_res1"],
      target = links_dict["novae_domains_res1_spacial"],
      value = links_dict["value"],
  
  ))])
fig.update_layout(title_text="Novae using SOPA env. vs novae using NOVAE env.", font_size=10)
fig.show() 


#nuc vs allcell
s2_cell = sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample2/20260812_Sample2_xeniumSegCell_seqgeq.h5ad")
s2_nuc = sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample2/20260812_Sample2_xeniumSegNuc_seqgeq.h5ad")

sc.pl.embedding(s2_nuc, color=["n_counts","area"],basis='umapICPC')
sc.pl.embedding(s2_cell, color=["n_counts","cell_area","nucleus_area"],basis='umapICPC')


#sample4 segmentation tests
sdata_ar = sopa.io.xenium("/home/marine-louarn/Documents/Xenium_Calderaro/Xenium_exp1/Sample4",cells_boundaries=True,nucleus_labels= True, cells_labels=True, nucleus_boundaries= True)
sdata_ar.write("/home/marine-louarn/Documents/Xenium_Calderaro/Sample4/Sample4_Xenium_nuc.zarr")
sdata_xeniumNuc = sd.read_zarr("/home/marine-louarn/Documents/Xenium_Calderaro/Sample4/Sample4_Xenium_nuc.zarr")
sopa.aggregate(sdata_xeniumNuc,shapes_key='nucleus_boundaries')


sdata_baysor20 = sd.read_zarr("/home/marine-louarn/Documents/Xenium_Calderaro/Sample4/Sample4_baysor_min20.zarr")
sdata_cellpose = sd.read_zarr("/home/marine-louarn/Documents/Xenium_Calderaro/Sample4/Sample4_onlyCellpose.zarr")
sdata_xenium_cell = sd.read_zarr("/home/marine-louarn/Documents/Xenium_Calderaro/Sample4/Sample4_Xenium.zarr")
sdata_xenium_nuc = sd.read_zarr("/home/marine-louarn/Documents/Xenium_Calderaro/Sample4/Sample4_Xenium_nuc.zarr")

baysor10 = sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample4/20260720_Sample4.h5ad")
baysor20 = sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample4/20260817_Sample4_baysor20.h5ad")
cellpose = sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample4/20260817_Sample4_cellpose.h5ad")
xenium_cell = sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample4/20260817_Sample4_xenium_cell.h5ad")
xenium_nuc = sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample4/20260817_Sample4_xenium_nuc.h5ad")

list_adata =[baysor10,baysor20,cellpose,xenium_cell,xenium_nuc]
list_name=['baysor10','baysor20','cellpose','xenium_cell','xenium_nuc']
for i in range(0,5):
    print(list_name[i])
    print('nb cell:',len(list_adata[i].obs['leiden_ICPC']))
    print('nb cluster:',max(list_adata[i].obs['leiden_ICPC'].astype(int)))
    print('mean n_count',np.mean(list_adata[i].obs['n_counts']))
    print('min n_count',min(list_adata[i].obs['n_counts']))
    print('max n_count',max(list_adata[i].obs['n_counts']))
    if 'area' in list_adata[i].obs.columns :
        print('median area:',statistics.median(list_adata[i].obs['area']))
        print('mean area:',np.mean(list_adata[i].obs['area']))
        print('min area:',min(list_adata[i].obs['area']))
        print('max area:',max(list_adata[i].obs['area']))
    else :
        print('median area:',statistics.median(list_adata[i].obs['cell_area']))
        print('mean area:',np.mean(list_adata[i].obs['cell_area']))
        print('min area:',min(list_adata[i].obs['cell_area']))
        print('max area:',max(list_adata[i].obs['cell_area']))
    #sc.pl.embedding(list_adata[i], color="leiden_ICPC",basis='umapICPC',legend_loc='on data')

signature = pd.read_csv("/home/marine-louarn/Documents/test/20260204_Signatures_Cell_populations_HUMAN.csv")
filtered_signature = signature[signature['gene'].isin(baysor20.var_names.tolist())]
marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Family']].append(row['gene'])

marker_gene_dictionary=dict(marker_gene_dictionary)

list_adata =[baysor10,baysor20,cellpose,xenium_cell,xenium_nuc]
list_name=['baysor10','baysor20','cellpose','xenium_cell','xenium_nuc']
for i in range(0,5):
    print(list_name[i])
    uc.compute_ucell_scores(list_adata[i], signatures=marker_gene_dictionary, chunk_size=500)
    sc.pl.embedding(list_adata[i], color="Myeloid_UCell",basis='umapICPC')

#sample 2
adata_s2=sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample2_issues.h5ad")
sc.pl.embedding(adata_s2, color="leiden_ICPC",basis='umapICPC',legend_loc='on data')
#(adata_s2.X!=adata_s2.layers['lognorm']).nnz==0 #check if X is layers lognorm
