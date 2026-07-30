#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 11:23:48 2026

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
#sopa.settings.auto_save_on_disk = False 

"""Spacial
"""
#sdata_ar = sopa.io.xenium("/home/marine-louarn/Documents/Xenium_Calderaro/Xenium_exp1/Sample6",cells_boundaries=True,nucleus_labels= True, cells_labels=True, nucleus_boundaries= True)
#sdata_ar.write("/home/marine-louarn/Documents/Xenium_Calderaro/Sample6.zarr")
sdata = sd.read_zarr("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1.zarr")

'''
cluster = LocalCluster(n_workers=8,
                       threads_per_worker=4,
                       memory_target_fraction=0.95,
                       memory_limit='16GB')
'''

sopa.segmentation.tissue(sdata)
sopa.make_image_patches(sdata, patch_width=1200, patch_overlap=50)
sopa.utils.get_channel_names(sdata)
sopa.segmentation.cellpose(sdata, channels=["DAPI"], diameter=35, gpu=True)
sopa.make_transcript_patches(sdata, patch_width=500, prior_shapes_key="cellpose_boundaries")

sopa.settings.parallelization_backend = "dask"
sopa.settings.dask_client_kwargs["n_workers"] = 2
sopa.settings.dask_client_kwargs["memory_limit"] = '16GiB'
sopa.settings.dask_client_kwargs["threads_per_worker"] = 4
sopa.segmentation.baysor(sdata, min_area=10)
sopa.aggregate(sdata)

#xenimum explorer
sopa.io.explorer.write('/home/marine-louarn/Documents/Xenium_Calderaro/Sample1/Xenium_explorer_s1', sdata)
sopa.io.explorer.write_cell_categories('/home/marine-louarn/Documents/Xenium_Calderaro/Xenium_explorer_s1', adata_s1)

sdata.pl.render_images("morphology_focus").pl.render_shapes(
    "cellpose_boundaries", outline_alpha=1, fill_alpha=0, outline_color="#ffffff"
).pl.show("global")
sdata.tables["table"].obs["region"] = "cell_boundaries"
sdata.set_table_annotates_spatialelement("table", region="cell_boundaries")

total_count = sdata.tables['table'].obs.total_counts
total_count.index= sdata.tables['table'].obs.cell_id
sdata['cell_boundaries']['total_counts'] = sdata.tables['table'].obs.total_counts
sdata.pl.render_shapes("cell_boundaries", color="total_counts").pl.show()
plt.show()

sdata.pl.render_points(
    "transcripts",
    color="feature_name",
    groups=["ERBB2", "A2ML1"],
    palette=["orange", "blue"],
).pl.show()

"""Transcript
"""
adata = sdata["table"]
adata.layers["counts"] = adata.X.copy()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True)
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts"],
    jitter=0.4,
    multi_panel=True,
)
#sc.pp.scrublet(adata, batch_key="sample")

# -----------------------------------------
# Spatial density (post-QC, log)
# -----------------------------------------
'''xy = adata.obsm["spatial"]
x, y = xy[:, 0], xy[:, 1]

fig, ax = plt.subplots(figsize=(6.5, 6.5))

hb = ax.hexbin(x, y, gridsize=260, bins="log", mincnt=1, cmap="magma")

ax.invert_yaxis()
ax.set_aspect("equal", "box")
ax.set_xlabel("x (µm)")  # colors handled by rcParams
ax.set_ylabel("y (µm)")
ax.set_title("Spatial Density — Post QC", fontsize=14, weight="bold", pad=12, color="black")

# Colorbar styling (not covered by rcParams)
cb = fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.04)
cb.set_label("cells / bin", color="black")
cb.ax.tick_params(colors="black")
cb.outline.set_edgecolor("black")'''


#norm
thres = np.quantile(adata.obs["total_counts"], 0.99)
sc.pp.filter_cells(adata, min_counts=20)
sc.pp.filter_cells(adata, max_counts=thres)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata, zero_center=False)
adata.layers["scaled"] = adata.X.copy()

sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, metric="cosine")
sc.tl.leiden(adata, flavor="igraph", n_iterations=-1, resolution=1,random_state=0)
sc.tl.umap(adata, min_dist=0.1)
sc.pl.umap(adata, color="leiden")


#####
# Back to spacial
#####
sdata.tables["table"].obs["region"] = "cell_boundaries"
sdata.set_table_annotates_spatialelement("table", region="cell_boundaries")


leinden = sdata.tables['table'].obs.leiden
leinden.index= sdata.tables['table'].obs.cell_id
sdata['cell_boundaries']['leiden'] = sdata.tables['table'].obs.leiden
sdata.pl.render_shapes("cell_boundaries", color="leiden").pl.show()
plt.show()

#add ensembl_id                                       
server = biomart.BiomartServer('http://www.ensembl.org/biomart')         
mart = server.datasets['hsapiens_gene_ensembl']
response = mart.search({'attributes': ['hgnc_symbol', 'ensembl_gene_id',]})
data = response.raw.data.decode('ascii')
ensembl_genesymbol = {}
for line in data.splitlines():                                              
    line = line.split('\t')                                                                                                
    gene_symbol = line[0]                                                   
    ensembl_gene = line[1]
    ensembl_genesymbol[gene_symbol] = ensembl_gene
adata.var['gene_id']=adata.var.index
adata.var['gene_id']=adata.var['gene_id'].replace(ensembl_genesymbol)

#Markers
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', key_added = "wilcoxon")
markers = sc.get.rank_genes_groups_df(adata,group=None,key='wilcoxon')
markers.to_csv('fullmarkers_leiden_res1.csv')
sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key="wilcoxon", groupby="leiden", show_gene_labels=True)
sc.pl.rank_genes_groups_dotplot(adata,n_genes=5, key="wilcoxon", groupby="leiden")

#need to filter bf doing dict
signature_marker = pd.read_csv("/home/marine-louarn/ref/20260317_Human_Signature_Genes_for_scRNAseq_data_Extraction.csv")
filtered_signature = signature_marker[signature_marker['Vizgen_Gene'].isin(adata.var_names.tolist())]

marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Population']].append(row['Vizgen_Gene'])
marker_gene_dictionary=dict(marker_gene_dictionary)
sc.pl.dotplot(adata, marker_gene_dictionary, 'leiden_ICPC', dendrogram=True)

adata.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1.h5ad")
adata= sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1_scconcept.h5ad")

#add hand-annotated
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


#novae
novae.spatial_neighbors(adata, radius=80)
novae.plot.connectivities(adata)
model = novae.Novae.from_pretrained("MICS-Lab/novae-brain-0")
model.compute_representations(adata, zero_shot=True)
model.fine_tune(adata)
model.compute_representations(adata)
model.assign_domains(adata, resolution=1)
novae.plot.domains(adata)
novae.plot.spatially_variable_genes(adata, top_k=3, vmax="p95", cell_size=20)
novae.plot.pathway_scores(adata, pathways="/home/marine-louarn/ref/h.all.v2026.1.Hs.json", figsize=(10, 7))
adata.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260528_Sample1_novae.h5ad")

sdata.tables["table"]=adata
sdata.tables["table"].obs["region"] = "cell_boundaries"
sdata.set_table_annotates_spatialelement("table", region="cell_boundaries")

leinden = sdata.tables['table'].obs.novae_domains_7
leinden.index= sdata.tables['table'].obs.cell_id
sdata['cell_boundaries']['novae_domains_7'] = sdata.tables['table'].obs.novae_domains_7
sdata.pl.render_shapes("cell_boundaries", color="novae_domains_7").pl.show()
plt.show()

#after cellTypist
adata_ct= sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/old/Sample1_celltypist.h5ad")
adata_ct.obs["region"] = "cell_boundaries"
sdata.set_table_annotates_spatialelement("table", region="cell_boundaries")

leinden = adata_ct.obs.predicted_labels
leinden.index= adata_ct.obs.cell_id
sdata['cell_boundaries']['predicted_labels'] = adata_ct.obs.predicted_labels
sdata.pl.render_shapes("cell_boundaries", color="predicted_labels").pl.show()
plt.show()

#Markers
sc.tl.rank_genes_groups(adata, 'predicted_labels', method='wilcoxon', key_added = "wilcoxon_celltypist")
markers = sc.get.rank_genes_groups_df(adata,group=None,key='wilcoxon_celltypist')
sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key="wilcoxon_celltypist", groupby="predicted_labels", show_gene_labels=True)
sc.pl.rank_genes_groups_dotplot(adata,n_genes=5, key="wilcoxon_celltypist", groupby="predicted_labels")

