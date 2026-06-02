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
import biomart
os.chdir("/home/marine-louarn/Documents/Spacial_Roussy/documents_20260410/SFTP_SHARING_ARC_PGA/")

def save_show_close(fig, dpi: int = 300):
    fig.tight_layout()
    plt.show()
    plt.close(fig)

signature = pd.read_csv("/home/marine-louarn/ref/20260317_Human_Signature_Genes_for_scRNAseq_data_Extraction.csv")

"""Spacial
"""
sdata = sopa.io.merscope("/home/marine-louarn/Documents/Spacial_Roussy/documents_20260410/SFTP_SHARING_ARC_PGA/region_20EN-1063")
sdata.write("region_20EN-1063_v1.zarr")
sdata = sd.read_zarr("region_20EN-1063_v1.zarr") 

adata_gr = sc.read_h5ad("region_20EN-1063/202508221546_TEST-RT03825-5_VMSC19302_region_20EN-1063.h5ad")
sopa.io.explorer.write(
    "region_24.explorer",
    sdata,
    image_key="SFTP_SHARING_ARC_PGA_region_24EN-0017_z3",
    shapes_key="cellpose_boundaries",
    points_key="SFTP_SHARING_ARC_PGA_region_24EN-0017_transcripts",
    gene_column="gene",
    lazy=True,
    ram_threshold_gb=16,
)

sopa.segmentation.tissue(sdata)
sopa.make_image_patches(sdata, patch_width=1200, patch_overlap=50)
sopa.utils.get_channel_names(sdata)
sopa.settings.parallelization_backend = "dask"
sopa.segmentation.cellpose(sdata, channels=["DAPI"], diameter=35)#, gpu=True)
sopa.make_transcript_patches(sdata, patch_width=500, prior_shapes_key="cellpose_boundaries")
sopa.segmentation.baysor(sdata, min_area=10) #does not run outside of terminal
sopa.aggregate(sdata)


#sd.sanitize_table(adata)

sdata.pl.render_images("SFTP_SHARING_ARC_PGA_region_20EN-1063_z3").pl.render_shapes(
    "baysor_boundaries", outline_alpha=1, fill_alpha=0, outline_color="#ffffff"
).pl.show("global")
plt.show()

"""Transcript
"""
adata = sdata["table"]
adata.layers["counts"] = adata.X.copy()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True)
sc.pl.violin(
    adata,
    ["n_genes_by_counts","total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

thres = np.quantile(adata.obs["total_counts"], 0.99)
sc.pp.filter_cells(adata, min_counts=20)
sc.pp.filter_cells(adata, max_counts=thres)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()
sc.pp.scale(adata, zero_center=False, max_value=10)

sc.pl.highest_expr_genes(adata, n_top=30, )

sc.pp.pca(adata, n_comps=50)
sc.pl.pca_variance_ratio(adata, log=True)
sc.pp.neighbors(adata, metric="cosine")
sc.tl.leiden(adata, flavor="igraph", n_iterations=-1, resolution=0.5)
sc.tl.umap(adata)
sc.pl.umap(adata, color="leiden", legend_loc='on data')

sdata.tables["table"].obs["region"] = "baysor_boundaries"
sdata.set_table_annotates_spatialelement("table", region="baysor_boundaries")

leinden = sdata.tables['table'].obs.leiden
leinden.index= sdata.tables['table'].obs.cell_id
sdata['baysor_boundaries']['leiden'] = sdata.tables['table'].obs.leiden
sdata.pl.render_shapes("baysor_boundaries", color="leiden").pl.show()
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
markers.to_csv('fullmarkers_leiden_res05.csv')
sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key="wilcoxon", groupby="leiden", show_gene_labels=True)
sc.pl.rank_genes_groups_dotplot(adata,n_genes=5, key="wilcoxon", groupby="leiden")

#need to filter bf doing dict
filtered_signature = signature[signature['Vizgen_Gene'].isin(adata.var_names.tolist())]

marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Population']].append(row['Vizgen_Gene'])
marker_gene_dictionary=dict(marker_gene_dictionary)
sc.pl.dotplot(adata, marker_gene_dictionary, 'leiden', dendrogram=True)

sdata.write("tuto.zarr",overwrite=True)

#novae
novae.spatial_neighbors(adata, radius=80)
novae.plot.connectivities(adata)
model = novae.Novae.from_pretrained("MICS-Lab/novae-brain-0")
model.compute_representations(adata, zero_shot=True)
model.fine_tune(adata)
model.compute_representations(adata)
model.assign_domains(adata, level=7)
novae.plot.domains(adata)
novae.plot.spatially_variable_genes(adata, top_k=3, vmax="p95", cell_size=20)
novae.plot.pathway_scores(adata, pathways="/home/marine-louarn/ref/h.all.v2026.1.Hs.json", figsize=(10, 7))

adata.write_h5ad("region_20EN-1063.h5ad")
adata_sc =sc.read_h5ad("region_20EN-1063_scconcept.h5ad")
