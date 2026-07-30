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
adata_s2.obs.insert(0, 'Sample', 'Sample2')
adata_s1=sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260603_Sample1_all_annot.h5ad")
adata_s1.obs.insert(0, 'Sample', 'Sample1')

adata_myelo_2= sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260616_Sample2myelo_ICPC.h5ad")
adata_myelo_1= sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260529_Sample1_myelo_ICPC.h5ad")

#extract
func.h5ad_to_tsv(adata_s1,"/home/marine-louarn/Documents/Xenium_Calderaro/Xenium_exp1/Sample1")

adata_concat = anndata.concat([adata_s1, adata_s2], label="Sample", keys=["Sample1", "Sample2"])
adata_concat.obs_names_make_unique()

'''#bbknn: no good batch correction
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key="Sample") 
sc.tl.umap(adata_concat,key_added = 'bbknn')
sc.pl.embedding(adata_concat, color=["Sample"],basis='bbknn')
'''
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

harmony_int2 = integrate_harmony(adata_concat,'Sample')

sc.pl.embedding(harmony_int2, color=["leiden_harmony","Sample"],basis="harmony")

harmony_int2.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260616_harmonyInts1s2.h5ad")
harmony_int = sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260616_harmonyInts1s2.h5ad")


#ingest see to change umap
overlap = list(set(adata_myelo_1.var_names.tolist())&set(adata_myelo_2.var_names.tolist()))
adata_ingest_s2=adata_myelo_2[:,overlap].copy()
adata_ingest_s2.obsm['X_umap']=adata_ingest_s2.obsm['UMAP_IC2_PC10'].copy()
adata_ingest_s2.uns['umap']=adata_ingest_s2.uns['UMAP_IC2_PC10'].copy()
adata_ingest_s1=adata_myelo_1[:,overlap].copy()
ingest = sc.tl.ingest(adata_ingest_s1, adata_ingest_s2, obs="leiden_IC2_PC10",inplace=False)
ingest.uns["leiden_IC2_PC10_colors"] = adata_ingest_s2.uns["leiden_IC2_PC10_colors"]  # fix colors
sc.pl.umap(ingest, color=["leiden_IC2_PC10"], wspace=0.5)
sc.pl.umap(adata_ingest_s2, color=["leiden_IC2_PC10"], wspace=0.5)
sc.pl.embedding(adata_ingest_s2, color=["leiden_IC2_PC10"],basis='UMAP_IC2_PC10')


s1s2_concat = anndata.concat([adata_ingest_s2, ingest], label="sample", keys=["sample2", "sample1"])
s1s2_concat.obs["leiden_IC2_PC10"] = (
    s1s2_concat.obs["leiden_IC2_PC10"].astype("category").cat.reorder_categories(adata_ingest_s2.obs["leiden_IC2_PC10"].cat.categories)
)
sc.pl.umap(s1s2_concat, color=["sample", "leiden_IC2_PC10"])
sc.pl.embedding(s1s2_concat, color=["sample", "leiden_IC2_PC10"],basis='UMAP_IC2_PC10')


s1s2_concat.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260617_Ingests1s2_myelo.h5ad")
s1s2_concat=sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260617_Ingests1s2_myelo.h5ad")

s1s2_concat.obs['Sample_nb']= 'unassigned'
s1s2_concat.obs.loc[s1s2_concat.obs['sample'] =='sample1', 'Sample_nb'] = '1'
s1s2_concat.obs.loc[s1s2_concat.obs['sample']=='sample2', 'Sample_nb'] = '2'

signature_v2 = pd.read_csv("/home/marine-louarn/ref/20251001  JDD_BreastK.csv")
list_tokeep= signature_v2['gene'].tolist()
func.to_seqgeq_ingest(s1s2_concat,header_sg,"/home/marine-louarn/Documents/Xenium_Calderaro/20260617_IntegrationIngest_s1s2_myelo_seqgeq.txt",list_tokeep,True)
func.to_seqgeq_ingest(s1s2_concat,header_sg,"/home/marine-louarn/Documents/Xenium_Calderaro/20260617__IntegrationIngest_s1s2_myelo_seqgeq_full_v2.txt",s1s2_concat.var_names.tolist(),True)
