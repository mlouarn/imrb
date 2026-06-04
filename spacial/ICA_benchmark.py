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
import pyucell as uc
from collections import defaultdict

adata_myelo_test= sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260529_Sample1_myelo_ICPC.h5ad")

def ica(adata, n_components, inplace=True, **kwargs): 
    from sklearn.decomposition import FastICA 
    ica_transformer = FastICA(n_components=n_components, **kwargs) 
    x_ica = ica_transformer.fit_transform(adata.X.toarray()) 
    if inplace:
        adata.obsm["X_ica"] = x_ica 
        adata.varm["ICs"] = ica_transformer.components_.T 
    else:
        return ica_transformer 
    
def ic_tokeep(adata, signature,signature_level=str,topX=2):
  ic_tokeeps = []
  for cluster in signature[signature_level].unique():
      genes = signature.loc[signature[signature_level]==cluster]
      genes = genes['gene'].unique().tolist()
      genes_keeps = list(set(genes) & set(ics.index))
      top_ics=[]
      if genes_keeps!=[]:
        ics_genes = ics.loc[genes_keeps]
        sum_ics = ics_genes.sum().abs().tolist()
        top_ics = sorted(range(len(sum_ics)), key=lambda i: sum_ics[i])[-topX:] 
      ic_tokeeps.append(top_ics)
  ic_tokeeps = sum(ic_tokeeps,[])
  ic_tokeeps = list(set(ic_tokeeps))
  return ic_tokeeps

mnp_markers = pd.read_csv("/home/marine-louarn/ref/MNP_Verse_DEG_megaclusters.csv")
mnp_adata = sc.read_h5ad("/home/marine-louarn/ref/2021_MNP_Verse.h5ad")
signature = pd.read_csv("/home/marine-louarn/Documents/test/20260204_Signatures_Cell_populations_HUMAN.csv")
signature_v2 = pd.read_csv("/home/marine-louarn/ref/20251001  JDD_BreastK.csv")

adata_myelo_test=adata_myelo.copy()
top=2
ic_tokeeps = ic_tokeep(adata_myelo_test, signature_v2[signature_v2['Family']=='Myeloid'],"Cell_Subset",topX=top)
max_pc=10
adata_myelo_test.obsm["X_ic_pca" ] = np.concatenate((adata_myelo_test.obsm["X_ica" ][:,ic_tokeeps], adata_myelo_test.obsm["X_pca" ][:,list(range(0,max_pc))]),axis=1)
adata_myelo_test.varm["IC_PCs"] = np.concatenate((adata_myelo_test.varm["ICs"][:,ic_tokeeps], adata_myelo_test.varm["PCs" ][:,list(range(0,max_pc))]),axis=1)

sc.pp.neighbors(adata_myelo_test, metric="cosine",use_rep="X_ic_pca",key_added="neigh_IC"+str(top)+"_PC"+str(max_pc))
sc.tl.leiden(adata_myelo_test, flavor="igraph", n_iterations=-1, resolution=1.5,neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc),key_added="leiden_res1.5_IC"+str(top)+"_PC"+str(max_pc))
sc.tl.umap(adata_myelo_test,key_added = "UMAP_IC"+str(top)+"_PC"+str(max_pc),neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc), min_dist=0.1)
#sc.pl.embedding(adata_myelo_test, color="leiden_IC"+str(top)+"_PC"+str(max_pc),basis="UMAP_IC"+str(top)+"_PC"+str(max_pc))
sc.pl.embedding(adata_myelo_test, color="leiden_ICPC",basis="umapICPC")
sc.pl.embedding(adata_myelo, color="leiden_ICPC",basis="umapICPC")

sc.pl.embedding(adata_myelo_test, color=["leiden_res1.5_IC"+str(top)+"_PC"+str(max_pc)],basis="UMAP_IC"+str(top)+"_PC"+str(max_pc))

adata_myelo_test.obs['novae_domain_nb']= '0'
adata_myelo_test.obs.loc[adata_myelo_test.obs['novae_domains_7']=='D1006', 'novae_domain_nb'] = '1'
adata_myelo_test.obs.loc[adata_myelo_test.obs['novae_domains_7']=='D1014', 'novae_domain_nb'] = '2'
adata_myelo_test.obs.loc[adata_myelo_test.obs['novae_domains_7']=='D1015', 'novae_domain_nb'] = '3'
adata_myelo_test.obs.loc[adata_myelo_test.obs['novae_domains_7']=='D1016', 'novae_domain_nb'] = '4'
adata_myelo_test.obs.loc[adata_myelo_test.obs['novae_domains_7']=='D983', 'novae_domain_nb'] = '5'
adata_myelo_test.obs.loc[adata_myelo_test.obs['novae_domains_7']=='D984', 'novae_domain_nb'] = '6'

to_seqgeq_2(adata_myelo_test,header_sg,"/home/marine-louarn/Documents/Xenium_Calderaro/20260529_Sample1_onlyMyelo_seqgeq_full.txt",adata.var_names.tolist())
adata_myelo_test.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260529_Sample1_myelo_ICPC.h5ad")

def to_seqgeq_2(adata, header_sg=str,file=str,genelist=list):
    with open(file, "w") as f:
        # write header
        with open(header_sg) as header:
            f.write(header.read())
        # write barcodes
        f.write("\t".join(adata.obs.cell_id)+"\n")
        # keep only genes in signatures
        genes_absent = [gene for gene in genelist if gene not in adata.var_names]
        genes_present = [gene for gene in genelist if gene in adata.var_names]
        genes_index = [adata.var_names.get_loc(gene) for gene in genes_present]
        print(f'genes absent from adata {genes_absent}')

        # write normalized counts
        mat = adata.layers["scaled"].T # gene x cell change to "counts" 
        for index in genes_index:
            gene_name = adata.var_names[index]
            row = mat[index,:].toarray()[0]
            row_str = [str(x) for x in row]
            f.write(gene_name + "\t" + "\t".join(row_str)+"\n")

        # write UMAP
        umap_x = adata.obsm['X_umap'][:,0]
        umap_x_str = [str(x) for x in umap_x]
        umap_y = adata.obsm['X_umap'][:,1]
        umap_y_str = [str(y) for y in umap_y]
        umap_icpc_x = adata.obsm['umapICPC'][:,0]
        umap_icpc_x_str = [str(x) for x in umap_icpc_x]
        umap_icpc_y = adata.obsm['umapICPC'][:,1]
        umap_icpc_y_str = [str(y) for y in umap_icpc_y]
        #umap_icpc2_10_x = adata.obsm['UMAP_IC2_PC10'][:,0]
        #umap_icpc2_10_x_str = [str(x) for x in umap_icpc2_10_x]
        #umap_icpc2_10_y = adata.obsm['UMAP_IC2_PC10'][:,1]
        #umap_icpc2_10_y_str = [str(y) for y in umap_icpc2_10_y]
        spatial_x = adata.obsm['spatial'][:,0]
        spatial_x_str = [str(x) for x in spatial_x]
        spatial_y = adata.obsm['spatial'][:,1]
        spatial_y_str = [str(y) for y in spatial_y]
        n_counts = adata.obs['n_counts']
        n_counts_str = [str(y) for y in n_counts]
        leiden_myelo = adata.obs['leiden_res1.5_IC2_PC10']#_myelo']
        leiden_myelo_str = [str(y) for y in leiden_myelo]
        f.write("umap_PC1" + "\t" + "\t".join(umap_x_str)+"\n")
        f.write("umap_PC2" + "\t" + "\t".join(umap_y_str)+"\n")
        f.write("umap_ICPC1" + "\t" + "\t".join(umap_icpc_x_str)+"\n")
        f.write("umap_ICPC2" + "\t" + "\t".join(umap_icpc_y_str)+"\n")
        #f.write("umap_IC2PC10_1" + "\t" + "\t".join(umap_icpc2_10_x_str)+"\n")
        #f.write("umap_IC2PC10_2" + "\t" + "\t".join(umap_icpc2_10_y_str)+"\n")
        f.write("Spatial_X" + "\t" + "\t".join(spatial_x_str)+"\n")
        f.write("Spatial_Y" + "\t" + "\t".join(spatial_y_str)+"\n")
        f.write("leiden_clusters_PC" + "\t" + "\t".join(adata.obs['leiden']) + "\n")
        f.write("leiden_clusters_ICPC" + "\t" + "\t".join(adata.obs['leiden_ICPC']) + "\n")
        f.write("leiden_res1.5_IC2PC10" + "\t" + "\t".join(leiden_myelo_str) + "\n")
        f.write("nCount" + "\t" + "\t".join(n_counts_str) + "\n")
        f.write("novae_domain_nb" + "\t" + "\t".join(adata.obs['novae_domain_nb']) + "\n")
        f.write("CellTypist_nb" + "\t" + "\t".join(adata.obs['CellTypist_nb']) + "\n")


#pseudobulk
deg_mnp = pd.read_csv("/home/marine-louarn/ref/MNP_Verse_DEG_megaclusters.csv")
sc.tl.rank_genes_groups(adata_myelo_test, 'leiden_res1.5_IC2_PC10', method='wilcoxon', key_added = "wilcoxon")
markers_myelo = sc.get.rank_genes_groups_df(adata_myelo_test,group=None,key='wilcoxon')
markers_myelo_pval = markers_myelo[markers_myelo['pvals']<0.05]
deg_mnp_pval = deg_mnp[deg_mnp['p_val']<0.05]

myelo_pseudobulk = sc.get.aggregate(adata_myelo_test, by=["leiden_res1.5_IC2_PC10"], func="sum", layer="counts")
ref_mnp = sc.read_h5ad("/home/marine-louarn/ref/2021_MNP_Verse.h5ad")
mnp_pseudobulk = sc.get.aggregate(ref_mnp, by=["MegaCluster"], func="sum", layer="counts")

genes_to_keep = list(set(deg_mnp['gene']) & set(markers_myelo['names']))
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
sns.heatmap(cor, annot=True,cmap="crest")
plt.show()

matrix = np.zeros((len(mnp_pseudobulk_mat.columns), len(myelo_pseudobulk_mat.columns)))
matrix=pd.DataFrame(matrix)
matrix.index=mnp_pseudobulk_mat.columns
matrix.columns=myelo_pseudobulk_mat.columns
for i in mnp_pseudobulk_mat.columns :
    for j in myelo_pseudobulk_mat.columns :
        mnp_gene=deg_mnp_pval[deg_mnp_pval['cluster']==i]['gene']
        myelo_gene=markers_myelo_pval[markers_myelo_pval['group']==j]['names']
        overlap = set(mnp_gene) & set(myelo_gene)
        matrix[j][i]= len(overlap)

sns.heatmap(matrix, annot=True,cmap="crest")
plt.show()