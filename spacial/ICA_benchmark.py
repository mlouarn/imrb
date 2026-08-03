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

sc.pp.pca(adata, n_comps=50)
ica(adata,n_components=100)
ics = pd.DataFrame(adata.varm['ICs'])
ics.index = (adata.var['gene_ids']).index

for top in range(2,5):
    ic_tokeeps = ic_tokeep(adata_myelo_test, signature_v2[signature_v2['Family']=='Myeloid'],"Cell_Subset",topX=top)
    for max_pc in range(10,50,10):
        adata_myelo_test.obsm["X_ic_pca" ] = np.concatenate((adata_myelo_test.obsm["X_ica" ][:,ic_tokeeps], adata_myelo_test.obsm["X_pca" ][:,list(range(0,max_pc))]),axis=1)
        adata_myelo_test.varm["IC_PCs"] = np.concatenate((adata_myelo_test.varm["ICs"][:,ic_tokeeps], adata_myelo_test.varm["PCs" ][:,list(range(0,max_pc))]),axis=1)

        sc.pp.neighbors(adata_myelo_test, metric="cosine",use_rep="X_ic_pca",key_added="neigh_IC"+str(top)+"_PC"+str(max_pc))
        sc.tl.leiden(adata_myelo_test, flavor="igraph", n_iterations=-1, resolution=1.5,neighbors_key="neigh_IC"+str(top)+"_PC"+str(max_pc),key_added="leiden_IC"+str(top)+"_PC"+str(max_pc))
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

#Validation/scoring in sample 1 res1.5
adata_myelo_1= sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260629_Sample1_myelo_noB_ICPC.h5ad")
#sum exp gene of op in cluster
import pyucell as uc

signature_v2 = pd.read_csv("/home/marine-louarn/ref/20251001  JDD_BreastK.csv")
signature_v2 = signature_v2[signature_v2['Family']=='Myeloid']

filtered_signature = signature_v2[signature_v2['gene'].isin(adata_myelo_1.var_names.tolist())]
marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Cell_Subset']].append(row['gene']) #change for sub

marker_gene_dictionary=dict(marker_gene_dictionary)
uc.compute_ucell_scores(adata_myelo_1, signatures=marker_gene_dictionary, chunk_size=500)

a = list(marker_gene_dictionary.keys())
a = [x for x in a if str(x) != 'nan']
a = [x+"_UCell" for x in a]
tmp=adata_myelo_1.obs[a].copy()
tmp['Unsure']=0.001
adata_myelo_1.obs['U_cell_signature']=tmp.idxmax(axis=1)
sc.pl.embedding(adata_myelo_1, color="U_cell_signature",basis='UMAP_IC2_PC10')
sc.pl.embedding(adata_myelo_1, color=['Macro FOLR2_UCell', 'Macro IL4I1_UCell', 'Macro M1_UCell', 'Macro M2_UCell', 'Macro pressure_UCell', 'Macro TREM2_UCell', 'MGC_UCell', 'pan-Macro_UCell', 'CCR7 mDC_UCell', 'DC CD207_UCell', 'DC1_UCell', 'DC2+3_UCell', 'Mono_UCell', 'leiden_res1.5_IC2_PC10'],basis='UMAP_IC2_PC10')
pd.crosstab(adata_myelo_1.obs["MNP_cellType"], adata_myelo_1.obs["U_cell_signature"])

df = pd.DataFrame(columns=a,index=adata_myelo_1.obs['leiden_res1.5_IC2_PC10'].unique())
for j in adata_myelo_1.obs['leiden_res1.5_IC2_PC10'].unique():
    i_index= adata_myelo_1.obs[adata_myelo_1.obs['leiden_res1.5_IC2_PC10']==j]
    mat=i_index[a].copy()
    mean=mat.mean(axis=0)
    df.loc[j]=pd.to_numeric(mean)


df = df[['Macro FOLR2_UCell', 'Macro IL4I1_UCell', 'Macro TREM2_UCell', 'CCR7 mDC_UCell', 'DC CD207_UCell', 'DC1_UCell', 'DC2+3_UCell', 'Mono_UCell']].apply(pd.to_numeric)
df=df.T
df_norm = df.div(df.sum(axis=1), axis=0)

sns.heatmap(df_norm,cmap="crest")
plt.show()



#compared to Seqgeq signature
mean_signature = pd.read_csv("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1/SyntParam.csv", skiprows=5,index_col='Gene')
mean_signature=mean_signature.loc[['DC1Mean','MacroFOLR2Mean','MacroIL4I1Mean','MacroTREM2Mean','Mast cellsMean','mregDCMean','pan-BMean','pan-KupfferCellsMean','PMNMean','cMoMean']]
#mean_signature.columns =mean_signature.columns[1:].append(pd.Index(['aaaajiih-1']))
mean_signature.columns=pd.Index(['aaaajiih-1']).append(mean_signature.columns[1:])
df = pd.DataFrame(columns=mean_signature.index,index=adata_myelo_1.obs['leiden_res1.5_IC2_PC10'].unique())
for j in adata_myelo_1.obs['leiden_res1.5_IC2_PC10'].unique():
    i_index= adata_myelo_1.obs[adata_myelo_1.obs['leiden_res1.5_IC2_PC10']==j]
    tmp = mean_signature[i_index['cell_id']]
    tmp=tmp.T
    mean=tmp.mean()
    df.loc[j]=pd.to_numeric(mean)

df[mean_signature.index] = df[mean_signature.index].apply(pd.to_numeric)

sns.heatmap(df)
plt.show()

#pd.crosstab(adata_myelo_1.obs["MNP_cellType"], adata_myelo_1.obs["leiden_IC2_PC10"])

#distance umap
import numpy as np
from scipy.spatial.distance import pdist
dist = pd.DataFrame(columns=adata_myelo_test.obs['MNP_cellType'].unique(),index=list(adata_myelo_test.obsm.keys()))
j=0
for k in list(adata_myelo_test.obsm.keys()):
	for i in adata_myelo_test.obs['MNP_cellType'].unique():
		tmp = adata_myelo_test[adata_myelo_test.obs['MNP_cellType']==i].copy()
		dist[i].values[j] = pdist(tmp.obsm[k]).mean() 
	j+=1


#jaccard annotated

df = pd.DataFrame(columns=adata_myelo_1.obs['MNP_cellType'].unique(),index=adata_myelo_1.obs['leiden_res1.5_IC2_PC10'].unique())
for i in adata_myelo_1.obs['MNP_cellType'].unique():
    row=[]
    for j in adata_myelo_1.obs['leiden_res1.5_IC2_PC10'].unique():
        i_index= adata_myelo_1.obs[adata_myelo_1.obs['MNP_cellType']==i].index
        j_index= adata_myelo_1.obs[adata_myelo_1.obs['leiden_res1.5_IC2_PC10']==j].index
        intersection = len(i_index.intersection(j_index))
        union = len(i_index.union(j_index))
        row.append(intersection/union)
    df[i]=row


def score_jaccard_mat(df):
    count_i = 0
    count_j = 0
    for i in df.index:
        sum_i = (df.iloc[[i]]>0.3).sum()
        count_i = count_i+(sum_i==1).sum()
    for j in df.columns:
        sum_j=(df[j]>0.3).sum()
        count_j = count_j+(sum_j>=1).sum()
    score = (count_i+count_j)/(len(df.index)+len(df.columns))
    return score


sns.heatmap(df)
plt.show()
