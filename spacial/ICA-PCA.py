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
import plotly.graph_objects as go


def ica(adata, n_components, inplace=True, **kwargs): 
    from sklearn.decomposition import FastICA 
    ica_transformer = FastICA(n_components=n_components, **kwargs) 
    x_ica = ica_transformer.fit_transform(adata.X.toarray()) 
    if inplace:
        adata.obsm["X_ica"] = x_ica 
        adata.varm["ICs"] = ica_transformer.components_.T 
    else:
        return ica_transformer 

def ic_tokeep(adata,ics, signature,signature_level=str,topX=2):
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

def ic_tokeep_deg(adata,ics, file_deg,topX=2):
  ic_tokeeps = []
  for cluster in file_deg['cluster'].unique():
      genes = file_deg.loc[file_deg['cluster']==cluster]
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

signature = pd.read_csv("/home/marine-louarn/Documents/test/20260204_Signatures_Cell_populations_HUMAN.csv")
signature_v2 = pd.read_csv("/home/marine-louarn/ref/20251001  JDD_BreastK.csv")
adata = sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1.h5ad")

sc.pp.pca(adata, n_comps=50)
ica(adata,n_components=50)
ics = pd.DataFrame(adata.varm['ICs'])
ics.index = (adata.var['gene_ids']).index


ic_tokeeps = ic_tokeep(adata,ics, signature,"Major_Cell_Populations")
adata.obsm["X_ic_pca" ] = np.concatenate((adata.obsm["X_ica" ][:,ic_tokeeps], adata.obsm["X_pca" ]),axis=1)
adata.varm["IC_PCs"] = np.concatenate((adata.varm["ICs"][:,ic_tokeeps], adata.varm["PCs" ]),axis=1)

sc.pp.neighbors(adata, metric="cosine",use_rep="X_ic_pca",key_added='neigh_ICPC')
sc.tl.leiden(adata, flavor="igraph", n_iterations=-1, resolution=1,neighbors_key='neigh_ICPC',key_added='leiden_ICPC')
sc.tl.umap(adata,key_added = 'umapICPC',neighbors_key='neigh_ICPC', min_dist=0.1)
#colors=['#023fa5','#7d87b9','#11c638','#ef9708','#0fcfc0','#9cded6','#d5eae7','#f3e1eb','#f6c4e1','#f79cd4','#bec1d4','#d6bcc0','#bb7784','#8e063b','#4a6fe3','#8595e1','#b5bbe3','#e6afb9','#e07b91','#d33f6a','#8dd593','#c6dec7','#ead3c6','#f0b98d']
sc.pl.embedding(adata, color="leiden_ICPC",basis='umapICPC',legend_loc='on data')


#####
# Back to spacial
#####
sdata.tables["table"]=adata
sdata.tables["table"].obs["region"] = "cell_boundaries"
sdata.set_table_annotates_spatialelement("table", region="cell_boundaries")

leinden = sdata.tables['table'].obs.leiden_ICPC
leinden.index= sdata.tables['table'].obs.cell_id
sdata['cell_boundaries']['leiden_ICPC'] = sdata.tables['table'].obs.leiden_ICPC
sdata.pl.render_shapes("cell_boundaries", color="leiden_ICPC").pl.show()
plt.show()

#markers
sc.tl.rank_genes_groups(adata, 'leiden_ICPC', method='wilcoxon', key_added = "wilcoxon_ICPC")
markers = sc.get.rank_genes_groups_df(adata,group=None,key='wilcoxon_ICPC')
markers.to_csv('fullmarkers_ICPC_res1.csv')
sc.tl.dendrogram(adata, groupby="leiden_ICPC")
sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key="wilcoxon_ICPC", groupby="leiden_ICPC", show_gene_labels=True)
sc.pl.rank_genes_groups_dotplot(adata,n_genes=5, key="wilcoxon_ICPC", groupby="leiden_ICPC")

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', key_added = "wilcoxon")
sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key="wilcoxon", groupby="leiden", show_gene_labels=True)
sc.pl.rank_genes_groups_dotplot(adata,n_genes=5, key="wilcoxon", groupby="leiden")

test=adata.obs
adata.obs = test.join(adata_myelo_test.obs, rsuffix="_myelo")
test=adata.obsm.copy()
test.join(adata_myelo_test.obsm, rsuffix="_myelo")

adata.obs['novae_domain_nb']= '0'
adata.obs.loc[adata.obs['novae_domains_7']=='D1006', 'novae_domain_nb'] = '1'
adata.obs.loc[adata.obs['novae_domains_7']=='D1014', 'novae_domain_nb'] = '2'
adata.obs.loc[adata.obs['novae_domains_7']=='D1015', 'novae_domain_nb'] = '3'
adata.obs.loc[adata.obs['novae_domains_7']=='D1016', 'novae_domain_nb'] = '4'
adata.obs.loc[adata.obs['novae_domains_7']=='D983', 'novae_domain_nb'] = '5'
adata.obs.loc[adata.obs['novae_domains_7']=='D984', 'novae_domain_nb'] = '6'

adata.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260529_Sample1_ICPC.h5ad")
to_seqgeq_2(adata,header_sg,"/home/marine-louarn/Documents/Xenium_Calderaro/20260529_Sample1_seqgeq_test.txt",genelist)
to_seqgeq_2(adata,header_sg,"/home/marine-louarn/Documents/Xenium_Calderaro/20260529_Sample1_seqgeq_full.txt",adata.var_names.tolist())

#after_all_ref_simple
adata_simpleRef = sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260603_Sample1_all_annot.h5ad")
test=adata_simpleRef.obs.copy()
adata_simpleRef.obs = test.join(adata_myelo_test.obs, rsuffix="_myelo")
adata_simpleRef.obs['novae_domain_nb']= '0'
adata_simpleRef.obs.loc[adata_simpleRef.obs['novae_domains_7']=='D1006', 'novae_domain_nb'] = '1'
adata_simpleRef.obs.loc[adata_simpleRef.obs['novae_domains_7']=='D1014', 'novae_domain_nb'] = '2'
adata_simpleRef.obs.loc[adata_simpleRef.obs['novae_domains_7']=='D1015', 'novae_domain_nb'] = '3'
adata_simpleRef.obs.loc[adata_simpleRef.obs['novae_domains_7']=='D1016', 'novae_domain_nb'] = '4'
adata_simpleRef.obs.loc[adata_simpleRef.obs['novae_domains_7']=='D983', 'novae_domain_nb'] = '5'
adata_simpleRef.obs.loc[adata_simpleRef.obs['novae_domains_7']=='D984', 'novae_domain_nb'] = '6'

adata_simpleRef.obs['CellTypist_v1_nb']= '0'
adata_simpleRef.obs['CellTypist_v2_nb']= '0'
adata_simpleRef.obs['Tangram_nb']= '0'
adata_simpleRef.obs['SingleR_nb']= '0'
j=0
for i in adata_simpleRef.obs['Annotation_CellTypist'].unique().sort_values() :
   adata_simpleRef.obs.loc[adata_simpleRef.obs['Annotation_CellTypist']==i, 'CellTypist_v1_nb'] = str(j)
   adata_simpleRef.obs.loc[adata_simpleRef.obs['Annotation_CellTypist_v2']==i, 'CellTypist_v2_nb'] = str(j)
   adata_simpleRef.obs.loc[adata_simpleRef.obs['Annotation_Tangram']==i, 'Tangram_nb'] = str(j)
   adata_simpleRef.obs.loc[adata_simpleRef.obs['Annotation_SingleR']==i, 'SingleR_nb'] = str(j)
   j+=1

import functions_spatial as func
header_sg="/home/marine-louarn/Documents/Header_SeqGeq.txt"
signature = pd.read_csv("/home/marine-louarn/Documents/test/20260204_Signatures_Cell_populations_HUMAN.csv")
genelist= signature['gene'].tolist()
func.to_seqgeq_2(adata_simpleRef,header_sg,"/home/marine-louarn/Documents/Xenium_Calderaro/20260603_Sample1_seqgeq_refSimple.txt",genelist)

#adata_ct.obs=adata_ct.obs.drop(columns='neighborhood_valid_myelo') #issues with type when writing the file
#adata_ct.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260601_Sample1_ICPC_ct.h5ad")
#adata_ct.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260602_Sample1_ICPC_ct_hcc_dutertre.h5ad")

colors=['#919191','#919191','#919191','#919191','#919191','#ff0000']
sc.pl.embedding(adata_simpleRef, color=["Annotation_CellTypist","Annotation_CellTypist_v2","Annotation_Tangram","Annotation_SingleR"],basis="umapICPC",palette=colors)



#sankey chart
df= adata.obs[['leiden','leiden_ICPC']]
df1 = df.groupby(['leiden','leiden_ICPC']).size().reset_index()
df1.columns = ['source', 'target', 'value']
df1['source']=pd.to_numeric(df1['source'])
df1['target']=pd.to_numeric(df1['target'])+100
df1['value']=pd.to_numeric(df1['value'])
links_dict = df1.to_dict(orient='list')
unique_source_target = list(pd.unique(df1[['source', 'target']].values.ravel('K')))

fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = unique_source_target,
      color = "blue"
      
    ),
    link = dict(
      source = links_dict["source"],
      target = links_dict["target"],
      value = links_dict["value"],
  
  ))])
fig.update_layout(title_text="Sankey Diagram", font_size=10)
fig.show() 

df= adata.obs[['leiden','predicted_labels']]
df1 = df.groupby(['leiden','predicted_labels']).size().reset_index()
df1.columns = ['source', 'target', 'value']
df1['source']=pd.to_numeric(df1['source'])
#df1['target']=df1['target'].map({101:'B',102:'C1Q Macrophage -16',103:'CD16+ Monocytes -1',104:'CD16+ Monocytes -5',105:'CD16- Monocytes -12',106:'CD16- Monocytes -8',107:'DC2/DC3 -14 ',108:'Endothelial ',109:'FTL Macrophage -17',110:'Fibroblast',111:'HES1 Macrophage -2',112:'Hepatocyte',113:'IL1B Monocytes -15',114:'IL4I1 Macrophage -6',115:'ISG Monocytes -4',116:'Macrophage -11',117:'Macrophage -13',118:'Macrophage -7 ',119:'Myeloid ',120:'Proliferating cells -10',121:'T/NK',122:'TREM2 Macrophage -3',123:'Tcell Doublets -9 ',124:'mDC_10',125:'mDC_2 ',126:'mDC_6 ',127:'mDC_7 ',128:'mDC_9 ',129:'pDC '})
df1['value']=pd.to_numeric(df1['value'])

unique_source_target = list(pd.unique(df1[['source', 'target']].values.ravel('K')))
mapping_dict = {k: v for v, k in enumerate(unique_source_target)}
df1['target'] = df1['target'].map(mapping_dict)
links_dict = df1.to_dict(orient='list')

fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = unique_source_target,
      color = "blue"
      
    ),
    link = dict(
      source = links_dict["source"],
      target = links_dict["target"],
      value = links_dict["value"],
  
  ))])
fig.update_layout(title_text="Sankey Diagram: leiden to predicted (CellTypist)", font_size=10)
fig.show() 


df= adata.obs[['leiden_ICPC','predicted_labels']]
df1 = df.groupby(['leiden_ICPC','predicted_labels']).size().reset_index()
df1.columns = ['source', 'target', 'value']
df1['source']=pd.to_numeric(df1['source'])
#df1['target']=df1['target'].map({101:'B',102:'C1Q Macrophage -16',103:'CD16+ Monocytes -1',104:'CD16+ Monocytes -5',105:'CD16- Monocytes -12',106:'CD16- Monocytes -8',107:'DC2/DC3 -14 ',108:'Endothelial ',109:'FTL Macrophage -17',110:'Fibroblast',111:'HES1 Macrophage -2',112:'Hepatocyte',113:'IL1B Monocytes -15',114:'IL4I1 Macrophage -6',115:'ISG Monocytes -4',116:'Macrophage -11',117:'Macrophage -13',118:'Macrophage -7 ',119:'Myeloid ',120:'Proliferating cells -10',121:'T/NK',122:'TREM2 Macrophage -3',123:'Tcell Doublets -9 ',124:'mDC_10',125:'mDC_2 ',126:'mDC_6 ',127:'mDC_7 ',128:'mDC_9 ',129:'pDC '})
df1['value']=pd.to_numeric(df1['value'])

unique_source_target = list(pd.unique(df1[['source', 'target']].values.ravel('K')))
mapping_dict = {k: v for v, k in enumerate(unique_source_target)}
df1['target'] = df1['target'].map(mapping_dict)
links_dict = df1.to_dict(orient='list')

fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = unique_source_target,
      color = "blue"
      
    ),
    link = dict(
      source = links_dict["source"],
      target = links_dict["target"],
      value = links_dict["value"],
  
  ))])
fig.update_layout(title_text="Sankey Diagram: leiden using IC to predicted (CellTypist)", font_size=10)
fig.show() 

