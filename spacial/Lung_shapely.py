import spatialdata as sd
import scanpy as sc
import os
from shapely import box, Point, STRtree
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns 
import squidpy as sq

'''
get c_id of neighbour
'''
sdata_p1 = sd.read_zarr("/home/marine-louarn/Documents/Lung_Spatial/Patient1_proseg.zarr") 
adata_patient1= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20261015_Patient1_withAllSubtypes.h5ad")

geom_p1 = sdata_p1.shapes['proseg_boundaries']['geometry']
tree_p1 = STRtree(geom_p1)
adata_patient1.obs['Neighbour_shapely'] = 'nan'
adata_patient1.obs['Neighbour_shapely'] = adata_patient1.obs['Neighbour_shapely'].astype(object)
for cell in adata_patient1.obs['cell_id'].tolist():
    indices = tree_p1.query(geom_p1[cell],predicate='intersects').tolist()
    cell_touching = geom_p1.index[indices].tolist()
    cell_touching.remove(cell)
    adata_patient1.obs.at[cell, 'Neighbour_shapely'] = cell_touching

sdata_p2 = sd.read_zarr("/home/marine-louarn/Documents/Lung_Spatial/Patient2_proseg.zarr") 
adata_p2 = sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20260922_Patient2_cellposeProseg.h5ad")

geom_p2 = sdata_p2.shapes['proseg_boundaries']['geometry']
tree_p2 = STRtree(geom_p2)
adata_p2.obs['Neighbour_shapely'] = 'nan'
adata_p2.obs['Neighbour_shapely'] = adata_p2.obs['Neighbour_shapely'].astype(object)
for cell in adata_p2.obs['cell_id'].tolist():
    indices = tree_p2.query(geom_p2[cell],predicate='intersects').tolist()
    cell_touching = geom_p2.index[indices].tolist()
    cell_touching.remove(cell)
    adata_p2.obs.at[cell, 'Neighbour_shapely'] = cell_touching

'''
T next to plasma or B
'''
#patient1

list_plasma_cell = adata_patient1.obs[adata_patient1.obs['Cell_Population'].isin(['Plasma','B_Plasma','Plasma.CXCL10'])].index.tolist()
adata_B= sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20261008_Patient1_Bplasma.h5ad")
list_BandPlasma = adata_B.obs.index.tolist()
list_B = adata_patient1.obs[adata_patient1.obs['Cell_Population'].isin(['B_stroma','B.CD83','B_endo','B_macro','B.prolf','B_perivascular'])].index.tolist()

adata_patient1.obs['Voisin_Plasma']= '0'
adata_patient1.obs['Voisin_B'] = '0'
adata_patient1.obs['Voisin_B_and_Plasma'] = '0'

for cell in adata_patient1.obs['cell_id'].tolist():
    val_plasma = len(set(adata_patient1.obs['Neighbour_shapely'][cell]).intersection(list_plasma_cell))
    val_B = len(set(adata_patient1.obs['Neighbour_shapely'][cell]).intersection(list_B))
    val_Bandplasma = len(set(adata_patient1.obs['Neighbour_shapely'][cell]).intersection(list_BandPlasma))
    adata_patient1.obs.at[cell, 'Voisin_Plasma'] = str(val_plasma)
    adata_patient1.obs.at[cell, 'Voisin_B'] = str(val_B)
    adata_patient1.obs.at[cell, 'Voisin_B_and_Plasma'] = str(val_Bandplasma)

adata_patient1.obs['Voisin_B_YN'] = ['Y' if int(i) >0 else 'N' for i in adata_patient1.obs['Voisin_B']]
adata_patient1.obs['Voisin_Plasma_YN'] = ['Y' if int(i) >0 else 'N' for i in adata_patient1.obs['Voisin_Plasma']]
adata_patient1.obs['Voisin_B_and_Plasma_YN'] = ['Y' if int(i) >0 else 'N' for i in adata_patient1.obs['Voisin_B_and_Plasma']]

adata_patient1.obs['Voisin'] ='None'
adata_patient1.obs.loc[adata_patient1.obs['Voisin_B_YN']== 'Y' , 'Voisin'] = 'B'
adata_patient1.obs.loc[adata_patient1.obs['Voisin_Plasma_YN']== 'Y' , 'Voisin'] = 'Plasma'
adata_patient1.obs.loc[(adata_patient1.obs['Voisin_Plasma_YN']== 'Y') & (adata_patient1.obs['Voisin_B_YN']== 'Y'), 'Voisin'] = 'Both'

df=adata_patient1.obs[['cell', 'Cell_Type', 'Cell_Population', 'Neighbour_shapely', 'Voisin_Plasma', 'Voisin_B', 'Voisin_B_and_Plasma', 'Voisin_B_YN', 'Voisin_Plasma_YN', 'Voisin_B_and_Plasma_YN', 'Voisin']].copy()
df.to_csv("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20261022_Patient1_Bneighbour_meta.csv")

with open('/home/marine-louarn/Documents/Lung_Spatial/Patient1/20261019_ALL_neighbour_B.txt', "w") as f:
	f.write("Voisin_Plasma" + "\t" + "\t".join(adata_patient1.obs['Voisin_Plasma']) + "\n")
	f.write("Voisin_B" + "\t" + "\t".join(adata_patient1.obs['Voisin_B']) + "\n")
	f.write("Voisin_B_and_Plasma" + "\t" + "\t".join(adata_patient1.obs['Voisin_B_and_Plasma']) + "\n")
	f.write("Subset_leiden" + "\t" + "\t".join(adata_patient1.obs['Subset_leiden'].astype(str)) + "\n")
	f.write("Subset_umap_1" + "\t" + "\t".join(adata_patient1.obs['Subset_umap_1'].astype(str)) + "\n")
	f.write("Subset_umap_2" + "\t" + "\t".join(adata_patient1.obs['Subset_umap_2'].astype(str)) + "\n")
	f.write("Subset_leiden_IC_PC" + "\t" + "\t".join(adata_patient1.obs['Subset_leiden_IC_PC'].astype(str)) + "\n")
	f.write("Subset_umap_ICPC_1" + "\t" + "\t".join(adata_patient1.obs['Subset_umap_ICPC_1'].astype(str)) + "\n")
	f.write("Subset_umap_ICPC_2" + "\t" + "\t".join(adata_patient1.obs['Subset_umap_ICPC_2'].astype(str)) + "\n")


adata_p1CD8=adata_patient1[adata_patient1.obs['Cell_Population']=='T.CD8']
sc.pl.dotplot(adata_p1CD8, ['TNFRSF17','STAP1','CD79A','MZB1','IFNAR2','CXCR6','GPX3','ANGPT1','CXCL12','PDGFRA','IL13','EOMES','CD79B','MMP2','HLA-DQA1','IL3RA','VSIR','TNFRSF13B','CD27','SNAI1','FOXP3','CMKLR1'], 
              groupby='Voisin_Plasma_YN',show=False,save='Patient1_BorPlasma.pdf')



#patient2

list_plasma_cell = adata_p2.obs[adata_p2.obs['Cell_Type'].isin(['Plasma cells'])].index.tolist()
list_B = adata_p2.obs[adata_p2.obs['Cell_Type'].isin(['B cells'])].index.tolist()
list_BandPlasma = adata_p2.obs[adata_p2.obs['Cell_Type'].isin(['B cells','Plasma cells'])].index.tolist()

adata_p2.obs['Voisin_Plasma']= '0'
adata_p2.obs['Voisin_B'] = '0'
adata_p2.obs['Voisin_B_and_Plasma'] = '0'
adata_p2.obs['Voisin_Plasma_YN']= 'N'
adata_p2.obs['Voisin_B_YN'] = 'N'
adata_p2.obs['Voisin_B_and_Plasma_YN'] = 'N'

for cell in adata_p2.obs['cell_id'].tolist():
    val_plasma = len(set(adata_p2.obs['Neighbour_shapely'][cell]).intersection(list_plasma_cell))
    val_B = len(set(adata_p2.obs['Neighbour_shapely'][cell]).intersection(list_B))
    val_Bandplasma = len(set(adata_p2.obs['Neighbour_shapely'][cell]).intersection(list_BandPlasma))
    adata_p2.obs.at[cell, 'Voisin_Plasma'] = str(val_plasma)
    adata_p2.obs.at[cell, 'Voisin_Plasma_YN'] = 'Y' if val_plasma > 0 else 'N'
    adata_p2.obs.at[cell, 'Voisin_B'] = str(val_B)
    adata_p2.obs.at[cell, 'Voisin_B_YN'] = 'Y' if val_B > 0 else 'N'
    adata_p2.obs.at[cell, 'Voisin_B_and_Plasma'] = str(val_Bandplasma)
    adata_p2.obs.at[cell, 'Voisin_B_and_Plasma_YN'] = 'Y' if val_Bandplasma > 0 else 'N'


adata_p2.obs['Voisin'] ='None'
adata_p2.obs.loc[adata_p2.obs['Voisin_B_YN']== 'Y' , 'Voisin'] = 'B'
adata_p2.obs.loc[adata_p2.obs['Voisin_Plasma_YN']== 'Y' , 'Voisin'] = 'Plasma'
adata_p2.obs.loc[(adata_p2.obs['Voisin_Plasma_YN']== 'Y') & (adata_p2.obs['Voisin_B_YN']== 'Y'), 'Voisin'] = 'Both'

adata_p2.obs['Voisin_nb'] ='0'
adata_p2.obs.loc[adata_p2.obs['Voisin_B_YN']== 'Y' , 'Voisin_nb'] = '1'
adata_p2.obs.loc[adata_p2.obs['Voisin_Plasma_YN']== 'Y' , 'Voisin_nb'] = '2'
adata_p2.obs.loc[(adata_p2.obs['Voisin_Plasma_YN']== 'Y') & (adata_p2.obs['Voisin_B_YN']== 'Y'), 'Voisin_nb'] = '3'

df=adata_p2.obs[['cell', 'Cell_Type', 'Cell_Population', 'Neighbour_shapely', 'Voisin_Plasma', 'Voisin_B', 'Voisin_B_and_Plasma', 'Voisin_B_YN', 'Voisin_Plasma_YN', 'Voisin_B_and_Plasma_YN', 'Voisin']].copy()
df.to_csv("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261022_Patient2_Bneighbour_meta.csv")


adata_p2.obs['Is_T.CD8'] = ['1' if i=='T.CD8' else '0' for i in adata_p2.obs['Cell_Population']]

adata_p2T = adata_p2[adata_p2.obs['Cell_Population']=='T.CD8'].copy()
#sc.tl.rank_genes_groups(adata_p2T, groupby='Voisin_Plasma_YN', key_added='group_DE_results')
#marker_B_neigh = sc.get.rank_genes_groups_df(adata_p2T,group=None,key='group_DE_results')
#marker_B_neigh.to_csv('/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261020_marker_CD8voisinPlasma.csv')
#sc.pl.rank_genes_groups_heatmap(adata_p2T, groupby='Voisin_Plasma_YN',key='group_DE_results')

sc.pl.dotplot(adata_p2T, ['TNFRSF17','STAP1','CD79A','MZB1','IFNAR2','CXCR6','GPX3','ANGPT1','CXCL12','PDGFRA','IL13','EOMES','CD79B','MMP2','HLA-DQA1','IL3RA','VSIR','TNFRSF13B','CD27','SNAI1','FOXP3','CMKLR1'], 
              groupby='Voisin_Plasma_YN',show=False,save='Patient2_BorPlasma.png')

with open('/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261021_ALL_neighbour_B.txt', "w") as f:
    f.write("Is_T.CD8_0.No_1.Yes" + "\t" + "\t".join(adata_p2.obs['Is_T.CD8']) + "\n")
    f.write("Voisin_Plasma" + "\t" + "\t".join(adata_p2.obs['Voisin_Plasma']) + "\n")
    f.write("Voisin_B" + "\t" + "\t".join(adata_p2.obs['Voisin_B']) + "\n")
    f.write("Voisin_B_and_Plasma" + "\t" + "\t".join(adata_p2.obs['Voisin_B_and_Plasma']) + "\n")
    f.write("Voisin_0.None_1.B_2.Plasma_3.both" + "\t" + "\t".join(adata_p2.obs['Voisin_nb']) + "\n")
    f.write("Subset_leiden" + "\t" + "\t".join(adata_p2.obs['Subset_leiden'].astype(str)) + "\n")
    f.write("Subset_umap_1" + "\t" + "\t".join(adata_p2.obs['Subset_umap_1'].astype(str)) + "\n") 
    f.write("Subset_umap_2" + "\t" + "\t".join(adata_p2.obs['Subset_umap_2'].astype(str)) + "\n")
    f.write("Subset_leiden_IC_PC" + "\t" + "\t".join(adata_p2.obs['Subset_leiden_IC_PC'].astype(str)) + "\n")
    f.write("Subset_umap_ICPC_1" + "\t" + "\t".join(adata_p2.obs['Subset_umap_ICPC_1'].astype(str)) + "\n")
    f.write("Subset_umap_ICPC_2" + "\t" + "\t".join(adata_p2.obs['Subset_umap_ICPC_2'].astype(str)) + "\n")

'''
Keeping all cells at X of T in contact of plasma
'''
x = 2
#p1
sq.gr.spatial_neighbors_radius(adata_patient1, radius=x)
sq.tl.var_by_distance(adata_patient1,cluster_key="Voisin_Plasma_YN",groups="Y",spatial_key="spatial")
distances = adata_patient1.obsm["design_matrix"]["Y_raw"]
adata_subset_p1 = adata_patient1[distances <= x].copy()

sc.pl.embedding(adata_subset_p1, color="Cell_Type",basis='spatial')

adata_patient1.obs['Domain_Tplasma'] = 'N'
adata_patient1.obs.loc[adata_patient1.obs['cell_id'].isin(adata_subset_p1.obs['cell_id']), 'Domain_Tplasma'] = 'Y'

df=adata_patient1.obs[['Cell_Population','Domain_Tplasma']]
df_new = df.pivot_table(columns="Cell_Population", index="Domain_Tplasma", aggfunc=len, fill_value=0)
ax = df_new.plot.bar(stacked=True,cmap='tab20')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Type", framealpha=0);
plt.show()

#p2
sq.gr.spatial_neighbors_radius(adata_p2, radius=x)
sq.tl.var_by_distance(adata_p2,cluster_key="Voisin_Plasma_YN",groups="Y",spatial_key="spatial")
distances = adata_p2.obsm["design_matrix"]["Y_raw"]
adata_subset = adata_p2[distances <= x].copy()

sc.pl.embedding(adata_subset, color="Cell_Type",basis='spatial')

adata_p2.obs['Domain_Tplasma'] = 'N'
adata_p2.obs.loc[adata_p2.obs['cell_id'].isin(adata_subset.obs['cell_id']), 'Domain_Tplasma'] = 'Y'

df=adata_p2.obs[['Cell_Subset','Domain_Tplasma']]
df_new = df.pivot_table(columns="Cell_Subset", index="Domain_Tplasma", aggfunc=len, fill_value=0)
ax = df_new.plot.bar(stacked=True,cmap='gist_rainbow')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Type", framealpha=0);
plt.show()
df_per = df.pivot_table(columns="Cell_Subset", index="Domain_Tplasma", aggfunc=len, fill_value=0).apply(lambda x: x/sum(x), axis=1)
ax = df_per.plot.bar(stacked=True,cmap='gist_rainbow')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Population", framealpha=0);
plt.show()


filter_col = [col for col in df_new if col.startswith('Fibroblasts')]
df_fibro = df_new[filter_col]
ax = df_fibro.plot.bar(stacked=True,cmap='tab20b')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Type", framealpha=0);
plt.show()

df_per_fibro = df_fibro.apply(lambda x: x/sum(x), axis=1)
ax = df_per_fibro.plot.bar(stacked=True,cmap='gist_rainbow')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Population", framealpha=0);
plt.show()

df_per_fibro_T = df_per_fibro.T
df_per_fibro_T['fold']=df_per_fibro_T['Y']/df_per_fibro_T['N']

cxcl12= list(pd.read_csv("Patient2/leiden_clusters_Fibro_CXCL12+.csv", skiprows=5).columns)

adata_p2.obs['Cell_Subset'] = adata_p2.obs['Cell_Subset'].astype('str')
adata_p2.obs.loc[adata_p2.obs['cell_id'].isin(cxcl12), 'Cell_Subset'] = 'Fibroblasts.CXCL12'

'''
FOLR2 environment
'''

adata_p1_folr2 = adata_patient1[adata_patient1.obs['Cell_Population']=='Mac.FOLR2'].copy()

p1_folrd_df = pd.DataFrame(index = adata_p1_folr2.obs['cell_id'].tolist(),columns =adata_patient1.obs['Cell_Type'].cat.categories.tolist())
p1_folrd_df = p1_folrd_df.replace(np.NaN, 0)
#p1_folrd_df.fillna(0, inplace=True)

for cell in adata_p1_folr2.obs['cell_id'].tolist():
    voisins = adata_p1_folr2.obs['Neighbour_shapely'][cell]
    for i in voisins :
        if i in adata_patient1.obs['cell_id'].tolist() :
            ct = adata_patient1.obs['Cell_Type'][i]
            p1_folrd_df[ct].loc[cell]
            p1_folrd_df.loc[cell, ct] = p1_folrd_df[ct].loc[cell]+1


sns.clustermap(p1_folrd_df)
plt.show()