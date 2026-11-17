import spatialdata as sd
import scanpy as sc
import os
from shapely import box, Point, STRtree
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns 
import squidpy as sq
import pandas as pd
from sklearn.cluster import KMeans

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
adata_p2 = sc.read_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient2/20261103_Patient2_withAllSubtypes.h5ad")

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
FOLR2 environment
'''
#voisin direct
adata_p1_folr2 = adata_patient1[adata_patient1.obs['Cell_Population']=='Mac.FOLR2'].copy()

p1_folrd_df = pd.DataFrame(index = adata_p1_folr2.obs['cell_id'].tolist(),columns =adata_patient1.obs['Cell_Type'].cat.categories.tolist())
p1_folrd_df = p1_folrd_df.replace(np.NaN, 0)
#p1_folrd_df.fillna(0, inplace=True)

for cell in adata_p1_folr2.obs['cell_id'].tolist():
    voisins = adata_p1_folr2.obs['Neighbour_shapely'][cell]
    for i in voisins :
        if i in adata_patient1.obs['cell_id'].tolist() :
            ct = adata_patient1.obs['Cell_Type'][i]
            #p1_folrd_df[ct].loc[cell]
            p1_folrd_df.loc[cell, ct] = p1_folrd_df[ct].loc[cell]+1


sns.clustermap(p1_folrd_df)
plt.savefig('Alexandre_Paper_FOLR2/Patient1/Patient1_neighbour_FOLR2_heatmap.png', bbox_inches='tight')
plt.close('all')
sns.violinplot(data=p1_folrd_df)
plt.savefig('Alexandre_Paper_FOLR2/Patient1/Patient1_neighbour_FOLR2_violin.png', bbox_inches='tight')
plt.close('all')

#voisin +radius x
x= 100 #mean diameter tumour c = 11
name = 'Neighbour_'+str(x)
geom_p1 = sdata_p1.shapes['proseg_boundaries']['geometry']
tree_p1 = STRtree(geom_p1)
adata_patient1.obs[name] = 'nan'
adata_patient1.obs[name] = adata_patient1.obs[name].astype(object)
for cell in adata_patient1.obs['cell_id'].tolist():
    indices = tree_p1.query(geom_p1[cell],predicate='dwithin',distance=x).tolist()
    cell_touching = geom_p1.index[indices].tolist()
    cell_touching.remove(cell)
    adata_patient1.obs.at[cell, name] = cell_touching

adata_p1_folr2 = adata_patient1[adata_patient1.obs['Cell_Population']=='Mac.FOLR2'].copy()
p1_folr2_cell_type = pd.DataFrame(index = adata_p1_folr2.obs['cell_id'].tolist(),columns =adata_patient1.obs['Cell_Type'].cat.categories.tolist())
p1_folr2_cell_type = p1_folr2_cell_type.replace(np.NaN, 0)
p1_folr2_cell_subset = pd.DataFrame(index = adata_p1_folr2.obs['cell_id'].tolist(),columns =adata_patient1.obs['Cell_Subset'].unique().tolist())
p1_folr2_cell_subset = p1_folr2_cell_subset.replace(np.NaN, 0)

for cell in adata_p1_folr2.obs['cell_id'].tolist():
    voisins = adata_p1_folr2.obs[name][cell]
    for i in voisins :
        if i in adata_patient1.obs['cell_id'].tolist() :
            ct = adata_patient1.obs['Cell_Type'][i]
            p1_folr2_cell_type.loc[cell, ct] = p1_folr2_cell_type[ct].loc[cell]+1
            ss = adata_patient1.obs['Cell_Subset'][i]
            p1_folr2_cell_subset.loc[cell, ss] = p1_folr2_cell_subset[ss].loc[cell]+1

sns.clustermap(p1_folr2_cell_type)
plt.savefig('Alexandre_Paper_FOLR2/Patient1_neighbour_'+str(x)+'_FOLR2_heatmap.png')
plt.close('all')
sns.violinplot(p1_folr2_cell_type)
plt.savefig('Alexandre_Paper_FOLR2/Patient1_neighbour_'+str(x)+'_FOLR2_violin.png', bbox_inches='tight')
plt.close('all')
p1_folr2_cell_type.to_csv('Alexandre_Paper_FOLR2/Patient1/Patient1_neighbour_'+str(x)+'_CellType_FOLR2.csv')
p1_folr2_cell_subset.to_csv('Alexandre_Paper_FOLR2/Patient1/Patient1_neighbour_'+str(x)+'_CellSubset_FOLR2.csv')

#flat_list = [ x for xs in adata_p1_folr2.obs['Neighbour_10'] for x in xs]
#sc.pl.embedding(adata_patient1[adata_patient1.obs['cell_id'].isin(flat_list)], color="Cell_Type",basis='spatial')

p1_all_folr2_df = pd.DataFrame(index =adata_patient1.obs['Cell_Type'].cat.categories.tolist(),columns = ['0','10','20','30','40','50','60','70','80','90','100'])
p1_all_folr2_df = p1_all_folr2_df.replace(np.NaN, 0)
p1_new_folr2_df = pd.DataFrame(index =adata_patient1.obs['Cell_Type'].cat.categories.tolist(),columns = ['0','10','20','30','40','50','60','70','80','90','100'])
p1_new_folr2_df = p1_all_folr2_df.replace(np.NaN, 0)

for x in range(0,101,10):
    df = pd.read_csv('Alexandre_Paper_FOLR2/Patient1/Patient1_neighbour_'+str(x)+'_CellType_FOLR2.csv')
    sum = df.sum()
    p1_all_folr2_df[str(x)]=sum[1:].to_list()
    if int(x) >10:
        p1_new_folr2_df[str(x)] = p1_all_folr2_df[str(x)]-p1_all_folr2_df[str(x-10)]
    else:
        p1_new_folr2_df[str(x)] = p1_all_folr2_df[str(x)]

ax = p1_new_folr2_df.T.plot.bar(stacked=True,cmap='tab20')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Type", framealpha=0);
plt.savefig('Alexandre_Paper_FOLR2/Patient1/AllCells_at_X__foreach_FOLR2_p1_v2.png', bbox_inches='tight')
p1_new_folr2_df_norm = (p1_new_folr2_df/(p1_new_folr2_df.sum())).T
ax = p1_new_folr2_df_norm.plot.bar(stacked=True,cmap='tab20')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Type", framealpha=0);
plt.savefig('Alexandre_Paper_FOLR2/Patient1/AllCells_at_X_norm_foreach_FOLR2_p1_v2.png', bbox_inches='tight')
plt.close()
p1_new_folr2_df_norm.plot()
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Type", framealpha=0);
plt.savefig('Alexandre_Paper_FOLR2/Patient1/AllCells_at_X_normLine_foreach_FOLR2_p1_v2.png', bbox_inches='tight')

#idem but radius
sq.gr.spatial_neighbors_radius(adata_patient1, radius=10)
adj_matrix = adata_patient1.obsp["spatial_connectivities"]

p1_folr2_10_df = pd.DataFrame(index = adata_p1_folr2.obs['cell_id'].tolist(),columns =adata_patient1.obs['Cell_Type'].cat.categories.tolist())
p1_folr2_10_df = p1_folr2_10_df.replace(np.NaN, 0)
for cell in adata_p1_folr2.obs['cell_id'].tolist():
    target_indices = np.where(adata_patient1.obs['cell_id'] == cell)[0]
    neighbor_counts = np.array(adj_matrix[:, target_indices].sum(axis=1))
    voisins = adata_patient1[neighbor_counts > 0].obs['cell_id'].tolist()
    for i in voisins :
        if i in adata_patient1.obs['cell_id'].tolist() :
            ct = adata_patient1.obs['Cell_Type'][i]
            p1_folr2_10_df.loc[cell, ct] = p1_folr2_10_df[ct].loc[cell]+1

sns.clustermap(p1_folr2_10_df)
plt.show()

p = gpd.GeoSeries(geom_p1[voisins]) #if necessary to plot
p.plot() #.boundary.plot()
plt.show()
sc.pl.embedding(adata_patient1[adata_patient1.obs['cell_id'].isin(voisins)], color="Cell_Type",basis='spatial')
#adata_filtered = adata_patient1[neighbor_counts >0].copy()


#radius p1
x=10
sq.gr.spatial_neighbors_radius(adata_patient1, radius=x)
sq.tl.var_by_distance(adata_patient1,cluster_key="Cell_Population",groups="Mac.FOLR2",spatial_key="spatial")
distances_folr2 = adata_patient1.obsm["design_matrix"]["Mac.FOLR2_raw"]
adata_subset_p1 = adata_patient1[distances_folr2 <= x].copy()

sc.pl.embedding(adata_subset_p1, color="Cell_Type",basis='spatial')

df_dist= adata_patient1.obs[['Cell_Type','Cell_Subset','area']]
df_dist= df_dist.join(adata_patient1.obsm["design_matrix"]["Mac.FOLR2_raw"].to_frame())
df_dist['Distance'] = pd.cut(df_dist['Mac.FOLR2_raw'], [0,10,20,30,40,50,60,70,80,90,100])
dfn=df_dist[['Distance','area']].groupby(['Distance']).sum()
df_used = df_dist[['Cell_Type','Distance']]
df_dist_new = df_used.pivot_table(columns="Cell_Type", index="Distance", aggfunc=len, fill_value=0)
df_dist_density = df_dist_new.T.div(dfn.T.iloc[0], axis= 'columns').T


ax = df_dist_new.plot.bar(stacked=True,cmap='tab20')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Type", framealpha=0);
plt.savefig('Alexandre_Paper_FOLR2/AllCells_at_X_FOLR2_p1.png', bbox_inches='tight')
df_dist_new_norm = df_dist_new.apply(lambda x: x/sum(x), axis=1)
ax = df_dist_new_norm.plot.bar(stacked=True,cmap='tab20')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Type", framealpha=0);
plt.savefig('Alexandre_Paper_FOLR2/AllCells_at_X_norm_FOLR2_p1.png', bbox_inches='tight')
ax = df_dist_density.plot.bar(stacked=True,cmap='tab20')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Type", framealpha=0);
plt.savefig('Alexandre_Paper_FOLR2/AllCells_at_X_density_FOLR2_p1.png', bbox_inches='tight')


for i in adata_patient1.obs['Cell_Type'].unique():
    dfn_sub=df_dist[['Distance','area']].groupby(['Distance']).sum()
    df_used = df_dist[['Cell_Subset','Distance']]
    df_dist_new = df_used.pivot_table(columns="Cell_Subset", index="Distance", aggfunc=len, fill_value=0)
    filter_col = [col for col in df_dist_new if col.startswith(i)]
    df_fibro = df_dist_new[filter_col]
    df_fibro_density = df_fibro.T.div(dfn_sub.T.iloc[0], axis= 'columns').T
    ax = df_fibro.plot.bar(stacked=True,cmap='gist_rainbow')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    plt.style.use('ggplot')
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Population", framealpha=0);
    plt.savefig('Alexandre_Paper_FOLR2/Cell_subset_'+i+'_distance_to_FOLR2_p1.png', bbox_inches='tight')
    df_fibro_norm = df_dist_new[filter_col].apply(lambda x: x/sum(x), axis=1)
    ax = df_fibro_norm.plot.bar(stacked=True,cmap='gist_rainbow')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    plt.style.use('ggplot')
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Population", framealpha=0);
    plt.savefig('Alexandre_Paper_FOLR2/Cell_subset_'+i+'_norm_distance_to_FOLR2_p1.png', bbox_inches='tight')
    plt.close("all")
    ax = df_fibro_density.plot.bar(stacked=True,cmap='gist_rainbow')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    plt.style.use('ggplot')
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Population", framealpha=0);
    plt.savefig('Alexandre_Paper_FOLR2/Cell_subset_'+i+'_density_distance_to_FOLR2_p1.png', bbox_inches='tight')


#radius p2
x=10
sq.gr.spatial_neighbors_radius(adata_p2, radius=x)
sq.tl.var_by_distance(adata_p2,cluster_key="Cell_Subset",groups="Myeloid_#2",spatial_key="spatial")
distances_folr2 = adata_p2.obsm["design_matrix"]["Myeloid_#2_raw"]
adata_subset_p2 = adata_p2[distances_folr2 <= x].copy()

sc.pl.embedding(adata_subset_p2, color="Cell_Type",basis='spatial')

df_dist= adata_p2.obs[['Cell_Type','Cell_Subset']]
df_dist= df_dist.join(adata_p2.obsm["design_matrix"]["Myeloid_#2_raw"].to_frame())
df_dist['Distance'] = pd.cut(df_dist['Myeloid_#2_raw'], [0,10,20,30,40,50,60,70,80,90,100])
df_used = df_dist[['Cell_Type','Distance']]
df_dist_new = df_used.pivot_table(columns="Cell_Type", index="Distance", aggfunc=len, fill_value=0)
ax = df_dist_new.plot.bar(stacked=True,cmap='tab20')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Type", framealpha=0);
plt.savefig('Alexandre_Paper_FOLR2/AllCells_at_X_FOLR2_p2.png', bbox_inches='tight')
df_dist_new_norm = df_dist_new.apply(lambda x: x/sum(x), axis=1)
ax = df_dist_new_norm.plot.bar(stacked=True,cmap='tab20')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.style.use('ggplot')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Type", framealpha=0);
plt.savefig('Alexandre_Paper_FOLR2/AllCells_at_X_norm_FOLR2_p2.png', bbox_inches='tight')

for i in adata_patient1.obs['Cell_Type'].unique():
    df_used = df_dist[['Cell_Subset','Distance']]
    df_dist_new = df_used.pivot_table(columns="Cell_Subset", index="Distance", aggfunc=len, fill_value=0)
    filter_col = [col for col in df_dist_new if col.startswith(i)]
    df_fibro = df_dist_new[filter_col]
    ax = df_fibro.plot.bar(stacked=True,cmap='gist_rainbow')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    plt.style.use('ggplot')
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Population", framealpha=0);
    plt.savefig('Alexandre_Paper_FOLR2/Cell_subset_'+i+'_distance_to_FOLR2_p2.png', bbox_inches='tight')
    df_fibro_norm = df_dist_new[filter_col].apply(lambda x: x/sum(x), axis=1)
    ax = df_fibro_norm.plot.bar(stacked=True,cmap='gist_rainbow')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    plt.style.use('ggplot')
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="Cell Population", framealpha=0);
    plt.savefig('Alexandre_Paper_FOLR2/Cell_subset_'+i+'_norm_distance_to_FOLR2_p2.png', bbox_inches='tight')
    plt.close("all")


#Asign cluster to folr2 based on neigbors
#p1
df = pd.read_csv('Alexandre_Paper_FOLR2/Patient1/Patient1_neighbour_'+str(x)+'_FOLR2.csv', index_col=0)
kmeans = KMeans(n_clusters=6, random_state=42)
df['Cluster_ID'] = kmeans.fit_predict(df[adata_patient1.obs['Cell_Type'].cat.categories.tolist()])

network_pal = sns.color_palette('pastel', len(df['Cluster_ID'].unique()))
lut = dict(zip(df['Cluster_ID'].unique(), network_pal))
row_colors = df['Cluster_ID'].map(lut)
df_sort = df.sort_values(by=['Cluster_ID'])
sns.clustermap(df_sort.drop(columns=['Cluster_ID']), row_colors=row_colors, row_cluster=False, col_cluster=False, dendrogram_ratio=0.05, cbar_pos=None)
plt.savefig('Alexandre_Paper_FOLR2/Patient1/Patient1_neighbour_'+str(x)+'_FOLR2_heatmap_cluster.png')
plt.close('all') 

adata_p1_folr2 = adata_patient1[adata_patient1.obs['Cell_Population']=='Mac.FOLR2'].copy()
adata_p1_folr2.obs.loc[adata_p1_folr2.obs['cell_id'].isin(df.index), 'Cluster_Folr2'] = df['Cluster_ID']
adata_p1_folr2.obs['Cluster_Folr2']=+adata_p1_folr2.obs['Cluster_Folr2'].astype("category")
adata_p1_folr2.obs['Cluster_Folr2_str']=('#'+adata_p1_folr2.obs['Cluster_Folr2'].astype(str)).astype("category")

sc.tl.rank_genes_groups(adata_p1_folr2, 'Cluster_Folr2', method='wilcoxon', key_added = "wilcoxon")
markers = sc.get.rank_genes_groups_df(adata_p1_folr2,group=None,key='wilcoxon')
markers.to_csv('Alexandre_Paper_FOLR2/Patient1/Patient1_neighbour_'+str(x)+'_FOLR2_markers.csv')
sc.tl.rank_genes_groups(adata_p1_folr2, 'Cluster_Folr2_str', method='wilcoxon', key_added = "wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata_p1_folr2,n_genes=5, key="wilcoxon", groupby="Cluster_Folr2_str",min_logfoldchange=2,show=False)
plt.savefig('Alexandre_Paper_FOLR2/Patient1/Patient1_neighbour_'+str(x)+'_FOLR2_markersDot.png', bbox_inches='tight')
plt.close()

#p2
df_neig_p2 = pd.read_csv('Alexandre_Paper_FOLR2/Patient2/Patient2_neighbour_'+str(x)+'_CellType_FOLR2.csv', index_col=0)
kmeans = KMeans(n_clusters=6, random_state=42)
df_neig_p2['Cluster_ID'] = kmeans.fit_predict(df_neig_p2[adata_p2.obs['Cell_Type'].cat.categories.tolist()])

network_pal = sns.color_palette('pastel', len(df_neig_p2['Cluster_ID'].unique()))
lut = dict(zip(df_neig_p2['Cluster_ID'].unique(), network_pal))
row_colors = df_neig_p2['Cluster_ID'].map(lut)
df_neig_p2_sort = df_neig_p2.sort_values(by=['Cluster_ID'])
sns.clustermap(df_neig_p2_sort.drop(columns=['Cluster_ID']), row_colors=row_colors, row_cluster=False, col_cluster=False, dendrogram_ratio=0.05, cbar_pos=None)
plt.savefig('Alexandre_Paper_FOLR2/Patient2/Patient2_neighbour_'+str(x)+'_FOLR2_heatmap_cluster.png')
plt.close('all') 

adata_p2_folr2 = adata_p2[adata_p2.obs['Cell_Subset']=='Myeloid_#2'].copy()
adata_p2_folr2.obs.loc[adata_p2_folr2.obs['cell_id'].isin(df_neig_p2.index), 'Cluster_Folr2'] = df_neig_p2['Cluster_ID']
adata_p2_folr2.obs['Cluster_Folr2']=adata_p2_folr2.obs['Cluster_Folr2'].astype("category")
adata_p2_folr2.obs['Cluster_Folr2_str']=('#'+adata_p2_folr2.obs['Cluster_Folr2'].astype(str)).astype("category")

sc.tl.rank_genes_groups(adata_p2_folr2, 'Cluster_Folr2', method='wilcoxon', key_added = "wilcoxon")
markers = sc.get.rank_genes_groups_df(adata_p2_folr2,group=None,key='wilcoxon')
#markers.to_csv('Alexandre_Paper_FOLR2/Patient2/Patient2_neighbour_'+str(x)+'_FOLR2_markers.csv')
sc.tl.rank_genes_groups(adata_p2_folr2, 'Cluster_Folr2_str', method='wilcoxon', key_added = "wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata_p2_folr2,n_genes=5, key="wilcoxon", groupby="Cluster_Folr2_str",min_logfoldchange=1.5,show=False)
plt.savefig('Alexandre_Paper_FOLR2/Patient2/Patient2_neighbour_'+str(x)+'_FOLR2_markersDot.png', bbox_inches='tight')
plt.close()

'''
Next to FOLR2 changes
'''
c_type = ['B cells', 'Fibroblasts', 'Plasma cells', 'T cells']
adata_p2_folr2_c_list = adata_p2[adata_p2.obs['Cell_Subset']=='Myeloid_#2'].obs['cell_id'].tolist()
for ct in c_type:
    adata = adata_p2[adata_p2.obs['Cell_Type']==ct].copy()
    adata.obs['Voisin_FOLR2']= '0'
    for cell in adata.obs['cell_id'].tolist():
        val_forl2 = len(set(adata.obs['Neighbour_shapely'][cell]).intersection(adata_p2_folr2_c_list))
        adata.obs.at[cell, 'Voisin_FOLR2'] = str(val_forl2)
    adata.obs['Voisin_FOLR2_YN'] = ['Y' if int(i) >0 else 'N' for i in adata.obs['Voisin_FOLR2']]
    sc.tl.rank_genes_groups(adata, groupby='Voisin_FOLR2_YN', key_added='group_DE_results')
    marker_folr2_neigh = sc.get.rank_genes_groups_df(adata,group=None,key='group_DE_results')
    marker_folr2_neigh.to_csv('Alexandre_Paper_FOLR2/Patient2/Patient2_'+ct+'_markers.csv')
    sc.pl.rank_genes_groups_dotplot(adata, groupby='Voisin_FOLR2_YN',key='group_DE_results',min_logfoldchange=1,show=False)
    plt.savefig('Alexandre_Paper_FOLR2/Patient2/Patient2_'+ct+'_markersDot.png', bbox_inches='tight')
    plt.close()

c_type = ['B cells and plasma cells', 'Fibroblasts','T cells']
adata_p1_folr2_c_list = adata_patient1[adata_patient1.obs['Cell_Population']=='Mac.FOLR2'].obs['cell_id'].tolist()
for ct in c_type:
    adata = adata_patient1[adata_patient1.obs['Cell_Type']==ct].copy()
    adata.obs['Voisin_FOLR2']= '0'
    for cell in adata.obs['cell_id'].tolist():
        val_forl2 = len(set(adata.obs['Neighbour_shapely'][cell]).intersection(adata_p2_folr2_c_list))
        adata.obs.at[cell, 'Voisin_FOLR2'] = str(val_forl2)
    adata.obs['Voisin_FOLR2_YN'] = ['Y' if int(i) >0 else 'N' for i in adata.obs['Voisin_FOLR2']]
    sc.tl.rank_genes_groups(adata, groupby='Voisin_FOLR2_YN', key_added='group_DE_results')
    marker_folr2_neigh = sc.get.rank_genes_groups_df(adata,group=None,key='group_DE_results')
    marker_folr2_neigh.to_csv('Alexandre_Paper_FOLR2/Patient1/Patient1_'+ct+'_markers.csv')
    sc.pl.rank_genes_groups_dotplot(adata, groupby='Voisin_FOLR2_YN',key='group_DE_results',min_logfoldchange=1,show=False)
    plt.savefig('Alexandre_Paper_FOLR2/Patient1/Patient1_'+ct+'_markersDot.png', bbox_inches='tight')
    plt.close()