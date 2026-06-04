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
        #f.write("leiden_res1.5_IC2PC10" + "\t" + "\t".join(leiden_myelo_str) + "\n")
        f.write("nCount" + "\t" + "\t".join(n_counts_str) + "\n")
        f.write("CellTypist_v1_nb" + "\t" + "\t".join(adata.obs['CellTypist_v1_nb']) + "\n")
        f.write("CellTypist_v2_nb" + "\t" + "\t".join(adata.obs['CellTypist_v2_nb']) + "\n")
        f.write("Tangram_nb" + "\t" + "\t".join(adata.obs['Tangram_nb']) + "\n")
        f.write("SingleR_nb" + "\t" + "\t".join(adata.obs['SingleR_nb']) + "\n")

def to_seqgeq(adata, header_sg=str,file=str,genelist=list):
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
        spatial_x = adata.obsm['spatial'][:,0]
        spatial_x_str = [str(x) for x in spatial_x]
        spatial_y = adata.obsm['spatial'][:,1]
        spatial_y_str = [str(y) for y in spatial_y]
        n_counts = adata.obs['n_counts']
        n_counts_str = [str(y) for y in n_counts]
        f.write("umap_PC1" + "\t" + "\t".join(umap_x_str)+"\n")
        f.write("umap_PC2" + "\t" + "\t".join(umap_y_str)+"\n")
        f.write("umap_ICPC1" + "\t" + "\t".join(umap_icpc_x_str)+"\n")
        f.write("umap_ICPC2" + "\t" + "\t".join(umap_icpc_y_str)+"\n")
        f.write("Spatial_X" + "\t" + "\t".join(spatial_x_str)+"\n")
        f.write("Spatial_Y" + "\t" + "\t".join(spatial_y_str)+"\n")
        f.write("leiden_clusters_PC" + "\t" + "\t".join(adata.obs['leiden']) + "\n")
        f.write("leiden_clusters_ICPC" + "\t" + "\t".join(adata.obs['leiden_ICPC']) + "\n")
        f.write("nCount" + "\t" + "\t".join(n_counts_str) + "\n")
        f.write("novae_domain_nb" + "\t" + "\t".join(adata.obs['novae_domain_nb']) + "\n")
