import anndata as ad
import scanpy as sc
import sys

header_sg="/home/marine-louarn/Documents/Header_SeqGeq.txt"
signature = pd.read_csv("/home/marine-louarn/Documents/test/20260204_Signatures_Cell_populations_HUMAN.csv")
genelist= signature['gene'].tolist()
file = "/home/marine-louarn/Documents/Xenium_Calderaro/full_seqgeq.txt"

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

if __name__=="__main__":
    h5ad = sys.argv[1]
    outfilename = sys.argv[2]
    genelist_file = sys.argv[3]
    header_sg = sys.argv[4]

    adata = sc.read_h5ad(h5ad)
    adata.X = adata.layers['data']

    if genelist_file=='all':
        genelist=adata.var_names
    else:
        with open(genelist_file, 'r') as file:
            genelist = [line.strip() for line in file]
        genelist = genelist[1:]


to_seqgeq(adata, file, genelist, header_sg)

#test to look at matrix
tmp = sdata.tables['table']
tmp.layers["counts"] = tmp.X.copy()

test.X = tmp.layers["counts"]
test = tmp.to_df().T
mnp_pseudobulk_mat = mnp_pseudobulk_mat[mnp_pseudobulk_mat.index.isin(genes_to_keep)]