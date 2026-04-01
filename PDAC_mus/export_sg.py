import anndata as ad
import scanpy as sc
import sys

def write_seqgeq(adata:ad.AnnData, file:str, genelist:list, header_sg:str):
    with open(file, "w") as f:
        # write header
        with open(header_sg) as header:
            f.write(header.read())
        # write barcodes
        f.write("\t".join(adata.obs_names)+"\n")
        # keep only genes in signatures
        genes_absent = [gene for gene in genelist if gene not in adata.var_names]
        genes_present = [gene for gene in genelist if gene in adata.var_names]
        genes_index = [adata.var_names.get_loc(gene) for gene in genes_present]
        print(f'genes absent from adata {genes_absent}')
        # print(f'genes present from adata {genes_present}')

        # write normalized counts
        mat = adata.layers["data"].T # gene x cell
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
        f.write("umap_1" + "\t" + "\t".join(umap_x_str)+"\n")
        f.write("umap_2" + "\t" + "\t".join(umap_y_str)+"\n")

        # write clusters
        f.write("leiden_clusters" + "\t" + "\t".join(adata.obs['leiden']) + "\n")

        # GSE GSM=orig.ident en numérique et donner l'ordre
        # make dico for correspondance
        # dicogsm = {}
        # n=0
        # for GSM in sorted(set(adata.obs['orig.ident'])):
        #     dicogsm[GSM] = str(n)
        #     n+=1
        # adata.obs['GSM_num'] = [dicogsm[GSM] for GSM in adata.obs['orig.ident']]

        # dicogse = {}
        # n=0
        # for GSE in sorted(set(adata.obs['GSE'])):
        #     dicogse[GSE] = str(n)
        #     n+=1
        # adata.obs['GSE_num'] = [dicogse[GSE] for GSE in adata.obs['GSE']]

        # # print the correspondance tables
        # print("GSE table: ")
        # print(dicogse)
        # print("GSM table: ")
        # print(dicogsm)
        # with open(file+'_corresp.txt','w') as corresp:
        #     for gsm in dicogsm.keys():
        #         corresp.write(f'{gsm}\t{dicogsm[gsm]}\n')
        #     for gse in dicogse.keys():
        #         corresp.write(f'{gse}\t{dicogse[gse]}\n')

        # write to seqgeq the dataset (GSE), sample (GSM) and treatment
        f.write("dataset_pymt" + "\t" + "\t".join(adata.obs['dataset_pymt'])+"\n")
        f.write("dataset_pdac" + "\t" + "\t".join(adata.obs['dataset'])+"\n")
        f.write("sample" + "\t" + "\t".join(adata.obs['sample'])+"\n")
        f.write("treatment" + "\t" + "\t".join(adata.obs['treatment'])+"\n")

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

    # sc.pp.pca(adata)
    # sc.pp.neighbors(adata)
    # sc.tl.umap(adata)
    # sc.tl.leiden(adata)
    # print(adata)

    write_seqgeq(adata, outfilename, genelist, header_sg)