# usage : int_harmony.py anndata1.h5ad anndata2.h5ad [...] outputfile
import anndata as ad
import scanpy as sc
import sys
#import scanpy.external as sce
# there is a known problem with scanpy external harmony
import harmonypy as hm
import re
import pickle

def integrate_harmony(adata, batch_effect_key):
    sc.pp.pca(adata, n_comps=50)
    harmony_output = hm.run_harmony(adata.obsm['X_pca'], 
                                    adata.obs, batch_effect_key)
    adata.obsm['X_pca'] = harmony_output.Z_corr
    sc.pp.neighbors(adata, n_neighbors=30, 
            n_pcs=10, use_rep='X_pca')
    sc.tl.umap(adata, min_dist=0.05)
    sc.tl.leiden(adata)
    return(adata)

if __name__=="__main__":
    ad.settings.allow_write_nullable_strings=True # mandatory to write to h5ad

    h5ad_filenames = sys.argv[1:-1]
    output = sys.argv[-1]
    print(h5ad_filenames)

    # process
    list_adata = []
    for h5ad in h5ad_filenames:
        print(h5ad)
        adata = sc.read_h5ad(h5ad)
        print(adata)
        adata.X = adata.layers["data"]
        list_adata.append(adata)
    
    # concatenate and integrate
    GSE_id = [re.findall(r"GS[ME][0-9]{5,9}", h5ad)[0] for h5ad in h5ad_filenames]
    adata = ad.concat(list_adata, label="GSE", 
                    keys=GSE_id, index_unique="-") 
    del list_adata
    adata = integrate_harmony(adata, "GSE")

    sc.pl.umap(adata, color="GSE", save="_v3.png")
    # save
    with open(output+'.pkl', 'wb') as file:
        pickle.dump(adata, file)
    adata.write_h5ad(output)
