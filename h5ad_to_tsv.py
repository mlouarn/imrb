import scanpy as sc
import pandas as pd


adata = sc.read_h5ad("./pbmc3k.h5ad")
adata.X = adata.layers["data"]
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

output_basename = "pbmc3k"

with open(f"{output_basename}_matrix.tsv", "w") as file:
        mat = adata.X.T # gene x cell
        for index in range(len(adata.obs_names)):
            row = mat[index,:].toarray()[0]
            row_str = [str(x) for x in row]
            file.write("\t".join(row_str)+"\n")

pd.DataFrame(adata.var_names).to_csv(f"{output_basename}_features.tsv", header=False, index = False)

pd.DataFrame(adata.obs_names).to_csv(f"{output_basename}_barcodes.tsv", header=False, index = False)

df_umap = pd.DataFrame(adata.obsm["X_umap"], columns = ["umap_1", "umap_2"])
df_umap.to_csv(f"{output_basename}_reduction.tsv", sep = "\t")

adata.obs.to_csv(f"{output_basename}_metadata.tsv", sep = "\t")