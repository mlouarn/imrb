# input one argument : rds filename
# writes the same file as h5ad

library(anndataR)
library(SeuratObject)
library(stringr)

filename_rds = commandArgs(trailingOnly = T)[1]
print(filename_rds)

seurat_obj = readRDS(filename_rds)
adata_obj = as_AnnData(seurat_obj, assay_name="RNA")
adata_obj$X = adata_obj$layers['data'][[1]]
write_h5ad(adata_obj, str_replace(filename_rds, 
                                  pattern = "\\.rds$",
                                  replacement = ".h5ad"))
