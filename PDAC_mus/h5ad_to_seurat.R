# input one argument : h5ad filename
# writes the same file as an rds

library(Seurat)
library(stringr)
library(anndataR)

seurat_analyse <- function(obj_seurat){
  obj_seurat[["percent.mt"]] <- PercentageFeatureSet(obj_seurat, pattern = "^mt-")
  obj_seurat <- subset(obj_seurat, subset = nFeature_RNA > 200& percent.mt <= 20)
  obj_seurat <- NormalizeData(obj_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  obj_seurat <- FindVariableFeatures(obj_seurat, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(obj_seurat)
  obj_seurat <- ScaleData(obj_seurat, features = all.genes)
  
  obj_seurat <- RunPCA(obj_seurat, features = VariableFeatures(object = obj_seurat))
  
  obj_seurat <- FindNeighbors(obj_seurat, dims = 1:10)
  obj_seurat <- FindClusters(obj_seurat, resolution = 1)
  obj_seurat <- RunUMAP(obj_seurat, dims = 1:10)
  return(obj_seurat)
}

filename_h5ad = commandArgs(trailingOnly = T)[1]
print(filename_h5ad)

adata = read_h5ad(filename_h5ad)
seurat_obj = adata$as_Seurat()
rownames(seurat_obj)
seurat_obj = seurat_analyse(seurat_obj) 
saveRDS(seurat_obj, str_replace(filename_h5ad, 
                                  pattern = "\\.h5ad$",
                                  replacement = ".rds"))
