# integrate_rds rds1 rds2 rds3 ... rds_output
library(Seurat)
library(stringr)

arguments = commandArgs(trailingOnly = TRUE)
# rds_files = c("GSE132582/GSE132582_myeloid.rds", "GSE150176/GSE150176_myeloid.rds", "GSE171602/GSE171602_myeloid.rds", 
#             "GSE197329/GSE197329_myeloid.rds", "GSE220959/GSE220959_myeloid.rds", "GSE246458/GSE246458_myeloid.rds", "GSE249540/GSE249540_myeloid.rds",
#             "GSE275785/GSE275785_myeloid.rds", "GSE276238/GSE276238_myeloid.rds", 
#             "GSE283198/GSE283198_myeloid.rds", "GSE289401/GSE289401_myeloid.rds", 
#             "GSE296677/GSE296677_myeloid.rds", "GSE307209/GSE307209_myeloid.rds")
rds_out = arguments[length(arguments)]
rds_files = arguments[1:length(arguments)-1]

seurat_integration <- function(list_seurat){
  int.anchors <- FindIntegrationAnchors(object.list = list_seurat, dims = 1:10)
  to_integrate <- Reduce(intersect, lapply(int.anchors@object.list, rownames))
  int.combined <- IntegrateData(anchorset = int.anchors,features.to.integrate=to_integrate, dims = 1:10)
  DefaultAssay(int.combined) <- "RNA"
  int.combined <- JoinLayers(int.combined, assay = "RNA")
  DefaultAssay(int.combined) <- "integrated"
  
  int.combined <- ScaleData(int.combined, verbose = FALSE)
  int.combined <- RunPCA(int.combined, npcs = 50, verbose = FALSE)
  int.combined <- RunUMAP(int.combined, reduction = "pca", dims = 1:10)
  int.combined <- FindNeighbors(int.combined, reduction = "pca", dims = 1:10)
  int.combined <- FindClusters(int.combined, resolution = 1)
  return(int.combined)
}

list_seurat_obj = lapply(rds_files, readRDS)
names(list_seurat_obj) = rds_files
print(list_seurat_obj)
print(lapply(list_seurat_obj, function(seu){
  return(rownames(seu)[1:30])
}))
seurat_obj_int = seurat_integration(list_seurat_obj)
saveRDS(seurat_obj_int, rds_out)