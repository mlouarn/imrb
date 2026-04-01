# int_gse gsedir
library(Seurat)
library(stringr)
library(future)
options(future.globals.maxSize = 8000 * 1024^2)

args = commandArgs(trailingOnly = TRUE)
gse_dir = args[1]
setwd(gse_dir)

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

add_origin = function(seurat_obj, GSM){
  seurat_obj$orig.ident <- GSM
  return(seurat_obj)
}

rds_files = list.files('.', pattern = 'GSM.*\\.rds')
gsm_list = str_extract(rds_files, 'GSM[0-9]{5,8}')
list_seurat_obj = lapply(rds_files, readRDS)
list_seurat_obj = mapply(add_origin, list_seurat_obj, gsm_list)
seurat_obj_int = seurat_integration(list_seurat_obj)
saveRDS(seurat_obj_int, paste0(gse_dir, '.rds'))