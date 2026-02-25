library(anndataR)
library(SeuratObject)
library(stringr)

list_rds = list.files(pattern = "\\.rds$")
for(filename_rds in list_rds){
  if(file.exists(str_replace(filename_rds, "rds$", "h5ad"))){
    next()
  }
  print(filename_rds)
  seurat_obj = readRDS(filename_rds)
  adata_obj = as_AnnData(seurat_obj)
  write_h5ad(adata_obj, str_replace(filename_rds,
                                    pattern = "\\.rds$",
                                    replacement = ".h5ad"))
  print("done")
}

