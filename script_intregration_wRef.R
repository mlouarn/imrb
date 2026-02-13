#!/usr/bin/env Rscript

library(colorRamp2)
library(dplyr)
library(ggplot2)
library(magrittr)
library(patchwork)
library(RColorBrewer)
library(Seurat)
library(UCell)


#######Fuctions
seurat_integration<- function(list_seurat){
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
signatures_mouse = read.csv('/home/marinelouarn/Documents/Alexandre_mouse/20260409_Mouse_cellPopulation_Signature_Genes_Conensus.csv')
signatures_mouse = signatures_mouse[signatures_mouse$Keep_for_initial_screening=='y',]
add_signatures_mouse <- function(obj_seurat, signature_file, level_signature){
  for(tissue in as.list(unique(signature_file[level_signature]))[[1]]){
    if (tissue!=""){
      genes_mouse = signature_file[signature_file[level_signature]==tissue,]$Signature_gene_mouse.style_symbol
      obj_seurat <- suppressWarnings(AddModuleScore_UCell(obj_seurat,
                                         assay= 'RNA',
                                         features = list(genes_mouse),
                                         name = tissue))
      obj_seurat[[as.character(tissue)]] <- obj_seurat[[paste0('signature_1',tissue)]]
    }
  }
  gc()
  return(obj_seurat)
}

####Main
folder_GSE = commandArgs(trailingOnly=TRUE)
reference = readRDS('../zenodo-15826764/zenodo-15826764.rds')
catch ='nice'

rds_files = list.files(path = folder_GSE, pattern = '*.rds$')
list_seurat = lapply(rds_files, readRDS)
list_seurat[[length(list_seurat)+1]] <- reference

integrated = seurat_integration(list_seurat)

tryCatch({
  integrated <- add_signatures_mouse(integrated,signatures_mouse,'LEVEL_2')
}, error = function(e) {
  print("An error occurred:", e$message)
  catch = 'error'
})


png("UMAP.png",width = 600, height = 500, unit = "px")
  DimPlot(integrated, reduction = "umap",label = T)
dev.off()

png("UMAP_ident.png",width = 600, height = 500, unit = "px")
  DimPlot(integrated, reduction = "umap", group.by = "orig.ident")
dev.off()

if (catch !='error'){
  png("UMAP_signaturesLvl2.png",width = 600, height = 500, unit = "px")
    FeaturePlot(integrated,features=unique(signatures_mouse$LEVEL_2))& 
      scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  dev.off()
    
  png("Vln_signaturesLvl2.png",width = 600, height = 500, unit = "px")
    VlnPlot(integrated,features=unique(signatures_mouse$LEVEL_2), group.by='seurat_clusters')
  dev.off()
}

saveRDS(integrated, file = "????.rds")

SG_Header <- read.table("Header_SeqGeq.txt",sep="\t",header=F,row.names=1,check.names = F)
SG_Header <- unname(SG_Header)
counts_matrix <- as.data.frame(integrated[["integrated"]]$data)
write.table(SG_Header, file = "test.txt",sep="\t")
suppressWarnings(write.table(counts_matrix, file = "test.txt",sep="\t",append=T))
