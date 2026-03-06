library(Seurat,future)
library(dplyr)
suppressMessages(library(tidyverse))
#plan(workers = 2)
options(future.globals.maxSize = 8000 * 1024^2)

myelo_1 <- readRDS( "Integrated_myelo_1.rds")
#myelo_1 <- readRDS( "Integrated_myelo_1.rds")
myelo_2[["CellName"]] <- colnames(myelo_2)

cell_tokeep = suppressMessages(as.data.frame(read_csv("Batch2_13datasets_total_DC.csv", skip = 5)))
cell_tokeep = colnames(cell_tokeep)
cell_tokeep = cell_tokeep[-c(1,2)]

myelo_2_s <- subset(myelo_2, subset = CellName %in% cell_tokeep)
myelo_2_s <- RunPCA(myelo_2_s, npcs = 100)
myelo_2_s <- RunICA(myelo_2_s, nics = 100)
myelo_2_s <- FindNeighbors(myelo_2_s, reduction = "harmony", dims = 1:30)
myelo_2_s <- RunUMAP(myelo_2_s, dims = 1:30, reduction = "harmony")
myelo_2_s <- FindClusters(myelo_2_s, resolution = 0.2)
saveRDS( myelo_2_s, "Integrated_myelo_2_total_DC.rds")


DefaultAssay(myelo_2) <- 'RNA'
DefaultAssay(myelo_1) <- 'RNA'

myelo_2 <- JoinLayers(myelo_2)
myelo_2 <- SCTransform(myelo_2)

myelo_1 <- JoinLayers(myelo_1)
myelo_1 <- SCTransform(myelo_1)

myelo_DC <- merge(myelo_1, myelo_2, project = "integrated_total")
myelo_DC <- SCTransform(myelo_DC, vars.to.regress = "percent.mt", verbose = FALSE)
gc()
myelo_DC <- RunPCA(myelo_DC, npcs = 100)
myelo_DC <- RunICA(myelo_DC, nics = 100)
myelo_DC <- FindNeighbors(myelo_DC, dims = 1:10, reduction = "pca")
myelo_DC <- RunUMAP(myelo_DC, dims = 1:10, reduction = "pca", reduction.name = "umap.unintegrated")

myelo_DC <- IntegrateLayers(object = myelo_DC, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",assay = "SCT", verbose = FALSE)
myelo_DC[["RNA"]] <- JoinLayers(myelo_DC[["RNA"]])
myelo_DC <- FindNeighbors(myelo_DC, reduction = "harmony", dims = 1:30)
myelo_DC <- RunUMAP(myelo_DC, dims = 1:30, reduction = "harmony")

saveRDS(myelo_DC, "Integrated_myelo_total_DC.rds")

G_Header <- read.table("../pymt_scrape/Header_SeqGeq.txt",sep="\t",header=F,row.names=1,check.names = F)
SG_Header <- unname(SG_Header)

counts_matrix <- as.data.frame(myelo_DC[["SCT"]]$data) %>%
  filter(row.names(.) %in% genes_signature)

umap = t(myelo_DC[["umap"]]@cell.embeddings)
sample=myelo_DC$orig.ident
dataset=myelo_DC$orig.ident
pc=t(myelo_DC@reductions$pca@cell.embeddings)
ica=t(myelo_DC@reductions$ica@cell.embeddings)

total = rbind(counts_matrix,umap,sample,dataset,pc,ica)

write.table(SG_Header, file = "Integrated_myelo_total_DC.txt",sep="\t")
suppressWarnings(write.table(total, file = "Integrated_myelo_total_DC.txt",sep="\t",append=T))

