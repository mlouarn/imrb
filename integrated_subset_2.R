library(Seurat,future)
library(dplyr)
suppressMessages(library(tidyverse))
#plan(workers = 2)
options(future.globals.maxSize = 8000 * 1024^2)

myelo_2 = readRDS("Integrated_myelo_2_DC2.rds")
myelo_1 =readRDS("Integrated_myelo_1_DC2.rds")
DefaultAssay(myelo_2) <- 'RNA'
DefaultAssay(myelo_1) <- 'RNA'

myelo_2 <- JoinLayers(myelo_2)
myelo_2 <- SCTransform(myelo_2)

myelo_1 <- JoinLayers(myelo_1)
myelo_1 <- SCTransform(myelo_1)

myelo_DC <- merge(myelo_1, myelo_2, project = "integrated_total")

frequencies <- table(myelo_DC$orig.ident)
sample_to_keep = names(frequencies)[frequencies > 30] 

myelo_DC <- subset(myelo_DC, subset = orig.ident %in% sample_to_keep)

DefaultAssay(myelo_DC) <- 'RNA'
myelo_DC[["RNA"]] <- JoinLayers(myelo_DC[["RNA"]])
myelo_DC[["RNA"]] <- split(myelo_DC[["RNA"]], f = myelo_DC$orig.ident)
myelo_DC <- NormalizeData(myelo_DC)
myelo_DC <- FindVariableFeatures(myelo_DC)
myelo_DC <- ScaleData(myelo_DC)
myelo_DC <- RunPCA(myelo_DC)

myelo_DC <- FindNeighbors(myelo_DC, dims = 1:30, reduction = "pca")
myelo_DC <- FindClusters(myelo_DC, resolution = 2, cluster.name = "unintegrated_clusters")

myelo_int <- IntegrateLayers(object = myelo_DC, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",k.weight=30, verbose = FALSE)
myelo_int <- RunICA(myelo_int,nics=100)
myelo_int <- RunPCA(myelo_int,npcs=100)
myelo_int <- FindNeighbors(myelo_int, reduction = "integrated.cca", dims = 1:30)
myelo_int <- FindClusters(myelo_int, resolution = 0.2)
myelo_int <- RunUMAP(myelo_int, dims = 1:30, reduction = "integrated.cca")
myelo_int[["RNA"]] <- JoinLayers(myelo_int[["RNA"]])


png("DCtotal_umap_cca.png",width = 1600, height = 500, unit = "px")
	DimPlot(myelo_int, reduction = "umap", group.by = c("orig.ident"))
dev.off()


saveRDS(myelo_int, "Integrated_cca_myelo_total_DC.rds")

SG_Header <- read.table("../pymt_scrape/Header_SeqGeq.txt",sep="\t",header=F,row.names=1,check.names = F)
SG_Header <- unname(SG_Header)

counts_matrix <- as.data.frame(myelo_int[["RNA"]]$data) %>%
  filter(row.names(.) %in% genes_signature)

umap = t(myelo_int[["umap"]]@cell.embeddings)
sample=myelo_int$orig.ident
dataset=myelo_int$orig.ident
samp_data=myelo_int$orig.ident
pc=t(myelo_int@reductions$pca@cell.embeddings)
ica=t(myelo_int@reductions$ica@cell.embeddings)
cluster = myelo_int$RNA_snn_res.0.2

total = rbind(counts_matrix,umap,sample,samp_data,dataset,pc,ica,cluster)
rownames(total)[rownames(total) == "205"] <- "Sample"
rownames(total)[rownames(total) == "206"] <- "Sample_in_Dataset"
rownames(total)[rownames(total) == "207"] <- "Dataset"
rownames(total)[rownames(total) == "408"] <- "Clusters.0.2"


write.table(SG_Header, file = "Integrated_cca_myelo_total_DC.txt",sep="\t")
suppressWarnings(write.table(total, file = "Integrated_cca_myelo_total_DC.txt",sep="\t",append=T))

