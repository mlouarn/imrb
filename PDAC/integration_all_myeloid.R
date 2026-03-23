suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

setwd('/home/marine-louarn/Documents/PDAC/')

rds_files = list.files(path = 'Datasets/', pattern = '*_myeloid.rds$')
list_seurat = lapply(paste0('Datasets/',rds_files), readRDS)

myeloid <- merge(list_seurat[[1]], list_seurat[-1], project = "PDAC_myeloid")
frequencies <- table(myeloid$orig.ident)
sample_to_keep = names(frequencies)[frequencies > 30] 
myeloid <- subset(myeloid, subset = orig.ident %in% sample_to_keep)

DefaultAssay(myeloid) <- 'RNA'
myeloid[["RNA"]] <- JoinLayers(myeloid[["RNA"]])
myeloid[["RNA"]] <- split(myeloid[["RNA"]], f = myeloid$orig.ident)
myeloid <- NormalizeData(myeloid)
myeloid <- FindVariableFeatures(myeloid)
myeloid <- ScaleData(myeloid)
myeloid <- RunPCA(myeloid)
myeloid <- FindNeighbors(myeloid, dims = 1:30, reduction = "pca")
myeloid <- FindClusters(myeloid, resolution = 2, cluster.name = "unintegrated_clusters")

PDAC_myelo <- IntegrateLayers(object = myeloid, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",k.weight=30, verbose = FALSE)
PDAC_myelo[["RNA"]] <- JoinLayers(PDAC_myelo[["RNA"]])
PDAC_myelo <- FindNeighbors(PDAC_myelo, reduction = "integrated.cca", dims = 1:50)
for (i in seq(0.2,1.2,0.2)){
  PDAC_myelo <- FindClusters(PDAC_myelo, resolution = i)
}
PDAC_myelo <- RunUMAP(PDAC_myelo, dims = 1:50, reduction = "integrated.cca")

saveRDS(PDAC_myelo,'humanPDAC_integrated.rds')
seqgeq_file(PDAC_myelo,signature_human,"humanPDAC_integrated.txt")

####DC
PDAC_myelo <- readRDS('humanPDAC_integrated.rds')
cell_tokeep = suppressMessages(as.data.frame(read_csv("20260323_human_PDAC_integrated_Total_DC_Extracted_CellIDs.csv", skip = 5)))
cell_tokeep = colnames(cell_tokeep)
cell_tokeep = cell_tokeep[-c(1)]
PDAC_myelo[["CellName"]] <- colnames(PDAC_myelo)
total_DC <- subset(PDAC_myelo, subset = CellName %in% cell_tokeep)
total_DC <- subset(total_DC, subset = !(orig.ident %in% c('Zenodo_N1','Zenodo_N2','Zenodo_N3','Zenodo_N4','Zenodo_N5','Zenodo_N6','Zenodo_N7','Zenodo_N8','Zenodo_N9','Zenodo_N10','Zenodo_N11','Zenodo_T1','Zenodo_T2','Zenodo_T3','Zenodo_T4','Zenodo_T5','Zenodo_T6','Zenodo_T7','Zenodo_T8','Zenodo_T9','Zenodo_T10','Zenodo_T11','Zenodo_T12','Zenodo_T13','Zenodo_T14','Zenodo_T15','Zenodo_T16','Zenodo_T17','Zenodo_T18','Zenodo_T19','Zenodo_T20','Zenodo_T21','Zenodo_T22','Zenodo_T23','Zenodo_T24')))
frequencies <- table(total_DC$orig.ident)
sample_to_keep = names(frequencies)[frequencies > 30] 
total_DC <- subset(total_DC, subset = orig.ident %in% sample_to_keep)


DefaultAssay(total_DC) <- 'RNA'
total_DC[["RNA"]] <- JoinLayers(total_DC[["RNA"]])
total_DC[["RNA"]] <- split(total_DC[["RNA"]], f = total_DC$orig.ident)
total_DC <- NormalizeData(total_DC)
total_DC <- FindVariableFeatures(total_DC)
total_DC <- ScaleData(total_DC)
total_DC <- RunPCA(total_DC)
total_DC <- FindNeighbors(total_DC, dims = 1:30, reduction = "pca")
total_DC <- FindClusters(total_DC, resolution = 2, cluster.name = "unintegrated_clusters")

PDAC_DC <- IntegrateLayers(object = total_DC, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",k.weight=30, verbose = FALSE)
PDAC_DC[["RNA"]] <- JoinLayers(PDAC_DC[["RNA"]])
PDAC_DC <- FindNeighbors(PDAC_DC, reduction = "integrated.cca", dims = 1:50)
for (i in seq(0.2,1.2,0.2)){
  PDAC_DC <- FindClusters(PDAC_DC, resolution = i)
}
PDAC_DC <- RunUMAP(PDAC_DC, dims = 1:50, reduction = "integrated.cca")
DimPlot(PDAC_DC, reduction = "umap", group.by = c("RNA_snn_res.1"),label=TRUE,cols=DiscretePalette(28,palette='polychrome'))

saveRDS(PDAC_DC,'humanPDAC_integrated_total_DC.rds')
seqgeq_file(PDAC_DC,signature_human,"humanPDAC_integrated_total_DC.txt")

counts_matrix_full <- as.data.frame(PDAC_DC[["RNA"]]$data) 
umap = t(PDAC_DC[["umap"]]@cell.embeddings)
dataset=PDAC_DC$orig.ident
patient=PDAC_DC$orig.ident
cluster1 = PDAC_DC$RNA_snn_res.1
tissu=PDAC_DC$orig.ident

total = rbind(counts_matrix_full,umap,dataset,patient,tissu,cluster1)
rownames(total)[rownames(total) == "44271"] <- "Dataset"
rownames(total)[rownames(total) == "44272"] <- "Patient"
rownames(total)[rownames(total) == "44273"] <- "Tissue"
rownames(total)[rownames(total) == "44274"] <- "Cluster.res1"
write.table(SG_Header, file = "humanPDAC_integrated_total_DC_full.txt",sep="\t")
suppressWarnings(write.table(total, file = "humanPDAC_integrated_total_DC_full.txt",sep="\t",append=T))
