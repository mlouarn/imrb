source('/home/marine-louarn/imrb/packages.R')
source('/home/marine-louarn/imrb/functions.R')

setwd('/home/marine-louarn/Documents/PDAC/')

ref_mreg <- readRDS('20211018_MigrDC_Verse_Reference_Prep.rds')
PDAC_DC <- saveRDS('humanPDAC_integrated_total_DC.rds')
PDAC_mreg <- subset(PDAC_DC, subset = RNA_snn_res.0.6 %in% c(5,8))

anchors <- FindTransferAnchors(
  reference = ref_mreg,
  query = PDAC_mreg,
  normalization.method = "LogNormalize",
  dims = 1:50,
  reference.reduction = 'pca'
)

PDAC_mreg <- MapQuery(
  anchorset = anchors,
  query = PDAC_mreg,
  reference = ref_mreg,
  refdata = list(
    Subset = "Subset"
  ),
)

DimPlot(PDAC_mreg, reduction = "umap", group.by = "predicted.Subset", label = TRUE, label.size = 3, repel = TRUE)

counts_matrix_full <- as.data.frame(PDAC_mreg[["RNA"]]$data) 

umap = t(PDAC_mreg[["umap"]]@cell.embeddings)
dataset=PDAC_mreg$orig.ident
patient=PDAC_mreg$orig.ident
cluster06 = PDAC_mreg$RNA_snn_res.0.6
tissu=PDAC_mreg$orig.ident
mrg.subset =PDAC_mreg$predicted.Subset

total = rbind(counts_matrix_full,umap,dataset,patient,tissu,cluster06,mrg.subset)
rownames(total)[rownames(total) == "44271"] <- "Dataset"
rownames(total)[rownames(total) == "44272"] <- "Patient"
rownames(total)[rownames(total) == "44273"] <- "Tissue"
rownames(total)[rownames(total) == "44274"] <- "Cluster.res0.6"
rownames(total)[rownames(total) == "44275"] <- "Cluster.res0.8"
write.table(SG_Header, file = "humanPDAC_integrated_total_DC_full.txt",sep="\t")
suppressWarnings(write.table(total, file = "humanPDAC_integrated_total_DC_full.txt",sep="\t",append=T))