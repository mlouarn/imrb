library(Seurat,future)
library(tidyverse)
library(dplyr)

#plan(workers = 2)
options(future.globals.maxSize = 8000 * 1024^2)

signatures_mouse = read.csv('signature.v1.csv')
signatures_mouse = signatures_mouse[signatures_mouse$Keep_for_initial_screening=='y',]
signatures_mouse_b = read.csv('signature_boissonnas.v1.csv')

GSE134904<-readRDS( "GSE134904_myeloid.rds")
GSE145502<-readRDS( "GSE145502_myeloid.rds")

list_datasets = c('GSE148445','GSE150675','GSE152674','GSE155011','GSE158677','GSE159478','GSE160641','GSE161769','GSE163508','GSE184096','GSE192935','GSE217847','GSE221528','GSE229355','GSE229765','GSE244582','GSE252083','GSE263963')#,'GSE264232','GSE275907','GSE283034','GSE283609','GSE297521','GSE297920','GSE304859','GSE307750')

print('itsgoing')
myelo <- merge(GSE134904, y = c(GSE145502), project = "integrated_myelo")
rm(GSE134904,GSE145502)

for (dataset in list_datasets){
  print(dataset)
  seurat_dataset <- readRDS(paste0(dataset,'_myeloid.rds'))
  if (length(colnames(seurat_dataset)) > 20000){
    seurat_dataset = subset(seurat_dataset, cells.use = sample(x = colnames(seurat_dataset), size = 20000) )
  }
  myelo <- merge(myelo, seurat_dataset, project = "integrated_myelo")
  rm(seurat_dataset)
  gc()
}

print('notded')
saveRDS(myelo, "Merge_myelo.rds")
print('checkpoint,saved')
myelo <- readRDS("Merge_myelo.rds")
myelo[['percent.mt']] <- PercentageFeatureSet(myelo, pattern = "^mt-")

myelo <- SCTransform(myelo, vars.to.regress = "percent.mt", verbose = FALSE)
gc()
myelo <- RunPCA(myelo)
myelo <- FindNeighbors(myelo, dims = 1:10, reduction = "pca")
myelo <- RunUMAP(myelo, dims = 1:10, reduction = "pca", reduction.name = "umap.unintegrated")

myelo <- IntegrateLayers(object = myelo, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",assay = "SCT", verbose = FALSE)
myelo[["RNA"]] <- JoinLayers(myelo[["RNA"]])
myelo <- FindNeighbors(myelo, reduction = "harmony", dims = 1:30)
myelo <- RunUMAP(myelo, dims = 1:30, reduction = "harmony")


DimPlot(myelo, reduction = "umap",label = T)
DimPlot(myelo, reduction = "umap", group.by = "orig.ident")

saveRDS(myelo, "Integrated_myelo.rds")
#myelo = readRDS("Integrated_myelo.rds")


SG_Header <- read.table("../pymt_scrape/Header_SeqGeq.txt",sep="\t",header=F,row.names=1,check.names = F)
SG_Header <- unname(SG_Header)

genes = unique(c(signatures_mouse$Signature_gene_mouse.style_symbol,signatures_mouse_b$Signature_genes))
counts_matrix <- as.data.frame(myelo[["SCT"]]$data) %>%
  filter(row.names(.) %in% genes)

umap = t(myelo[["umap"]]@cell.embeddings)
cluster = t(myelo[["seurat_clusters"]])
pc=t(myelo@reductions$pca@cell.embeddings)
myelo <- RunICA(myelo)
ica=t(myelo@reductions$ica@cell.embeddings)

sample=myelo$orig.ident
dataset=myelo$orig.ident

total = rbind( counts_matrix , umap, cluster,pc,ica,sample,dataset)
colnames(total)[colnames(total) == '297'] <- 'sample'
colnames(total)[colnames(total) == '298'] <- 'dataset'


write.table(SG_Header, file = "Integrated_myelo.txt",sep="\t")
suppressWarnings(write.table(total, file = "Integrated_myelo.txt",sep="\t",append=T))
