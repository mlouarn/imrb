source('/home/marine-louarn/imrb/packages.R')
setwd("/home/marine-louarn/Documents/PDAC/")

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

rename_org <- function(seurat_obj, rds_file){
  seurat_obj$orig.ident = str_extract(rds_file, ".*(?=\\.rds$)")
  return(seurat_obj)
}

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

seqgeq_file <- function(seurat_obj, list_genes, file_out){
  SG_Header <- read.table("~/Documents/Header_SeqGeq.txt",sep="\t",header=F,row.names=1,check.names = F)
  SG_Header <- unname(SG_Header)
  
  counts_matrix <- as.data.frame(seurat_obj[["RNA"]]$data) %>%
    filter(row.names(.) %in% list_genes)
  umap = t(seurat_obj[["umap"]]@cell.embeddings)
  total = rbind(counts_matrix,umap)
  
  write.table(SG_Header, file = file_out,sep="\t")
  suppressWarnings(write.table(total, file = file_out,sep="\t",append=T))
}

signature_human = read.csv("~/Documents/PDAC/20260317_Human_Signature_Genes_for_scRNAseq_data_Extraction.csv")
signature_human = signature_human$Vizgen_Gene

#Disco
disco = readRDS("Disco/Disco_integrated.rds")
disco$orig.ident = paste0("Disco_",disco$sample_id)
saveRDS(disco,"Disco/Disco_integrated.rds")
seqgeq_file(disco,signature_human,"Disco/Disco_rna_data.txt")

#scBroad
scBroad = readRDS("Datasets/scBroad/PDAC_RPCA_integrated_refined_annotation-001.rds")
scBroad = UpdateSeuratObject(scBroad)
seqgeq_file(scBroad,signature_human,"scBroad/scBroad_integrated.txt")
saveRDS(scBroad,"Datasets/scBroad/scBroad_integrated.rds")

#GSE134355
for (i in c('GSM4008637_Adult-Pancreas1','GSM4008704_Fetal-Pancreas1','GSM4008705_Fetal-Pancreas2','GSM4008706_Fetal-Pancreas3')){
  name = substr(i, 1, 10)
  df <- as.data.frame(fread(paste0('GSE134355/',i,'_dge.txt.gz')))
  rownames(df) <- df$GENE
  df$GENE <- NULL
  obj <- CreateSeuratObject(counts = df, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE134355/",name,".rds"))
}
rm(name,i,obj,df)
gc()

rds_files = list.files(path = 'GSE134355/', pattern = '*.rds$')
list_seurat = lapply(paste0('GSE134355/',rds_files), readRDS)
integrated = seurat_integration(list_seurat)
saveRDS(integrated,'GSE134355/GSE134355_integrated.rds')
seqgeq_file(integrated,signature_human,"GSE134355/GSE134355_integrated.txt")

#GSE217845
for (i in c('GSM6727546_PDAC_47','GSM6727547_PDAC_48','GSM6727548_PDAC_50','GSM6727549_PDAC_51','GSM6727550_PDAC_55','GSM6727551_PDAC_60','GSM6727542_LPDAC_15','GSM6727543_LPDAC_25','GSM6727544_LPDAC_26','GSM6727545_LPDAC_30')){
  obj <- ReadMtx(  mtx = paste0("GSE217845/",i,"_Tumor_matrix.mtx.gz"), features = paste0("GSE217845/",i,"_Tumor_features.tsv.gz"),
                   cells = paste0("GSE217845/",i,"_Tumor_barcodes.tsv.gz"),feature.column = 2)
  name = substr(i, 1, 10)
  colnames(obj) <- make.unique(colnames(obj))
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE217845/",name,".rds"))
}
rm(name,i,obj)
gc()
rds_files = list.files(path = 'GSE217845/', pattern = '*.rds$')
list_seurat = lapply(paste0('GSE217845/',rds_files), readRDS)
integrated = seurat_integration(list_seurat)
saveRDS(integrated,'GSE217845/GSE217845_integrated.rds')
seqgeq_file(integrated,signature_human,"GSE217845/GSE217845_integrated.txt")

#GSE247635
for (i in c('GSM7898197_Sample_4_8','GSM7898198_Sample_1_23','GSM7898199_Sample_5_30','GSM7898200_Sample_2_34','GSM7898201_Sample_6_62','GSM7898202_Sample_3_64')){
  obj <- ReadMtx(  mtx = paste0("GSE247635/",i,"_matrix.mtx.gz"), features = paste0("GSE247635/",i,"_features.tsv.gz"),
                   cells = paste0("GSE247635/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, 10)
  colnames(obj) <- make.unique(colnames(obj))
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE247635/",name,".rds"))
}
rm(name,i,obj)
gc()
rds_files = list.files(path = 'GSE247635/', pattern = '*.rds$')
list_seurat = lapply(paste0('GSE247635/',rds_files), readRDS)
integrated = seurat_integration(list_seurat)
saveRDS(integrated,'GSE247635/GSE247635_integrated.rds')
seqgeq_file(integrated,signature_human,"GSE247635/GSE247635_integrated.txt")

#GSE311788
for (i in c('GSM9332044_SC01','GSM9332045_SC02','GSM9332046_SC03','GSM9332047_SC04','GSM9332048_SC05','GSM9332049_SC06')){
  obj <- ReadMtx(  mtx = paste0("GSE311788/",i,"_matrix.mtx.gz"), features = paste0("GSE311788/",i,"_features.tsv.gz"),
                   cells = paste0("GSE311788/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, 10)
  colnames(obj) <- make.unique(colnames(obj))
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE311788/",name,".rds"))
}
rm(name,i,obj)
gc()
rds_files = list.files(path = 'GSE311788/', pattern = '*.rds$')
list_seurat = lapply(paste0('GSE311788/',rds_files), readRDS)
integrated = seurat_integration(list_seurat)
saveRDS(integrated,'GSE311788/GSE311788_integrated.rds')
seqgeq_file(integrated,signature_human,"GSE311788/GSE311788_integrated.txt")

#GSE212966
for (i in c('GSM6567165_ADJ1','GSM6567166_ADJ2','GSM6567167_ADJ3','GSM6567169_ADJ4','GSM6567170_ADJ5','GSM6567171_ADJ6','GSM6567157_PDAC1','GSM6567159_PDAC2','GSM6567160_PDAC3','GSM6567161_PDAC4','GSM6567163_PDAC5','GSM6567164_PDAC6')){
  obj <- ReadMtx(  mtx = paste0("GSE212966/",i,"_matrix.mtx.gz"), features = paste0("GSE212966/",i,"_genes.tsv.gz"),
                   cells = paste0("GSE212966/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, 10)
  colnames(obj) <- make.unique(colnames(obj))
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE212966/",name,".rds"))
}
rm(name,i,obj)
gc()
rds_files = list.files(path = 'GSE212966/', pattern = '*.rds$')
list_seurat = lapply(paste0('GSE212966/',rds_files), readRDS)
integrated = seurat_integration(list_seurat)
saveRDS(integrated,'GSE212966/GSE212966_integrated.rds')
seqgeq_file(integrated,signature_human,"GSE212966/GSE212966_integrated.txt")

#GSE231535
for (i in c("GSM7289739_PDAC1",'GSM7289740_PDAC2')){
  obj <- ReadMtx(  mtx = paste0("GSE231535/",i,"_matrix.mtx.gz"), features = paste0("GSE231535/",i,"_genes.tsv.gz"),
                   cells = paste0("GSE231535/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, 10)
  colnames(obj) <- make.unique(colnames(obj))
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE231535/",name,".rds"))
}
rm(name,i,obj)
gc()
rds_files = list.files(path = 'GSE231535/', pattern = '*.rds$')
list_seurat = lapply(paste0('GSE231535/',rds_files), readRDS)
integrated = seurat_integration(list_seurat)
saveRDS(integrated,'GSE231535/GSE231535_integrated.rds')
seqgeq_file(integrated,signature_human,"GSE231535/GSE231535_integrated.txt")

#GSE283206
for (i in c('GSM8657816_HT543','GSM8657817_HT544','GSM8657818_HT555','GSM8657820_HT692','GSM8657813_HT400_1','GSM8657813_HT400_2','GSM8657813_HT400_4','GSM8657814_HT434_1','GSM8657814_HT434_4','GSM8657815_HT442_1','GSM8657815_HT442_3','GSM8657819_HT561_3','GSM8657819_HT561_4')){
  obj <- ReadMtx(  mtx = paste0("GSE283206/",i,"_matrix.mtx.gz"), features = paste0("GSE283206/",i,"_features.tsv.gz"),
                   cells = paste0("GSE283206/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, 10)
  colnames(obj) <- make.unique(colnames(obj))
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE283206/",name,".rds"))
}
rm(name,i,obj)
rm(list_seurat, integrated,rds_files)
gc()
rds_files = list.files(path = 'GSE283206/', pattern = '*.rds$')
list_seurat = lapply(paste0('GSE283206/',rds_files), readRDS)
integrated = seurat_integration(list_seurat)
saveRDS(integrated,'GSE283206/GSE283206_integrated.rds')
seqgeq_file(integrated,signature_human,"GSE283206/GSE283206_integrated.txt")

#Zenodo 
library(reticulate)
use_condaenv("~/miniconda3/envs/spacialPy")
Zeonodo_3969339 <- sceasy::convertFormat('Zeonodo_3969339/StdWf1_PRJCA001063_CRC_besca2.raw.h5ad', from="anndata", to="seurat",
                                   outFile='Zeonodo_3969339/Zeonodo_3969339.rds')
Zeonodo_3969339[["percent.mt"]] <- PercentageFeatureSet(Zeonodo_3969339, pattern = "^mt-")
#Zeonodo_3969339 <- subset(Zeonodo_3969339, subset = nFeature_RNA > 200& percent.mt <= 20)
Zeonodo_3969339 <- NormalizeData(Zeonodo_3969339, normalization.method = "LogNormalize", scale.factor = 10000)
Zeonodo_3969339 <- FindVariableFeatures(Zeonodo_3969339, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Zeonodo_3969339)
Zeonodo_3969339 <- ScaleData(Zeonodo_3969339, features = all.genes)
Zeonodo_3969339 <- RunPCA(Zeonodo_3969339, features = VariableFeatures(object = Zeonodo_3969339))
Zeonodo_3969339 <- FindNeighbors(Zeonodo_3969339, dims = 1:10)
Zeonodo_3969339 <- FindClusters(Zeonodo_3969339, resolution = 1)
Zeonodo_3969339 <- RunUMAP(Zeonodo_3969339, dims = 1:10)

Zeonodo_3969339$orig.ident = paste0("Zenodo_",Zeonodo_3969339$Patient)
saveRDS(Zeonodo_3969339, file = "Zeonodo_3969339/Zeonodo_3969339_integrated.rds")
seqgeq_file(Zeonodo_3969339,signature_human,"Zeonodo_3969339/Zeonodo_3969339.txt")

