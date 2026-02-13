source('/home/marinelouarn/imrb/packages.R')
source('/home/marinelouarn/imrb/functions.R')

setwd(dir = "/home/marinelouarn/Documents/Alexandre_mouse/")

#GSE152277/
GSM4610558 <- ReadMtx(  mtx = "GSE152277/GSM4610558_A1_matrix.mtx.gz", features = "GSE152277/GSM4610558_A1_features.tsv.gz",
  cells = "GSE152277/GSM4610558_A1_barcodes.tsv.gz")
GSM4610558 <- CreateSeuratObject(counts = GSM4610558, project = "GSM4610558", min.cells = 3)
GSM4610558 <- seurat_analyse(GSM4610558)
saveRDS(GSM4610558, file = "GSE152277/GSM4610558.rds")

GSM4610559 <- ReadMtx(  mtx = "GSE152277/GSM4610559_A2_matrix.mtx.gz", features = "GSE152277/GSM4610559_A2_features.tsv.gz",
                        cells = "GSE152277/GSM4610559_A2_barcodes.tsv.gz")
GSM4610559 <- CreateSeuratObject(counts = GSM4610559, project = "GSM4610559", min.cells = 3)
GSM4610559 <- seurat_analyse(GSM4610559)
saveRDS(GSM4610559, file = "GSE152277/GSM4610559.rds")

GSM4610560 <- ReadMtx(  mtx = "GSE152277/GSM4610560_B1_matrix.mtx.gz", features = "GSE152277/GSM4610560_B1_features.tsv.gz",
                        cells = "GSE152277/GSM4610560_B1_barcodes.tsv.gz")
GSM4610560 <- CreateSeuratObject(counts = GSM4610560, project = "GSM4610560", min.cells = 3)
GSM4610560 <- seurat_analyse(GSM4610560)

saveRDS(GSM4610560, file = "GSE152277/GSM4610560.rds")

GSM4610561 <- ReadMtx(  mtx = "GSE152277/GSM4610561_B3_matrix.mtx.gz", features = "GSE152277/GSM4610561_B3_features.tsv.gz",
                        cells = "GSE152277/GSM4610561_B3_barcodes.tsv.gz")
GSM4610561 <- CreateSeuratObject(counts = GSM4610561, project = "GSM4610561", min.cells = 3)
GSM4610561 <- seurat_analyse(GSM4610561)
saveRDS(GSM4610561, file = "GSE152277/GSM4610561.rds")

GSM4610562 <- ReadMtx(  mtx = "GSE152277/GSM4610562_Tumor_B1_matrix.mtx.gz", features = "GSE152277/GSM4610562_Tumor_B1_features.tsv.gz",
                        cells = "GSE152277/GSM4610562_Tumor_B1_barcodes.tsv.gz")
GSM4610562 <- CreateSeuratObject(counts = GSM4610562, project = "GSM4610562", min.cells = 3)
GSM4610562 <- seurat_analyse(GSM4610562)
saveRDS(GSM4610562, file = "GSE152277/GSM4610562.rds")

GSM4610563 <- ReadMtx(  mtx = "GSE152277/GSM4610563_Tumor_B2_matrix.mtx.gz", features = "GSE152277/GSM4610563_Tumor_B2_features.tsv.gz",
                        cells = "GSE152277/GSM4610563_Tumor_B2_barcodes.tsv.gz")
GSM4610563 <- CreateSeuratObject(counts = GSM4610563, project = "GSM4610563", min.cells = 3)
GSM4610563 <- seurat_analyse(GSM4610563)
saveRDS(GSM4610563, file = "GSE152277/GSM4610563.rds")

###GSE216805
GSM6693193 <- Read10X_h5("GSE216805/GSM6693193_6647_Control_sample_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6693193 <- CreateSeuratObject(counts = GSM6693193, project = "GSM6693193", min.cells = 3)
GSM6693193 <- seurat_analyse(GSM6693193)
saveRDS(GSM6693193, file = "GSE216805/GSM6693193.rds")

GSM6693194 <- Read10X_h5("GSE216805/GSM6693194_6615_NrasG12D_sample_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6693194 <- CreateSeuratObject(counts = GSM6693194, project = "GSM6693193", min.cells = 3)
GSM6693194 <- seurat_analyse(GSM6693194)
saveRDS(GSM6693194, file = "GSE216805/GSM6693194.rds")


###GSE201247
GSM6056008 <- Read10X_h5("GSE201247/GSM6056008_1_ATTAC_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6056008 <- CreateSeuratObject(counts = GSM6056008, project = "GSM6056008", min.cells = 3)
GSM6056008 <- seurat_analyse(GSM6056008)
saveRDS(GSM6056008, file = "GSE201247/GSM6056008.rds")

GSM6056009 <- Read10X_h5("GSE201247/GSM6056009_2_Kras_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6056009 <- CreateSeuratObject(counts = GSM6056009, project = "GSM6056009", min.cells = 3)
GSM6056009 <- seurat_analyse(GSM6056009)
saveRDS(GSM6056009, file = "GSE201247/GSM6056009.rds")

GSM6056010 <- Read10X_h5("GSE201247/GSM6056010_3_ATTAC_Kras_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6056010 <- CreateSeuratObject(counts = GSM6056010, project = "GSM6056010", min.cells = 3)
GSM6056010 <- seurat_analyse(GSM6056010)
saveRDS(GSM6056010, file = "GSE201247/GSM6056010.rds")

GSM6056011 <- Read10X_h5("GSE201247/GSM6056011_4_ATTAC_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6056011 <- CreateSeuratObject(counts = GSM6056011, project = "GSM6056011", min.cells = 3)
GSM6056011 <- seurat_analyse(GSM6056011)
saveRDS(GSM6056011, file = "GSE201247/GSM6056011.rds")

GSM6056012 <- Read10X_h5("GSE201247/GSM6056012_5_Kras_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6056012 <- CreateSeuratObject(counts = GSM6056012, project = "GSM6056012", min.cells = 3)
GSM6056012 <- seurat_analyse(GSM6056012)
saveRDS(GSM6056012, file = "GSE201247/GSM6056012.rds")

GSM6056013 <- Read10X_h5("GSE201247/GSM6056013_6_ATTAC_Kras_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6056013 <- CreateSeuratObject(counts = GSM6056013, project = "GSM6056013", min.cells = 3)
GSM6056013 <- seurat_analyse(GSM6056013)
saveRDS(GSM6056013, file = "GSE201247/GSM6056013.rds")

###GSE217368 : Not touched
GSE217368 <- readRDS('GSE217368/GSE217368_seurat_tdTom_allSamples_ForSubmission.rds')

###GSE127465
df <- as.data.frame(fread('GSE127465/GSE127465_mouse_cell_metadata_15939x12.tsv.gz'))
df$unique_barcode = paste0(df$Barcode,'_',df$Library)
write.csv(df$unique_barcode,'GSE127465/barcodes.csv',col.names=FALSE,row.names=FALSE)

GSE127465 <- ReadMtx(  mtx = "GSE127465/GSE127465_mouse_counts_normalized_15939x28205.mtx.gz", features = "GSE127465/GSE127465_gene_names_mouse_28205.tsv.gz",
                        cells ='GSE127465/barcodes.csv' , feature.column = 1,
                       mtx.transpose = TRUE,skip.cell=1)
GSE127465 <- CreateSeuratObject(counts = GSE127465, project = "GSE127465", min.cells = 3)

GSE127465 <-AddMetaData(GSE127465, df, col.name = NULL)
GSE127465 <- seurat_analyse(GSE127465)
saveRDS(GSE127465, file = "GSE127465/GSE127465.rds")

###GSE189856
library(reticulate)
use_virtualenv("~/my-venv")
GSE189856 <- sceasy::convertFormat('GSE189856/GSE189856_Full_Data_Processed_adata.h5ad', from="anndata", to="seurat",
                      outFile='GSE189856.rds')

GSE189856[["percent.mt"]] <- PercentageFeatureSet(GSE189856, pattern = "^mt-")
#GSE189856 <- subset(GSE189856, subset = nFeature_RNA > 200& percent.mt <= 20)
GSE189856 <- NormalizeData(GSE189856, normalization.method = "LogNormalize", scale.factor = 10000)
GSE189856 <- FindVariableFeatures(GSE189856, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GSE189856)
GSE189856 <- ScaleData(GSE189856, features = all.genes)
GSE189856 <- RunPCA(GSE189856, features = VariableFeatures(object = GSE189856))
GSE189856 <- FindNeighbors(GSE189856, dims = 1:10)
GSE189856 <- FindClusters(GSE189856, resolution = 1)
GSE189856 <- RunUMAP(GSE189856, dims = 1:10)
saveRDS(GSE189856, file = "GSE189856/GSE189856.rds")

###GSE217847
GSM6727558 <- ReadMtx(  mtx = "GSE217847/GSM6727558_Tumor_DAY_10_matrix.mtx.gz", features = "GSE217847/GSM6727558_Tumor_DAY_10_features.tsv.gz",
                        cells = "GSE217847/GSM6727558_Tumor_DAY_10_barcodes.tsv.gz")
GSM6727558 <- CreateSeuratObject(counts = GSM6727558, project = "GSM6727558", min.cells = 3)
GSM6727558 <- seurat_analyse(GSM6727558)
saveRDS(GSM6727558, file = "GSE217847/GSM6727558.rds")

GSM6727559 <- ReadMtx(  mtx = "GSE217847/GSM6727559_Tumor_DAY_20_matrix.mtx.gz", features = "GSE217847/GSM6727559_Tumor_DAY_20_features.tsv.gz",
                        cells = "GSE217847/GSM6727559_Tumor_DAY_20_barcodes.tsv.gz")
GSM6727559 <- CreateSeuratObject(counts = GSM6727559, project = "GSM6727559", min.cells = 3)
GSM6727559 <- seurat_analyse(GSM6727559)
saveRDS(GSM6727559, file = "GSE217847/GSM6727559.rds")

GSM6727560 <- ReadMtx(  mtx = "GSE217847/GSM6727560_Tumor_DAY_30_matrix.mtx.gz", features = "GSE217847/GSM6727560_Tumor_DAY_30_features.tsv.gz",
                        cells = "GSE217847/GSM6727560_Tumor_DAY_30_barcodes.tsv.gz")
GSM6727560 <- CreateSeuratObject(counts = GSM6727560, project = "GSM6727560", min.cells = 3)
GSM6727560 <- seurat_analyse(GSM6727560)
saveRDS(GSM6727560, file = "GSE217847/GSM6727560.rds")

GSM6727561 <- ReadMtx(  mtx = "GSE217847/GSM6727561_Healthy_Pancreas_matrix.mtx.gz", features = "GSE217847/GSM6727561_Healthy_Pancreas_features.tsv.gz",
                        cells = "GSE217847/GSM6727561_Healthy_Pancreas_barcodes.tsv.gz")
GSM6727561 <- CreateSeuratObject(counts = GSM6727561, project = "GSM6727561", min.cells = 3)
GSM6727561 <- seurat_analyse(GSM6727561)
saveRDS(GSM6727561, file = "GSE217847/GSM6727561.rds")

GSM6727566 <- ReadMtx(  mtx = "GSE217847/GSM6727566_Tumor_DAY_7_WT_matrix.mtx.gz", features = "GSE217847/GSM6727566_Tumor_DAY_7_WT_features.tsv.gz",
                        cells = "GSE217847/GSM6727566_Tumor_DAY_7_WT_barcodes.tsv.gz")
GSM6727566 <- CreateSeuratObject(counts = GSM6727566, project = "GSM6727566", min.cells = 3)
GSM6727566 <- seurat_analyse(GSM6727566)
saveRDS(GSM6727566, file = "GSE217847/GSM6727566.rds")

GSM6727567 <- ReadMtx(  mtx = "GSE217847/GSM6727567_Tumor_DAY_7_COX2_KO_matrix.mtx.gz", features = "GSE217847/GSM6727567_Tumor_DAY_7_COX2_KO_features.tsv.gz",
                        cells = "GSE217847/GSM6727567_Tumor_DAY_7_COX2_KO_barcodes.tsv.gz")
GSM6727567 <- CreateSeuratObject(counts = GSM6727567, project = "GSM6727567", min.cells = 3)
GSM6727567 <- seurat_analyse(GSM6727567)
saveRDS(GSM6727567, file = "GSE217847/GSM6727567.rds")

GSM6727568 <- ReadMtx(  mtx = "GSE217847/GSM6727568_KPC_GEMM_Sample1_matrix.mtx.gz", features = "GSE217847/GSM6727568_KPC_GEMM_Sample1_features.tsv.gz",
                        cells = "GSE217847/GSM6727568_KPC_GEMM_Sample1_barcodes.tsv.gz")
GSM6727568 <- CreateSeuratObject(counts = GSM6727568, project = "GSM6727568", min.cells = 3)
GSM6727568 <- seurat_analyse(GSM6727568)
saveRDS(GSM6727568, file = "GSE217847/GSM6727568.rds")

GSM6727569 <- ReadMtx(  mtx = "GSE217847/GSM6727569_KPC_GEMM_Sample2_matrix.mtx.gz", features = "GSE217847/GSM6727569_KPC_GEMM_Sample2_features.tsv.gz",
                        cells = "GSE217847/GSM6727569_KPC_GEMM_Sample2_barcodes.tsv.gz")
GSM6727569 <- CreateSeuratObject(counts = GSM6727569, project = "GSM6727569", min.cells = 3)
GSM6727569 <- seurat_analyse(GSM6727569)
saveRDS(GSM6727569, file = "GSE217847/GSM6727569.rds")

GSM6727570 <- ReadMtx(  mtx = "GSE217847/GSM6727570_KC_Ortho_Sample1_matrix.mtx.gz", features = "GSE217847/GSM6727570_KC_Ortho_Sample1_genes.tsv.gz",
                        cells = "GSE217847/GSM6727570_KC_Ortho_Sample1_barcodes.tsv.gz")
GSM6727570 <- CreateSeuratObject(counts = GSM6727570, project = "GSM6727570", min.cells = 3)
GSM6727570 <- seurat_analyse(GSM6727570)
saveRDS(GSM6727570, file = "GSE217847/GSM6727570.rds")

GSM6727571 <- ReadMtx(  mtx = "GSE217847/GSM6727571_KC_Ortho_Sample2_matrix.mtx.gz", features = "GSE217847/GSM6727571_KC_Ortho_Sample2_genes.tsv.gz",
                        cells = "GSE217847/GSM6727571_KC_Ortho_Sample2_barcodes.tsv.gz")
GSM6727571 <- CreateSeuratObject(counts = GSM6727571, project = "GSM6727571", min.cells = 3)
GSM6727571 <- seurat_analyse(GSM6727571)
saveRDS(GSM6727571, file = "GSE217847/GSM6727571.rds")

GSM6727572 <- ReadMtx(  mtx = "GSE217847/GSM6727572_KC_Hetero_Sample1_matrix.mtx.gz", features = "GSE217847/GSM6727572_KC_Hetero_Sample1_genes.tsv.gz",
                        cells = "GSE217847/GSM6727572_KC_Hetero_Sample1_barcodes.tsv.gz")
GSM6727572 <- CreateSeuratObject(counts = GSM6727572, project = "GSM6727572", min.cells = 3)
GSM6727572 <- seurat_analyse(GSM6727572)
saveRDS(GSM6727572, file = "GSE217847/GSM6727572.rds")

GSM6727573 <- ReadMtx(  mtx = "GSE217847/GSM6727573_KC_Hetero_Sample2_matrix.mtx.gz", features = "GSE217847/GSM6727573_KC_Hetero_Sample2_genes.tsv.gz",
                        cells = "GSE217847/GSM6727573_KC_Hetero_Sample2_barcodes.tsv.gz")
GSM6727573 <- CreateSeuratObject(counts = GSM6727573, project = "GSM6727573", min.cells = 3)
GSM6727573 <- seurat_analyse(GSM6727573)
saveRDS(GSM6727573, file = "GSE217847/GSM6727573.rds")

###GSE155275
dt = t(fread("GSE155275/GSE155275_C1_RSEC_MolsPerCell.csv.gz"))
colnames(dt) <- as.character(dt[1, ])
dt <- dt[-1,]
GSE155275 <- CreateSeuratObject(counts = dt, project = "GSE155275", min.cells = 3)
GSE155275 <- seurat_analyse(GSE155275)
saveRDS(GSE155275, file = "GSE155275/GSE155275.rds")

###GSE240148
GSE240148_a1 <- ReadMtx(  mtx = "GSE240148/GSE240148_sample_a1_matrix.mtx.gz", features = "GSE240148/GSE240148_sample_a1_features.tsv.gz",
                        cells = "GSE240148/GSE240148_sample_a1_barcodes.tsv.gz")
GSE240148_a1 <- CreateSeuratObject(counts = GSE240148_a1, project = "GSE240148_a1", min.cells = 3)
GSE240148_a1 <- seurat_analyse(GSE240148_a1)
saveRDS(GSE240148_a1, file = "GSE240148/GSE240148_a1.rds")

GSE240148_a4 <- ReadMtx(  mtx = "GSE240148/GSE240148_sample_a4_matrix.mtx.gz", features = "GSE240148/GSE240148_sample_a4_features.tsv.gz",
                          cells = "GSE240148/GSE240148_sample_a4_barcodes.tsv.gz")
GSE240148_a4 <- CreateSeuratObject(counts = GSE240148_a4, project = "GSE240148_a4", min.cells = 3)
GSE240148_a4 <- seurat_analyse(GSE240148_a4)
saveRDS(GSE240148_a4, file = "GSE240148/GSE240148_a4.rds")

GSE240148_c1_1 <- ReadMtx(  mtx = "GSE240148/GSE240148_sample_c1_1_matrix.mtx.gz", features = "GSE240148/GSE240148_sample_c1_1_features.tsv.gz",
                            cells = "GSE240148/GSE240148_sample_c1_1_barcodes.tsv.gz")
GSE240148_c1_1 <- CreateSeuratObject(counts = GSE240148_c1_1, project = "GSE240148_c1_1", min.cells = 3)
GSE240148_c1_1 <- seurat_analyse(GSE240148_c1_1)
saveRDS(GSE240148_c1_1, file = "GSE240148/GSE240148_c1_1.rds")

GSE240148_c1_2 <- ReadMtx(  mtx = "GSE240148/GSE240148_sample_c1_2_matrix.mtx.gz", features = "GSE240148/GSE240148_sample_c1_2_features.tsv.gz",
                            cells = "GSE240148/GSE240148_sample_c1_2_barcodes.tsv.gz")
GSE240148_c1_2 <- CreateSeuratObject(counts = GSE240148_c1_2, project = "GSE240148_c1_2", min.cells = 3)
GSE240148_c1_2 <- seurat_analyse(GSE240148_c1_2)
saveRDS(GSE240148_c1_2, file = "GSE240148/GSE240148_c1_2.rds")

GSE240148_d1 <- ReadMtx(  mtx = "GSE240148/GSE240148_sample_d1_matrix.mtx.gz", features = "GSE240148/GSE240148_sample_d1_features.tsv.gz",
                          cells = "GSE240148/GSE240148_sample_d1_barcodes.tsv.gz")
GSE240148_d1 <- CreateSeuratObject(counts = GSE240148_d1, project = "GSE240148_d1", min.cells = 3)
GSE240148_d1 <- seurat_analyse(GSE240148_d1)
saveRDS(GSE240148_d1, file = "GSE240148/GSE240148_d1.rds")

GSE240148_d4 <- ReadMtx(  mtx = "GSE240148/GSE240148_sample_d4_matrix.mtx.gz", features = "GSE240148/GSE240148_sample_d4_features.tsv.gz",
                          cells = "GSE240148/GSE240148_sample_d4_barcodes.tsv.gz")
GSE240148_d4 <- CreateSeuratObject(counts = GSE240148_d4, project = "GSE240148_d4", min.cells = 3)
GSE240148_d4 <- seurat_analyse(GSE240148_d4)
saveRDS(GSE240148_d4, file = "GSE240148/GSE240148_d4.rds")

GSE240148_f1 <- ReadMtx(  mtx = "GSE240148/GSE240148_sample_f1_matrix.mtx.gz", features = "GSE240148/GSE240148_sample_f1_features.tsv.gz",
                          cells = "GSE240148/GSE240148_sample_f1_barcodes.tsv.gz")
GSE240148_f1 <- CreateSeuratObject(counts = GSE240148_f1, project = "GSE240148_f1", min.cells = 3)
GSE240148_f1 <- seurat_analyse(GSE240148_f1)
saveRDS(GSE240148_f1, file = "GSE240148/GSE240148_f1.rds")

GSE240148_f3 <- ReadMtx(  mtx = "GSE240148/GSE240148_sample_f3_matrix.mtx.gz", features = "GSE240148/GSE240148_sample_f3_features.tsv.gz",
                          cells = "GSE240148/GSE240148_sample_f3_barcodes.tsv.gz")
GSE240148_f3 <- CreateSeuratObject(counts = GSE240148_f3, project = "GSE240148_f3", min.cells = 3)
GSE240148_f3 <- seurat_analyse(GSE240148_f3)
saveRDS(GSE240148_f3, file = "GSE240148/GSE240148_f3.rds")

###GSE307143
for (i in c("GSM9217269_SS_20_0088","GSM9217270_SS_20_0089","GSM9217271_SS_20_0090","GSM9217272_SS_20_0335","GSM9217273_SS_20_0336","GSM9217274_SS_20_0337","GSM9217275_SS_20_0368","GSM9217276_SS_20_0369","GSM9217277_SS_20_0370","GSM9217278_SS_20_0449","GSM9217279_SS_20_0450","GSM9217280_SS_20_0451","GSM9217281_SS_20_0524","GSM9217282_SS_20_0525","GSM9217283_SS_20_0526","GSM9217284_SS_20_0527","GSM9217285_SS_20_0528","GSM9217286_SS_20_0529","GSM9217287_SS_20_0530","GSM9217288_SS_20_0531","GSM9217289_SS_20_0532","GSM9217290_SS_20_0578","GSM9217291_SS_20_0579","GSM9217292_SS_20_0580","GSM9217293_SS_20_0621","GSM9217294_SS_20_0622","GSM9217295_SS_20_0623","GSM9217296_SS_20_0744","GSM9217297_SS_20_0745","GSM9217298_SS_20_0746")){
  obj <- ReadMtx(  mtx = paste0("GSE307143/",i,"_matrix.mtx.gz"), features = paste0("GSE307143/",i,"_features.tsv.gz"),
                            cells = paste0("GSE307143/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, nchar(i)-11)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE307143/",name,".rds"))
}
rm(name,i)

###GSE291211
for (i in c("GSM8830804_T1","GSM8830805_T2","GSM8830806_T3","GSM8830807_T4","GSM8830808_T5","GSM8830809_T6","GSM8830810_T7","GSM8830811_T8","GSM8830812_T9","GSM8830813_T10")){
  obj <- ReadMtx(  mtx = paste0("GSE291211/",i,"_matrix.mtx.gz"), features = paste0("GSE291211/",i,"_features.tsv.gz"),
                   cells = paste0("GSE291211/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, nchar(i)-3)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE291211/",name,".rds"))
}
rm(name,i)
gc()

###GSE297920
for (i in c('GSM9001991_PyMT_Lung_NC_rep1','GSM9001992_PyMT_Lung_NC_rep2','GSM9001993_PyMT_Lung_CS_rep1','GSM9001994_PyMT_Lung_CS_rep2')){
  obj <- Read10X_h5(paste0("GSE297920/",i,".h5"), use.names = TRUE, unique.features = TRUE)
  name = substr(i, 1, nchar(i)-18)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE297920/",name,".rds"))
}
rm(name,i)
gc()

###GSE297521
for (i in c('GSM8993905_PyMT_Tumor_NC_rep1','GSM8993906_PyMT_Tumor_NC_rep2','GSM8993907_PyMT_Tumor_MS_rep1','GSM8993908_PyMT_Tumor_MS_rep2')){
  obj <- Read10X_h5(paste0("GSE297521/",i,"_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
  name = substr(i, 1, 10)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE297521/",name,".rds"))
}
rm(name,i)
gc()

###GSE275907
for (i in c('GSM8487918_PYMT_CXCR4_cKO_TUMOR','GSM8487919_PYMT_CXCR4_CTRL_TUMOR','GSM8767523_Normal_MG')){
  obj <- Read10X_h5(paste0("GSE275907/",i,"_filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
  name = substr(i, 1, 10)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE275907/",name,".rds"))
}
rm(name,i)
gc()


###GSE244582
for (i in c('GSM7821160_1','GSM7821161_8','GSM7821162_17','GSM7821163_21')){
  obj <- ReadMtx(  mtx = paste0("GSE244582/",i,"_matrix.mtx.gz"), features = paste0("GSE244582/",i,"_features.tsv.gz"),
                   cells = paste0("GSE244582/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, 10)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE244582/",name,".rds"))
}
rm(name,i)
gc()

###GSE283034
for (i in c('GSM8654897_s1.','GSM8654898_s2.','GSM8654899_s3.','GSM8654900_s4.','GSM8654901_s5.','GSM8654902_s6_')){
  df <- as.data.frame(t(fread(paste0('GSE283034/',i,'counts.csv.gz'))))
  colnames(df) <- as.character(df[1, ])
  df <- df[-1,]
  name = substr(i, 1, 10)  
  obj <- CreateSeuratObject(counts = df, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE283034/",name,".rds"))
}

###GSE264232
for (i in c('GSE264232_NTC-1','GSE264232_NTC-2','GSE264232_NTC-3','GSE264232_IL1a-2','GSE264232_IL1a-3')){
  obj <- ReadMtx(  mtx = paste0("GSE264232/",i,"_matrix.mtx.gz"), features = paste0("GSE264232/",i,"_features.tsv.gz"),
                   cells = paste0("GSE264232/",i,"_barcodes.tsv.gz"))
  #name = substr(i, 1, 10)
  obj <- CreateSeuratObject(counts = obj, project = i, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE264232/",i,".rds"))
}
rm(name,i)
gc()

###GSE229765
for (i in c("GSM7177495_PyMT_24h_BD_4186",  "GSM7177491_PyMT_Fresh_10x_rep1_2075",  "GSM7177492_PyMT_Fresh_10x_rep2_1925",  "GSM7177494_PyMT_Fresh_BD_3772")){
  obj <- ReadMtx(  mtx = paste0("GSE229765/",i,"cells_matrix.mtx.gz"), features = paste0("GSE229765/",i,"cells_genes.tsv.gz"),
                   cells = paste0("GSE229765/",i,"cells_barcodes.tsv.gz"),feature.column = 1)
  name = substr(i, 1, 10)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE229765/",name,".rds"))
}
rm(name,i,obj)
gc()

###GSE252083
for (i in c('GSM7993444_KO1','GSM7993445_KO2','GSM7993446_KO3','GSM7993447_WT1','GSM7993448_WT2')){
  obj <- ReadMtx(  mtx = paste0("GSE252083/",i,"_matrix.mtx.gz"), features = paste0("GSE252083/",i,"_genes.tsv.gz"),
                   cells = paste0("GSE252083/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, 10)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE252083/",name,".rds"))
}
rm(name,i,obj)
gc()

###GSE264124
for (i in c('GSM8212431_Sample1_Z46_filtered','GSM8212432_Sample2_Control_filtered')){
  obj <- ReadMtx(  mtx = paste0("GSE264124/",i,"_matrix.mtx.gz"), features = paste0("GSE264124/",i,"_features.tsv.gz"),
                   cells = paste0("GSE264124/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, 10)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE264124/",name,".rds"))
}
rm(name,i,obj)
gc()

###GSE270295 : Not touched
GSE270295 <- readRDS('GSE270295_all.tumor.sobj.rds')

###GSE159478
for (i in c('GSM4830541_T1','GSM4830542_T2','GSM4830543_T3','GSM4830544_T4')){
  obj <- ReadMtx(  mtx = paste0("GSE159478/",i,"_matrix.mtx.gz"), features = paste0("GSE159478/",i,"_genes.tsv.gz"),
                   cells = paste0("GSE159478/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, 10)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE159478/",name,".rds"))
}
rm(name,i,obj)
gc()

###GSE229355
for (i in c('GSM7159211_ABT_CD45_pos','GSM7159213_INK_CD45_pos','GSM7159215_Veh_CD45_pos','GSM7159217_WT_CD45_pos')){
  obj <- ReadMtx(  mtx = paste0("GSE229355/",i,"_matrix.mtx.gz"), features = paste0("GSE229355/",i,"_features.tsv.gz"),
                   cells = paste0("GSE229355/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, 10)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE229355/",name,".rds"))
}
rm(name,i,obj)
gc()

###GSE263963
for (i in c('GSM8208011_pymt_fresh_l1','GSM8208012_pymt_fresh_l2','GSM8208013_pymt_fresh_l3','GSM8208014_pymt_fresh_l4','GSM8208015_pymt_fresh_l5','GSM8208016_pymt_fresh_l6')){
  obj <- ReadMtx(  mtx = paste0("GSE263963/",i,"_matrix.mtx.gz"), features = paste0("GSE263963/",i,"_features.tsv.gz"),
                   cells = paste0("GSE263963/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, 10)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE263963/",name,".rds"))
}
rm(name,i,obj)
gc()

###GSE222524 : Not touched
GSM6925423 <- readRDS('GSE222524/GSM6925423_A_seurat.rds')
GSM6925424 <- readRDS('GSE222524/GSM6925424_B_seurat.rds')

###GSE161769
for (i in c('GSE161769_P-J','GSE161769_PYMT')){
  obj <- ReadMtx(  mtx = paste0("GSE161769/",i,"-matrix.mtx.gz"), features = paste0("GSE161769/",i,"-features.tsv.gz"),
                   cells = paste0("GSE161769/",i,"-barcodes.tsv.gz"))
  obj <- CreateSeuratObject(counts = obj, project = i, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE161769/",i,".rds"))
}
rm(i,obj)
gc()

###GSE155011
for (i in c('GSM4693144_16175X1','GSM4693145_16175X2')){
  obj <- ReadMtx(  mtx = paste0("GSE155011/",i,"_matrix.mtx.gz"), features = paste0("GSE155011/",i,"_features.tsv.gz"),
                   cells = paste0("GSE155011/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, 10)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE155011/",name,".rds"))
}
rm(name,i,obj)
gc()

###GSE150675
for (i in c('GSM4556249_GC_1','GSM4556250_GC_2','GSM4556251_Lit_1','GSM4556252_Lit_2','GSM4556245_Ctrl_1','GSM4556246_Ctrl_2','GSM4556247_Laser_1','GSM4556248_Laser_2')){
  obj <- ReadMtx(  mtx = paste0("GSE150675/",i,"_matrix.mtx.gz"), features = paste0("GSE150675/",i,"_features.tsv.gz"),
                   cells = paste0("GSE150675/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, 10)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE150675/",name,".rds"))
}
rm(name,i,obj)
gc()

###GSE160641
for (i in c('GSM4876355_pymt_tam_1','GSM4876356_pymt_tam_2','GSM4876357_pymt_tam_3')){
  obj <- ReadMtx(  mtx = paste0("GSE160641/",i,"_matrix.mtx.gz"), features = paste0("GSE160641/",i,"_features.tsv.gz"),
                   cells = paste0("GSE160641/",i,"_barcodes.tsv.gz"))
  name = substr(i, 1, 10)
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE160641/",name,".rds"))
}
rm(name,i,obj)
gc()

###GSE158677
for (i in c('GSM4805459_WT.T1','GSM4805460_WT.T2','GSM4805462_WT.T4','GSM4805454_ELF5.T2','GSM4805455_ELF5.T3','GSM4805456_ELF5.T4','GSM4805457_ELF5.T5','GSM4805458_ELF5.T6','GSM4805461_WT.T3','GSM4805463_WT.T5','GSM4805453_ELF5.T1')){
  obj <- ReadMtx(  mtx = paste0("GSE158677/",i,"_matrix.mtx.gz"), features = paste0("GSE158677/",i,"_features.tsv.gz"),
                   cells = paste0("GSE158677/",i,"_barcodes.tsv.gz"),feature.column = 1)
  name = substr(i, 1, 10)
  colnames(obj) <- make.unique(colnames(obj))
  obj <- CreateSeuratObject(counts = obj, project = name, min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = paste0("GSE158677/",name,".rds"))
}
rm(name,i,obj)
gc()

###GSE163508
  obj <- ReadMtx(  mtx = 'GSE163508/GSE163508_PYMT_aggr_matrix.mtx.gz', features = 'GSE163508/GSE163508_PYMT_aggr_features.tsv.gz',
                   cells = 'GSE163508/GSE163508_PYMT_aggr_barcodes.tsv.gz')
  obj <- CreateSeuratObject(counts = obj, project = 'GSE163508', min.cells = 3)
  obj <- seurat_analyse(obj)
  saveRDS(obj, file = 'GSE163508/GSE163508_PYMT.rds')

  ###GSE152674
  for (i in c('GSM4623139_WT','GSM4623140_KO')){
    name = substr(i, 1, 10)
    df <- as.data.frame(fread(paste0('GSE152674/',i,'_counts.txt.gz')))
    rownames(df) <- df$V1
    df$V1 <- NULL
    obj <- CreateSeuratObject(counts = df, project = name, min.cells = 3)
    obj <- seurat_analyse(obj)
    saveRDS(obj, file = paste0("GSE152674/",name,".rds"))
  }
  rm(name,i,obj)
  gc()
  
  ###GSE148445
  for (i in c('GSM4471590_PyMT1_CB','GSM4471591_PyMT2_CB','GSM4471592_PyMT3_CB','GSM4471593_PyMT4_CB','GSM4471594_WT_5','GSM4471595_WT6','GSM4471596_WT7','GSM4471597_WT8','GSM4471598_HK21_CB','GSM4471599_HK22_CB','GSM4471600_HK23_CB','GSM4471601_HK2_4','GSM4471602_HK2_5','GSM4471603_HK2_6','GSM4471604_HK27')){
    name = substr(i, 1, 10)
    df <- as.data.frame(fread(paste0('GSE148445/',i,'.dge.txt.gz')))
    rownames(df) <- df$GENE
    df$GENE <- NULL
    obj <- CreateSeuratObject(counts = df, project = name, min.cells = 3)
    obj <- seurat_analyse(obj)
    saveRDS(obj, file = paste0("GSE148445/",name,".rds"))
  }
  rm(name,i,obj)
  gc()
  
  ###GSE135096
  for (i in c('GSM3985845_Primary','GSM3985846_PYMT2','GSM3985847_PYMT3','GSM3985848_PYMT4','GSM3985849_PYMT6','GSM3985850_Secondary','GSM3985851_pMet1','GSM3985852_pMet2','GSM3985853_pAkt1KO','GSM3985854_Akt1B','GSM3985855_Akt1C','GSM3985856_pAkt2KO','GSM3985857_Akt2_B','GSM3985858_Akt2_C')){
    name = substr(i, 1, 10)
    df <- as.data.frame(fread(paste0('GSE135096/',i,'.dge.txt.gz')))
    rownames(df) <- df$GENE
    df$GENE <- NULL
    obj <- CreateSeuratObject(counts = df, project = name, min.cells = 3)
    obj <- seurat_analyse(obj)
    saveRDS(obj, file = paste0("GSE135096/",name,".rds"))
  }
  rm(name,i,obj,df)
  gc()
