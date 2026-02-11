source('/home/marinelouarn/imrb/packages.R')
source('/home/marinelouarn/imrb/functions.R')

setwd(dir = "/home/marinelouarn/Documents/Alexandre_mouse/")

#GSE152277_RAW/
GSM4610558 <- ReadMtx(  mtx = "GSE152277_RAW/GSM4610558_A1_matrix.mtx.gz", features = "GSE152277_RAW/GSM4610558_A1_features.tsv.gz",
  cells = "GSE152277_RAW/GSM4610558_A1_barcodes.tsv.gz")
GSM4610558 <- CreateSeuratObject(counts = GSM4610558, project = "GSM4610558", min.cells = 3)
GSM4610558 <- seurat_analyse(GSM4610558)
saveRDS(GSM4610558, file = "GSE152277_RAW/GSM4610558.rds")

GSM4610559 <- ReadMtx(  mtx = "GSE152277_RAW/GSM4610559_A2_matrix.mtx.gz", features = "GSE152277_RAW/GSM4610559_A2_features.tsv.gz",
                        cells = "GSE152277_RAW/GSM4610559_A2_barcodes.tsv.gz")
GSM4610559 <- CreateSeuratObject(counts = GSM4610559, project = "GSM4610559", min.cells = 3)
GSM4610559 <- seurat_analyse(GSM4610559)
saveRDS(GSM4610559, file = "GSE152277_RAW/GSM4610559.rds")

GSM4610560 <- ReadMtx(  mtx = "GSE152277_RAW/GSM4610560_B1_matrix.mtx.gz", features = "GSE152277_RAW/GSM4610560_B1_features.tsv.gz",
                        cells = "GSE152277_RAW/GSM4610560_B1_barcodes.tsv.gz")
GSM4610560 <- CreateSeuratObject(counts = GSM4610560, project = "GSM4610560", min.cells = 3)
GSM4610560 <- seurat_analyse(GSM4610560)

saveRDS(GSM4610560, file = "GSE152277_RAW/GSM4610560.rds")

GSM4610561 <- ReadMtx(  mtx = "GSE152277_RAW/GSM4610561_B3_matrix.mtx.gz", features = "GSE152277_RAW/GSM4610561_B3_features.tsv.gz",
                        cells = "GSE152277_RAW/GSM4610561_B3_barcodes.tsv.gz")
GSM4610561 <- CreateSeuratObject(counts = GSM4610561, project = "GSM4610561", min.cells = 3)
GSM4610561 <- seurat_analyse(GSM4610561)
saveRDS(GSM4610561, file = "GSE152277_RAW/GSM4610561.rds")

GSM4610562 <- ReadMtx(  mtx = "GSE152277_RAW/GSM4610562_Tumor_B1_matrix.mtx.gz", features = "GSE152277_RAW/GSM4610562_Tumor_B1_features.tsv.gz",
                        cells = "GSE152277_RAW/GSM4610562_Tumor_B1_barcodes.tsv.gz")
GSM4610562 <- CreateSeuratObject(counts = GSM4610562, project = "GSM4610562", min.cells = 3)
GSM4610562 <- seurat_analyse(GSM4610562)
saveRDS(GSM4610562, file = "GSE152277_RAW/GSM4610562.rds")

GSM4610563 <- ReadMtx(  mtx = "GSE152277_RAW/GSM4610563_Tumor_B2_matrix.mtx.gz", features = "GSE152277_RAW/GSM4610563_Tumor_B2_features.tsv.gz",
                        cells = "GSE152277_RAW/GSM4610563_Tumor_B2_barcodes.tsv.gz")
GSM4610563 <- CreateSeuratObject(counts = GSM4610563, project = "GSM4610563", min.cells = 3)
GSM4610563 <- seurat_analyse(GSM4610563)
saveRDS(GSM4610563, file = "GSE152277_RAW/GSM4610563.rds")

###GSE216805
GSM6693193 <- Read10X_h5("GSE216805_RAW/GSM6693193_6647_Control_sample_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6693193 <- CreateSeuratObject(counts = GSM6693193, project = "GSM6693193", min.cells = 3)
GSM6693193 <- seurat_analyse(GSM6693193)
saveRDS(GSM6693193, file = "GSE216805_RAW/GSM6693193.rds")

GSM6693194 <- Read10X_h5("GSE216805_RAW/GSM6693194_6615_NrasG12D_sample_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6693194 <- CreateSeuratObject(counts = GSM6693194, project = "GSM6693193", min.cells = 3)
GSM6693194 <- seurat_analyse(GSM6693194)
saveRDS(GSM6693194, file = "GSE216805_RAW/GSM6693194.rds")


###GSE201247
GSM6056008 <- Read10X_h5("GSE201247_RAW/GSM6056008_1_ATTAC_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6056008 <- CreateSeuratObject(counts = GSM6056008, project = "GSM6056008", min.cells = 3)
GSM6056008 <- seurat_analyse(GSM6056008)
saveRDS(GSM6056008, file = "GSE201247_RAW/GSM6056008.rds")

GSM6056009 <- Read10X_h5("GSE201247_RAW/GSM6056009_2_Kras_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6056009 <- CreateSeuratObject(counts = GSM6056009, project = "GSM6056009", min.cells = 3)
GSM6056009 <- seurat_analyse(GSM6056009)
saveRDS(GSM6056009, file = "GSE201247_RAW/GSM6056009.rds")

GSM6056010 <- Read10X_h5("GSE201247_RAW/GSM6056010_3_ATTAC_Kras_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6056010 <- CreateSeuratObject(counts = GSM6056010, project = "GSM6056010", min.cells = 3)
GSM6056010 <- seurat_analyse(GSM6056010)
saveRDS(GSM6056010, file = "GSE201247_RAW/GSM6056010.rds")

GSM6056011 <- Read10X_h5("GSE201247_RAW/GSM6056011_4_ATTAC_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6056011 <- CreateSeuratObject(counts = GSM6056011, project = "GSM6056011", min.cells = 3)
GSM6056011 <- seurat_analyse(GSM6056011)
saveRDS(GSM6056011, file = "GSE201247_RAW/GSM6056011.rds")

GSM6056012 <- Read10X_h5("GSE201247_RAW/GSM6056012_5_Kras_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6056012 <- CreateSeuratObject(counts = GSM6056012, project = "GSM6056012", min.cells = 3)
GSM6056012 <- seurat_analyse(GSM6056012)
saveRDS(GSM6056012, file = "GSE201247_RAW/GSM6056012.rds")

GSM6056013 <- Read10X_h5("GSE201247_RAW/GSM6056013_6_ATTAC_Kras_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
GSM6056013 <- CreateSeuratObject(counts = GSM6056013, project = "GSM6056013", min.cells = 3)
GSM6056013 <- seurat_analyse(GSM6056013)
saveRDS(GSM6056013, file = "GSE201247_RAW/GSM6056013.rds")

###GSE217368 : Not touched
GSE217368 <- readRDS('GSE217368/GSE217368_seurat_tdTom_allSamples_ForSubmission.rds')

###GSE127465
df <- as.data.frame(fread('GSE127465_RAW/GSE127465_mouse_cell_metadata_15939x12.tsv.gz'))
df$unique_barcode = paste0(df$Barcode,'_',df$Library)
write.csv(df$unique_barcode,'GSE127465_RAW/barcodes.csv',col.names=FALSE,row.names=FALSE)

GSE127465 <- ReadMtx(  mtx = "GSE127465_RAW/GSE127465_mouse_counts_normalized_15939x28205.mtx.gz", features = "GSE127465_RAW/GSE127465_gene_names_mouse_28205.tsv.gz",
                        cells ='GSE127465_RAW/barcodes.csv' , feature.column = 1,
                       mtx.transpose = TRUE,skip.cell=1)
GSE127465 <- CreateSeuratObject(counts = GSE127465, project = "GSE127465", min.cells = 3)

GSE127465 <-AddMetaData(GSE127465, df, col.name = NULL)
GSE127465 <- seurat_analyse(GSE127465)
saveRDS(GSE127465, file = "GSE127465_RAW/GSE127465.rds")

###GSE189856

sceasy::convertFormat('GSE189856/GSE189856_Full_Data_Processed_adata.h5ad', from="anndata", to="seurat",
                      outFile='GSE189856.rds')
