source('/home/marinelouarn/imrb/packages.R')
source('/home/marinelouarn/imrb/functions.R')

setwd(dir = "/home/marinelouarn/Documents/Alexandre_mouse/PYMT\ scRNAseq/boissonas_2022_tam-breast_GSE184096/")

tam1.data <- Read10X(data.dir = "TAM1")
tam1 <- CreateSeuratObject(counts = tam1.data, project = "tam1", min.cells = 3)
tam1[["percent.mt"]] <- PercentageFeatureSet(tam1, pattern = "^mt-")
mtNan=NULL
mtNan[is.na(tam1$percent.mt)]='NaN_T'
mtNan[!is.na(tam1$percent.mt)]='NaN_F'
tam1 <- AddMetaData(tam1,mtNan,'mtNan')
tam1$percent.mt[is.na(tam1$percent.mt)]=0
VlnPlot(tam1, features = c('nFeature_RNA','percent.mt'),ncol = 2)

tam2.data <- Read10X(data.dir = "TAM2")
tam2 <- CreateSeuratObject(counts = tam2.data, project = "tam2", min.cells = 3)
tam2[["percent.mt"]] <- PercentageFeatureSet(tam2, pattern = "^mt-")
mtNan=NULL
mtNan[is.na(tam2$percent.mt)]='NaN_T'
mtNan[!is.na(tam2$percent.mt)]='NaN_F'
tam2 <- AddMetaData(tam2,mtNan,'mtNan')
tam2$percent.mt[is.na(tam2$percent.mt)]=0
VlnPlot(tam2, features = c('nFeature_RNA','percent.mt'),ncol = 2)

myelo.data <- Read10X(data.dir = "Bulk Myelo")
myelo <- CreateSeuratObject(counts = myelo.data, project = "myelo", min.cells = 3)
myelo[["percent.mt"]] <- PercentageFeatureSet(myelo, pattern = "^mt-")
VlnPlot(myelo, features = c('nFeature_RNA','percent.mt'),ncol = 2)

lympho.data <- Read10X(data.dir = "Bulk Lympho2")
lympho <- CreateSeuratObject(counts = lympho.data, project = "lympho", min.cells = 3)
lympho[["percent.mt"]] <- PercentageFeatureSet(lympho, pattern = "^mt-")
mtNan=NULL
mtNan[is.na(lympho$percent.mt)]='NaN_T'
mtNan[!is.na(lympho$percent.mt)]='NaN_F'
lympho <- AddMetaData(lympho,mtNan,'mtNan')
lympho$percent.mt[is.na(lympho$percent.mt)]=0
VlnPlot(lympho, features = c('nFeature_RNA','percent.mt'),ncol = 2)

#normalization
tam1_s <- subset(tam1, subset = nFeature_RNA > 200& percent.mt <= 20)
tam1_s <- NormalizeData(tam1_s, normalization.method = "LogNormalize", scale.factor = 10000)
tam1_s <- FindVariableFeatures(tam1_s, selection.method = "vst", nfeatures = 2000)

tam2_s <- subset(tam2, subset = nFeature_RNA > 200& percent.mt <= 20)
tam2_s <- NormalizeData(tam2_s, normalization.method = "LogNormalize", scale.factor = 10000)
tam2_s <- FindVariableFeatures(tam2_s, selection.method = "vst", nfeatures = 2000)


myelo_s <- subset(myelo, subset = nFeature_RNA > 200 & percent.mt <= 20)
myelo_s <- NormalizeData(myelo_s, normalization.method = "LogNormalize", scale.factor = 10000)
myelo_s <- FindVariableFeatures(myelo_s, selection.method = "vst", nfeatures = 2000)


lympho_s <- subset(lympho, subset = nFeature_RNA > 200& percent.mt <= 20)
lympho_s <- NormalizeData(lympho_s, normalization.method = "LogNormalize", scale.factor = 10000)
lympho_s <- FindVariableFeatures(lympho_s, selection.method = "vst", nfeatures = 2000)

#scaling
all.genes <- rownames(tam1_s)
tam1_s <- ScaleData(tam1_s, features = all.genes)

all.genes <- rownames(tam2_s)
tam2_s <- ScaleData(tam2_s, features = all.genes)

all.genes <- rownames(myelo_s)
myelo_s <- ScaleData(myelo_s, features = all.genes)

all.genes <- rownames(lympho_s)
lympho_s <- ScaleData(lympho_s, features = all.genes)

#PCA
tam1_s <- RunPCA(tam1_s, features = VariableFeatures(object = tam1_s))
tam2_s <- RunPCA(tam2_s, features = VariableFeatures(object = tam2_s))
myelo_s <- RunPCA(myelo_s, features = VariableFeatures(object = myelo_s))
lympho_s <- RunPCA(lympho_s, features = VariableFeatures(object = lympho_s))

#cluster
tam1_s <- FindNeighbors(tam1_s, dims = 1:10)
tam1_s <- FindClusters(tam1_s, resolution = 1)
tam1_s <- RunUMAP(tam1_s, dims = 1:10)
DimPlot(tam1_s, reduction = "umap", label = TRUE)

FeaturePlot(tam2_s,features=c('Epcam','Itgax','Cd74','Zbtb46','Cadm1','Xcr1'))
FeaturePlot(tam2_s,features=c('Lyz2','Folr2','Spp1','Trem2'))

tam2_s <- FindNeighbors(tam2_s, dims = 1:10)
tam2_s <- FindClusters(tam2_s, resolution = 1)
tam2_s <- RunUMAP(tam2_s, dims = 1:10)
DimPlot(tam2_s, reduction = "umap", label = TRUE)

FeaturePlot(tam2_s,features=c('Epcam','Itgax','Cd74','Zbtb46','Cadm1','Xcr1'))
FeaturePlot(tam2_s,features=c('Lyz2','Folr2','Spp1','Trem2'))

myelo_s <- FindNeighbors(myelo_s, dims = 1:10)
myelo_s <- FindClusters(myelo_s, resolution = 1)
myelo_s <- RunUMAP(myelo_s, dims = 1:10)
DimPlot(myelo_s, reduction = "umap", label = TRUE)

FeaturePlot(myelo_s,features=c('Epcam','Itgax','Cd74','Zbtb46','Cadm1','Xcr1'))
FeaturePlot(myelo_s,features=c('Lyz2','Folr2','Spp1','Trem2'))

lympho_s <- FindNeighbors(lympho_s, dims = 1:10)
lympho_s <- FindClusters(lympho_s, resolution = 1)
lympho_s <- RunUMAP(lympho_s, dims = 1:10)
DimPlot(lympho_s, reduction = "umap", label = TRUE)
#DimPlot(lympho_s, reduction = "umap",group.by = 'mtNan', label = TRUE)

FeaturePlot(lympho_s,features=c('Epcam','Itgax','Cd74','Zbtb46','Cadm1','Xcr1','Cd79a','Ms4a1','Cd19'))
FeaturePlot(lympho_s,features=c('Lyz2','Folr2','Spp1','Trem2'))
FeaturePlot(lympho_s,features=c('Cd3e','Cd4','Cd8a','Cd8b'))



#JOIN DATASETS
#integration layers
int.merge = merge(lympho, y=c(tam1,tam2,myelo),add.cell.ids=c('lympho','tam1','tam2','myelo'))
#int.combined <- SCTransform(int.merge, vars.to.regress = "percent.mt", verbose = FALSE)
int.merge <- subset(int.merge, subset = nFeature_RNA > 200& percent.mt <= 20)
int.merge <- NormalizeData(int.merge, normalization.method = "LogNormalize", scale.factor = 10000)
int.merge <- FindVariableFeatures(int.merge,, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(int.merge)
int.merge <- ScaleData(int.merge, features = all.genes)
int.merge <- RunPCA(int.merge)
int.integrated <- IntegrateLayers(object = int.merge, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
int.integrated <- FindNeighbors(int.integrated, reduction = "harmony", dims = 1:30)
int.integrated <- FindClusters(int.integrated, resolution = 1)
int.integrated <- RunUMAP(int.integrated, reduction = "harmony", dims = 1:10)
int.integrated <- JoinLayers(int.integrated, assay = "RNA")
DimPlot(int.integrated, reduction = "umap",label=T)

int.integrated <- add_signatures_mouse(int.integrated,signatures_mouse,'LEVEL_1')
FeaturePlot(int.integrated,features=unique(signatures_mouse$LEVEL_1))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
VlnPlot(int.integrated,features=unique(signatures_mouse$LEVEL_1), group.by='seurat_clusters')
int.integrated <- add_signatures_mouse(int.integrated,signatures_mouse,'LEVEL_2')
FeaturePlot(int.integrated,features=unique(signatures_mouse$LEVEL_2))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
VlnPlot(int.integrated,features=unique(signatures_mouse$LEVEL_2), group.by='seurat_clusters')
int.integrated <- add_signatures_mouse(int.integrated,signatures_mouse,'LEVEL_3')
FeaturePlot(int.integrated,features=unique(signatures_mouse$LEVEL_3))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
VlnPlot(int.integrated,features=unique(signatures_mouse$LEVEL_3), group.by='seurat_clusters')

saveRDS(int.integrated, file = "int_integrated.rds")
int.integrated = readRDS(file = "int_integrated.rds")

#integration data
your_list <- list(tam1_s, tam2_s,myelo_s,lympho_s)
int.anchors <- FindIntegrationAnchors(object.list = your_list, dims = 1:10)
to_integrate <- Reduce(intersect, lapply(int.anchors@object.list, rownames))
int.combined <- IntegrateData(anchorset = int.anchors,features.to.integrate=to_integrate, dims = 1:10)
DefaultAssay(int.combined) <- "RNA"
int.combined <- JoinLayers(int.combined, assay = "RNA")
DefaultAssay(int.combined) <- "integrated"

int.combined <- ScaleData(int.combined, verbose = FALSE)
int.combined <- RunPCA(int.combined, npcs = 50, verbose = FALSE)
int.combined <- RunICA(int.combined, verbose = FALSE)
int.combined <- RunUMAP(int.combined, reduction = "pca", dims = 1:10)
int.combined <- FindNeighbors(int.combined, reduction = "pca", dims = 1:10)
int.combined <- FindClusters(int.combined, resolution = 1)
DimPlot(int.combined, reduction = "umap", group.by = "orig.ident")
DimPlot(int.combined, reduction = "umap",label=T)

saveRDS(int.combined, file = "int_combined.rds")
int.combined = readRDS(file = "int_combined.rds")

FeaturePlot(int.combined,features='nFeature_RNA')
VlnPlot(int.combined,features='nFeature_RNA', group.by='seurat_clusters')


FeaturePlot(int.combined,features=c('Epcam','Itgax','Cd74','Zbtb46','Cadm1','Xcr1','Cd79a','Ms4a1','Cd19','Ptprc'))
FeaturePlot(int.combined,features=c('Lyz2','Folr2','Spp1','Trem2'))
FeaturePlot(int.combined,features=c('Cd3e','Cd4','Cd8a','Cd8b'))
FeaturePlot(int.combined,features=c('Ptprc','Cd3e','Cd79a','C1qc','Folr2','Trem2','Nkg7','Ly6g','Itgax','Zbtb46','Clec9a','Wdfy4','C5ar1'))
FeaturePlot(int.combined,features=c('H2-Aa','H2-Eb1','H2-Ab1','Lyz2','Plac8','Cd14','Sirpa'))
FeaturePlot(int.combined,features=c('Lyz1','Lyz2'))
VlnPlot(int.combined,features=c('Lyz1','Lyz2'), group.by='seurat_clusters')


int.combined <- add_signatures_ortho(int.combined,signatures_human,'Family')
FeaturePlot(int.combined,features=unique(signatures_human$Family))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
VlnPlot(int.combined,features=c('Endothelial','Fibroblast','Immune','Lympho','Myeloid'), group.by='seurat_clusters')

int.combined <- add_signatures_mouse(int.combined,signatures_mouse,'LEVEL_1')
FeaturePlot(int.combined,features=unique(signatures_mouse$LEVEL_1))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
VlnPlot(int.combined,features=unique(signatures_mouse$LEVEL_1), group.by='seurat_clusters')
int.combined <- add_signatures_mouse(int.combined,signatures_mouse,'LEVEL_2')
FeaturePlot(int.combined,features=unique(signatures_mouse$LEVEL_2))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
VlnPlot(int.combined,features=unique(signatures_mouse$LEVEL_2), group.by='seurat_clusters')
FeaturePlot(int.combined,features=unique(signatures_mouse$LEVEL_3))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
VlnPlot(int.combined,features=unique(signatures_mouse$LEVEL_3), group.by='seurat_clusters')

find_best_resolution(int.combined,'LEVEL_1')

#DEG
int.combined.markers <- FindAllMarkers(int.combined,only.pos = TRUE)
write.csv(int.combined.markers,'Markers_integrated.csv')
int.combined.markers %>% 
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC >1) %>%
  slice_head(n=30) %>%
  ungroup() -> top30
png("top5_marker_int.png",width = 550, height =600,unit="px")
  DoHeatmap(int.combined,features =top5$gene) +NoLegend()
dev.off()

###ICA
int.combined <- RunICA(int.combined)


###Processed
load('processed/PYMT2Sel.Robj')
PYMT2Sel <- SeuratObject::UpdateSeuratObject(PYMT2Sel)
DimPlot(PYMT2Sel, reduction = "tsne")

### heift
dt = fread("../heift_2022_FOLR2_GSE192935/GSE192935_Raw_Counts_SC_Mouse_Fcgr1pos_FigS3.csv.gz")
rownames(dt) <- dt$V1
dt$V1 <- NULL
Fcgr1pos_heift <- CreateSeuratObject(  dt,  meta.data = NULL,  project = "Fcgr1pos_heift")

meta = data.frame(fread("../heift_2022_FOLR2_GSE192935/GSE192935_Metadata_SC_Mouse_Fcgr1pos_FigS3.csv.gz"))
meta$V1 = gsub( "-", ".",meta$V1)
rownames(meta) <- meta$V1
meta$V1 <- NULL
colnames(meta) <- paste0(colnames(meta),'_paper')
Fcgr1pos_heift <-AddMetaData(Fcgr1pos_heift, meta, col.name = NULL)

Fcgr1pos_heift_s = seurat_analyse(Fcgr1pos_heift)

dt = fread("../heift_2022_FOLR2_GSE192935/GSE192935_Raw_Counts_SC_Mouse_TAMs_Fig3.csv.gz")
rownames(dt) <- dt$V1
dt$V1 <- NULL
tam_heift <- CreateSeuratObject(  dt,  meta.data = NULL,  project = "TAM_heift")

meta = data.frame(fread("../heift_2022_FOLR2_GSE192935/GSE192935_Metadata_SC_Mouse_TAMs_Fig3.csv.gz"))
meta$V1 = gsub( "-", ".",meta$V1)
rownames(meta) <- meta$V1
meta$V1 <- NULL
colnames(meta) <- paste0(colnames(meta),'_paper')
tam_heift <-AddMetaData(tam_heift, meta, col.name = NULL)

tam_heift_s = seurat_analyse(tam_heift)
#integration
heift_int = seurat_integration(list(Fcgr1pos_heift_s,tam_heift_s))
DimPlot(heift_int, reduction = "umap", group.by = "orig.ident")
DimPlot(heift_int, reduction = "umap",label=T)

saveRDS(heift_int, file = "../heift_2022_FOLR2_GSE192935/heift_int.rds")

heift_int <- add_signatures_mouse(heift_int,signatures_mouse,'Major_cell_type')
FeaturePlot(heift_int,features=unique(signatures_mouse$Major_cell_type))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
VlnPlot(heift_int,features=unique(signatures_mouse$Major_cell_type), group.by='seurat_clusters')
heift_int <- add_signatures_mouse(heift_int,signatures_mouse,'Subset_or_state')
for(tissue in c('Mononuclear_phagocyte','Dendritic_cell','Mast_cell_or_basophil','NK_or_ILC')){
  subset = unique(signatures_mouse[signatures_mouse$Major_cell_type==tissue,]$Subset_or_state)
  ft_plot <- FeaturePlot(heift_int,features=subset)& 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  vln_plot <- VlnPlot(heift_int,features=subset, group.by='seurat_clusters')
  plot(ft_plot)
  plot(vln_plot)
}

