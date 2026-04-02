library(Seurat,future)
library(dplyr)
suppressMessages(library(tidyverse))
#plan(workers = 2)
options(future.globals.maxSize = 8000 * 1024^2)

setwd("/home/marine-louarn/Documents/Alexandre_mouse/Done/")

myelo_2 = readRDS("Integrated_myelo_2_total_DC.rds")
myelo_1 =readRDS("Integrated_myelo_1_total_DC.rds")
myelo_3 =readRDS("GSE131957/GSE131957_DC.rds")
DefaultAssay(myelo_2) <- 'RNA'
DefaultAssay(myelo_1) <- 'RNA'

myelo_2 <- JoinLayers(myelo_2)
myelo_2 <- SCTransform(myelo_2)

myelo_1 <- JoinLayers(myelo_1)
myelo_1 <- SCTransform(myelo_1)

myelo_DC <- merge(myelo_1, list(myelo_2,myelo_3), project = "integrated_total")

frequencies <- table(myelo_DC$orig.ident)
sample_to_keep = names(frequencies)[frequencies > 30] 
sample_to_keep = sample_to_keep[sample_to_keep %in% c('GSE264232_NTC-1','GSE264232_NTC-2','GSE264232_NTC-3','GSE264232_IL1a-2','GSE264232_IL1a-3','GSM8208011','GSM8208012','GSM8208013','GSM8208014','GSM8208015','GSM8208016','GSM7159211','GSM7159213','GSM7159215','GSM7159217','GSM9001991','GSM9001992','GSM9001993','GSM9001994','GSM8993905','GSM8993906','GSM8993907','GSM8993908','GSM4623139','GSM4623140','GSM3978673','GSM3978675','GSE214518_TAMs','GSE214518_TIMs','GSM8212431','GSM8212432','GSM8830813_','GSM9217269','GSM9217270','GSM9217271','GSM9217272','GSM9217273','GSM9217274','GSM9217275','GSM9217276','GSM9217277','GSM9217278','GSM9217279','GSM9217280','GSM9217281','GSM9217282','GSM9217283','GSM9217284','GSM9217285','GSM9217286','GSM9217287','GSM9217288','GSM9217289','GSM9217290','GSM9217291','GSM9217292','GSM9217293','GSM9217294','GSM9217295','GSM9217296','GSM9217297','GSM9217298','GSM3832736','GSM3832738','GSM3832740','GSM3832742')]

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
myelo_int[["RNA"]] <- JoinLayers(myelo_int[["RNA"]])
myelo_int <- RunICA(myelo_int,nics=100)
myelo_int <- RunPCA(myelo_int,npcs=100)
myelo_int <- FindNeighbors(myelo_int, reduction = "integrated.cca", dims = 1:50)
for (i in seq(0.2,1.2,0.2)){
  myelo_int <- FindClusters(myelo_int, resolution = i)
}

#myelo_int <- RunUMAP(myelo_int, dims = c(1,3,5,6,7,9,10,12,16,19,22,32,35,40,44,47,48), reduction = "integrated.cca",min.dist = 0.05)
#myelo_int <- RunUMAP(myelo_int, dims = c(1:10,12,16,19,22,32,35,40,44,47,48), reduction = "integrated.cca",min.dist = 0.05)
myelo_int <- RunUMAP(myelo_int, dims = 1:50, reduction = "integrated.cca",min.dist = 0.05)
DimPlot(myelo_10, reduction = "umap", group.by = c("RNA_snn_res.1.2"),label=TRUE,cols=DiscretePalette(28,palette='polychrome'))

myelo_int <- FindNeighbors(myelo_int, reduction = "integrated.cca", dims = 1:12)
for (i in seq(0.2,1.2,0.2)){
  myelo_int <- FindClusters(myelo_int, resolution = i)
}
myelo_int <- RunUMAP(myelo_int, dims = 1:12, reduction = "integrated.cca",min.dist = 0.3,reduction.key = 'umap03_',reduction.name='umap03')
myelo_int <- RunUMAP(myelo_int, dims = 1:12, reduction = "integrated.cca",min.dist = 0.05,reduction.key = 'umap005_',reduction.name='umap005')

DimPlot(myelo_int, reduction = "umap005", group.by = c("RNA_snn_res.1"),label=TRUE,cols=DiscretePalette(28,palette='polychrome'))
DimPlot(myelo_int, reduction = "umap", group.by = c("orig.ident"))

myelo_int[["percent.rp"]] <- PercentageFeatureSet(myelo_int, pattern ="^Rp[ls]")
FeaturePlot(myelo_int, features = c('percent.rp'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(myelo_int, features = c('nFeature_RNA'))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

myelo_int.markers <- FindAllMarkers(myelo_int,group.by ="RNA_snn_res.1", only.pos = TRUE)
write_csv(myelo_int.markers, "Total_DC_cluster_res1_Markers.csv")
markers <- FindAllMarkers(myelo_int,group.by ="RNA_snn_res.1", only.pos = TRUE,logfc.threshold = 1)
markers <- markers[markers$p_val_adj<0.05,]
write_csv(markers, "Res_1_PYMT_DEG_pval005.csv")
markers001 <- markers[markers$p_val_adj<0.01,]
write_csv(markers001, "Res_1_PYMT_DEG_pval001.csv")
markers_1 <- markers[markers$p_val_adj<0.001,]
#markers_1 <- markers_1[markers_1$avg_log2FC>1,]
write_csv(markers_1, "Res_1_PYMT_DEG_pval0.001_logFC1.csv")

mono.de <- FindMarkers(myelo_int, ident.1 = "10", ident.2 = "12",group.by='Cluster_res1.2', verbose = FALSE)
write.csv(mono.de, "Res1.2_Markers_clust10vs12.csv")

#here
DEG_forheat=as.data.frame(read_excel('Res_1_PYMT_DEG_pval0.001_logFC1.xlsx', 1))
DEG_forheat = DEG_forheat[DEG_forheat$KEEP %in% c('x'),]
myelo_int$Cluster_res1_order <- factor(
  myelo_int$Cluster_res1,
  levels = c(4,1,11,17,20,16,12,19,14,8,15,18,6,10,2,5,3,7,9,13)
)
DoHeatmap(myelo_int, features= DEG_forheat$gene,group.by = 'Cluster_res1_order')+ scale_fill_gradientn(colors = c("blue", "white", "red"))
avPYMT <- AverageExpression(myelo_int,group.by ="RNA_snn_res.1")
avPYMT <- avPYMT$RNA[rowSums(avPYMT$RNA) != 0, ] 
avPYMT_a <- NormalizeData(avPYMT)
avPYMT_a <- ScaleData(avPYMT_a,features = rownames(avPYMT))
#avPYMT_scal = avPYMT_a@assays$RNA@layers$scale.data
colnames(avPYMT_a) <- as.character(seq(1,20))
avPYMT_a=avPYMT_a[DEG_forheat$gene,]
Heatmap(avPYMT_a, name = "mat",cluster_rows =FALSE, show_column_dend = FALSE,column_names_side ='top',row_names_side = c("left"), show_row_dend = FALSE,column_order =c(4,1,11,17,20,16,12,19,14,8,15,18,6,10,2,5,3,7,9,13))

dat_z = myelo_int.markers %>%
  group_by(cluter) %>%
  mutate(z_score = scale(y))
d
clustree(myelo_int, prefix = "RNA_snn_res.")

saveRDS(myelo_int, "Integrated_10_total_DC.rds")
myelo_int <- readRDS("Integrated_10_total_DC.rds")

SG_Header <- read.table("../Header_SeqGeq.txt",sep="\t",header=F,row.names=1,check.names = F)
SG_Header <- unname(SG_Header)

new_genes <-read.csv('New_genes_of interest.csv')
new_genes <- unique(new_genes$genes)

counts_matrix <- as.data.frame(myelo_int[["RNA"]]$data) %>%
  filter(row.names(.) %in% genes_signature)

counts_matrix_new <- as.data.frame(myelo_int[["RNA"]]$data) %>%
  filter(row.names(.) %in% unique(c(new_genes,genes_signature)))

counts_matrix_full <- as.data.frame(myelo_int[["RNA"]]$data) 
counts_matrix_full <-counts_matrix_full[rowSums(counts_matrix_full) != 0, ]

umap03 = t(myelo_int[["umap03"]]@cell.embeddings)
umap005 = t(myelo_int[["umap005"]]@cell.embeddings)
sample=myelo_int$orig.ident
dataset=myelo_int$orig.ident
samp_data=myelo_int$orig.ident
source('seqgeq_ref.R')
cca=t(myelo_int@reductions$integrated.cca@cell.embeddings)
cluster02 = myelo_int$RNA_snn_res.0.2
cluster04 = myelo_int$RNA_snn_res.0.4
cluster06 = myelo_int$RNA_snn_res.0.6
cluster08 = myelo_int$RNA_snn_res.0.8
cluster1 = myelo_int$RNA_snn_res.1
cluster1.2 = myelo_int$RNA_snn_res.1.2
PYMT=myelo_int$orig.ident
PYMT[!(PYMT %in% c("GSM3832736","GSM3832738","GSM3832740","GSM3832742"))]='1'
PYMT[PYMT %in% c("GSM3832736","GSM3832738","GSM3832740","GSM3832742")]='0'

total = rbind(counts_matrix_full,umap03,umap005,cca,sample,samp_data,dataset,cluster02,cluster04,cluster06,cluster08,cluster1,cluster1.2,PYMT,new_umap)
rownames(total)[rownames(total) == "78674"] <- "Sample"
rownames(total)[rownames(total) == "78675"] <- "Sample_in_Dataset"
rownames(total)[rownames(total) == "78676"] <- "Dataset"
rownames(total)[rownames(total) == "78677"] <- "Clusters.0.2"
rownames(total)[rownames(total) == "78678"] <- "Clusters.0.4"
rownames(total)[rownames(total) == "78679"] <- "Clusters.0.6"
rownames(total)[rownames(total) == "78680"] <- "Clusters.0.8"
rownames(total)[rownames(total) == "78681"] <- "Clusters.1"
rownames(total)[rownames(total) == "78682"] <- "Clusters.1.2"
rownames(total)[rownames(total) == "78683"] <- "is_PYMT"
write.table(SG_Header, file = "Integrated_10_total_DC_full_newUMAP.txt",sep="\t")
suppressWarnings(write.table(total, file = "Integrated_10_total_DC_full_newUMAP.txt",sep="\t",append=T))

total = rbind(counts_matrix_new,umap03,umap005,cca,sample,samp_data,dataset,cluster02,cluster04,cluster06,cluster08,cluster1,cluster1.2,PYMT,new_umap)
rownames(total)[rownames(total) == "292"] <- "Sample"
rownames(total)[rownames(total) == "293"] <- "Sample_in_Dataset"
rownames(total)[rownames(total) == "294"] <- "Dataset"
rownames(total)[rownames(total) == "295"] <- "Clusters.0.2"
rownames(total)[rownames(total) == "296"] <- "Clusters.0.4"
rownames(total)[rownames(total) == "297"] <- "Clusters.0.6"
rownames(total)[rownames(total) == "298"] <- "Clusters.0.8"
rownames(total)[rownames(total) == "299"] <- "Clusters.1"
rownames(total)[rownames(total) == "300"] <- "Clusters.1.2"
rownames(total)[rownames(total) == "301"] <- "is_PYMT"

write.table(SG_Header, file = "Integrated_10_total_DC_newUMAP.txt",sep="\t")
suppressWarnings(write.table(total, file = "Integrated_10_total_DC_newUMAP.txt",sep="\t",append=T))


###Get old umap
old_umap_coordinates_1 = suppressMessages(as.data.frame(read_delim("Leader_Lung.csv",delim=',', skip = 5)))
#old_umap_coordinates_1 = old_umap_coordinates_1[old_umap_coordinates_1[,1] %in% c('umap_1','umap_2'), ]
names <- colnames(old_umap_coordinates_1)
rownames(old_umap_coordinates_1) <- paste0('LUNG_',old_umap_coordinates_1[,1])
old_umap_coordinates_1[,1]<- NULL
colnames(old_umap_coordinates_1)<-c('wt_naive_AAACCTGAGGGATCTG',names[2:(length(names)-1)])

old_umap_coordinates_2 = suppressMessages(as.data.frame(read_delim("PyMT.csv",delim=',', skip = 5)))
#old_umap_coordinates_2 = old_umap_coordinates_2[old_umap_coordinates_2[,1] %in% c('umap_1','umap_2'), ]
names <- colnames(old_umap_coordinates_2)
rownames(old_umap_coordinates_2) <- paste0('PYMT_',old_umap_coordinates_2[,1])
old_umap_coordinates_2[,1]<- NULL
colnames(old_umap_coordinates_2)<-c('AAAGATGGTACTCGCG-1_1_1_1_1_1_1_1',names[3:(length(names))])

new_umap = as.data.frame(umap03)

for (i in colnames(new_umap)){
  if (i %in% colnames(old_umap_coordinates_1)){
    new_umap['LUNG_umap_1',i]=old_umap_coordinates_1['LUNG_umap_1',i]
    new_umap['LUNG_umap_2',i]=old_umap_coordinates_1['LUNG_umap_2',i]
    new_umap['LUNG_Cluster',i]=old_umap_coordinates_1['LUNG_seurat_clusters',i]
  }
  if (i %in% colnames(old_umap_coordinates_2)){
    new_umap['PYMT_umap_1',i]=old_umap_coordinates_2['PYMT_umap_1',i]
    new_umap['PYMT_umap_2',i]=old_umap_coordinates_2['PYMT_umap_2',i]
    new_umap['PYMT_Cluster',i]=old_umap_coordinates_2['PYMT_Clusters.1',i]
  }
}
new_umap <- new_umap[-c(1:2), ]

total = rbind(counts_matrix_full,umap03,umap005,cca,sample,samp_data,dataset,cluster02,cluster04,cluster06,cluster08,cluster1,cluster1.2,PYMT,new_umap)
write.table(SG_Header, file = "Integrated_10_total_DC_full_newUMAP.txt",sep="\t")
suppressWarnings(write.table(total, file = "Integrated_10_total_DC_full_newUMAP.txt",sep="\t",append=T))


#Dendrogam
DefaultAssay(myelo_int)<-'RNA'
cluster.averages <- AverageExpression(myelo_int,group.by ="RNA_snn_res.1")
cluster.averages<- CreateSeuratObject(cluster.averages$RNA)
myelo_int.markers <- FindAllMarkers(myelo_int,group.by ="RNA_snn_res.1.2", only.pos = TRUE)

Test <- NormalizeData(cluster.averages)
Test <- ScaleData(Test,features = rownames(myelo_int))
#Test <- ScaleData(Test,features = myelo_int.markers$gene)

clustertest = Test@assays$RNA@layers$scale.data

cor <- cor(clustertest, use = 'pairwise.complete.obs', method = 'spearman')
colnames(cor)<-as.character(seq(1,20))
rownames(cor)<-as.character(seq(1,20))

callback = function(hc, mat){
  sv = svd(t(mat))$v[,c(2)]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

ComplexHeatmap::pheatmap(cor, legend = TRUE, cluster_rows = T,
                         clustering_method = "ward.D", cluster_cols = T, heatmap_legend_param = list(title = "correlation"))
pheatmap(cor, legend = TRUE, cluster_rows = T,
         clustering_method = "ward.D", cluster_cols = T, heatmap_legend_param = list(title = "correlation"), clustering_callback=callback)

###Metacluster
myelo_int$MetaCluster =myelo_int$Cluster_res1
myelo_int$MetaCluster[myelo_int$MetaCluster %in% c(4,1,11)]="DC1"
myelo_int$MetaCluster[myelo_int$MetaCluster %in% c(17,20,16,12,19)]="Prolif"
myelo_int$MetaCluster[myelo_int$MetaCluster %in% c(14,8,15)]="mregDC"
myelo_int$MetaCluster[myelo_int$MetaCluster %in% c(18)]="ISG_DC"
myelo_int$MetaCluster[myelo_int$MetaCluster %in% c(6,10,2)]="DC3"
myelo_int$MetaCluster[myelo_int$MetaCluster %in% c(5)]="SiglecH_DC"
myelo_int$MetaCluster[myelo_int$MetaCluster %in% c(3,7)]="DC2"
myelo_int$MetaCluster[myelo_int$MetaCluster %in% c(9)]="Tdblets"
myelo_int$MetaCluster[myelo_int$MetaCluster %in% c(13)]="Mt_genes"

MetaCluster.markers <- FindAllMarkers(myelo_int,group.by ="MetaCluster", only.pos = TRUE,logfc.threshold = 1)
MetaCluster.markers <- MetaCluster.markers[MetaCluster.markers$p_val_adj<0.001,]
write_csv(MetaCluster.markers, "MetaCluster_PYMT_DEG_pval0.001_logFC1.csv")

saveRDS(myelo_int, "Integrated_10_total_DC.rds")

