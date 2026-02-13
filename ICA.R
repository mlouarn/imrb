###ICA
int.combined = readRDS(file = "PYMT\ scRNAseq/boissonas_2022_tam-breast_GSE184096/int_combined.rds")

int.combined <- RunICA(int.combined)

genes = signatures_mouse$Signature_gene_mouse.style_symbol

df = as.data.frame(int.combined@reductions$ica@feature.loadings)
df$id = rownames(df)
df_s =df[df$id %in% genes, ]
df_s$id <- NULL
sum = colSums(df_s)
sum = rank(-sum)
sum = sum[sum<6]
list_ICA = rownames(as.matrix(sum))


df_ICA =select(as.data.frame(int.combined@reductions$ica@cell.embeddings),all_of(list_ICA))
df_pca = as.data.frame(int.combined@reductions$pca@cell.embeddings)

new_dim = merge(df_ICA, df_pca[,1:10], by=0)
rownames(new_dim) = new_dim$Row.names
new_dim$Row.names <- NULL
colnames(new_dim) <- paste0("ICAnPCA_", 1:15)

int.combined[["ICAnPCA"]] <- CreateDimReducObject(embeddings =  as.matrix(new_dim), key = "ICAnPCA_", assay = DefaultAssay(int.combined)) #custom dimred
DimPlot(int.combined, reduction = "ICAnPCA", pt.size = 0.5)

int.combined <- FindNeighbors(int.combined, reduction = "ICAnPCA")
int.combined <- FindClusters(int.combined, resolution = 1)
int.combined <- RunUMAP(int.combined,dims = 1:15)
DimPlot(int.combined, reduction = "umap", label = TRUE)

FeaturePlot(int.combined,features=unique(signatures_mouse$LEVEL_1))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
VlnPlot(int.combined,features=unique(signatures_mouse$LEVEL_1), group.by='seurat_clusters')

FeaturePlot(int.combined,features=unique(signatures_mouse$LEVEL_2))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
VlnPlot(int.combined,features=unique(signatures_mouse$LEVEL_2), group.by='seurat_clusters')

FeaturePlot(int.combined,features=unique(signatures_mouse$LEVEL_3))& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
VlnPlot(int.combined,features=unique(signatures_mouse$LEVEL_3), group.by='seurat_clusters')

#2nd option
df$Signature = with(df, ifelse(rownames(df) %in% genes, 'In_signature', 'Not'))
wilcox.test(IC_3~ Signature, 
            data = df,
            exact = FALSE)
calculate_Uscore(int.combined@reductions$ica@feature.loadings, genes)


