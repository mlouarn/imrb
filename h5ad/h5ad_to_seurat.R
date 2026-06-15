i=basename

df <- as.data.frame(fread(paste0(i,"_matrix.tsv")))
genes <- fread(paste0(i,"_features.tsv"),header=F)
barcodes <- fread(paste0(i,"_barcodes.tsv"),header=F)
rownames(df) <- genes$V1
colnames(df) <- barcodes$V1
obj <- CreateSeuratObject(counts = df, project = i, min.cells = 3)
#add umap
mat <- fread(paste0(i,"_reduction.tsv"))
mat$V1 <- NULL
mat = as.matrix(mat)
rownames(mat)<- barcodes$V1
obj[['umap_python']] <- CreateDimReducObject(embeddings = mat, key = 'UMAPpython_', assay = 'RNA')
#add spatial
spatial <- fread(paste0(i,"_spatial.tsv"))
spatial$V1 <- NULL
spatial = as.matrix(spatial)
rownames(spatial)<- barcodes$V1
obj <- AddMetaData(obj, metadata = spatial, col.name = 'X')
obj <- AddMetaData(obj, metadata = spatial, col.name = 'Y')
#meta
metadata_df = as.data.frame(fread(paste0(i,"_metadata.tsv")))
rownames(metadata_df)= metadata_df$cell_id
obj <- AddMetaData(obj, metadata = metadata_df, col.name = 'leiden_ICPC')
obj <- seurat_analyse(obj)
saveRDS(obj, file = paste0(i,".rds"))
