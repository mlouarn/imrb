
#list of signature genes

signatures_human = read.csv('/home/marinelouarn/Documents/Alexandre_mouse/20260204_Signatures_Cell_populations_HUMAN.csv')
signatures_mouse = read.csv('/home/marinelouarn/Documents/Alexandre_mouse/20260405_Mouse_cellPopulation_Signature_Genes_ConensusNEW.csv')
signatures_mouse = signatures_mouse[signatures_mouse$Keep_for_initial_screening=='y',]

add_signatures_ortho <- function(obj_seurat, signature_file, level_signature){
  for(tissue in as.list(unique(signature_file[level_signature]))[[1]]){
    genes = signature_file[signature_file[level_signature]==tissue,]$gene
    ortho = findOrthologsHsMm(from_filters = "hgnc_symbol",
                              from_values = genes ,
                              to_attributes = "external_gene_name")
    genes_mouse = ortho$external_gene_name
    print(tissue)
    obj_seurat <- AddModuleScore_UCell(obj_seurat,
                                       assay= 'RNA',
                                      features = list(genes_mouse),
                                      name = tissue)
    obj_seurat[[as.character(tissue)]] <- obj_seurat[[paste0('signature_1',tissue,'_human',)]]
  }
  gc()
  return(obj_seurat)
}

add_signatures_mouse <- function(obj_seurat, signature_file, level_signature){
  for(tissue in as.list(unique(signature_file[level_signature]))[[1]]){
    genes_mouse = signature_file[signature_file[level_signature]==tissue,]$gene
    print(tissue)
    obj_seurat <- AddModuleScore_UCell(obj_seurat,
                                       assay= 'RNA',
                                       features = list(genes_mouse),
                                       name = tissue)
    obj_seurat[[as.character(tissue)]] <- obj_seurat[[paste0('signature_1',tissue)]]
  }
  gc()
  return(obj_seurat)
}

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
