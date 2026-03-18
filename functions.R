
#list of signature genes

signatures_mouse = read.csv("~/Documents/Alexandre_mouse/Done/signature.v1.csv")
signatures_mouse = signatures_mouse[signatures_mouse$Keep_for_initial_screening=='y',]
signatures_mouse_b = read.csv("~/Documents/Alexandre_mouse/signature_boissonnas.v1.csv")
genes_signature = unique(c(signatures_mouse$Signature_gene_mouse.style_symbol,signatures_mouse_b$Signature_genes))

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
    if (tissue!=""){
      genes_mouse = signature_file[signature_file[level_signature]==tissue,]$Signature_gene_mouse.style_symbol
      print(tissue)
      obj_seurat <- AddModuleScore_UCell(obj_seurat,
                                         assay= 'RNA',
                                         features = list(genes_mouse),
                                         name = tissue)
      obj_seurat[[as.character(tissue)]] <- obj_seurat[[paste0('signature_1',tissue)]]
    }
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

seurat_integration_sct<- function(merged_seurat){
  merged_seurat <- SCTransform(merged_seurat, vars.to.regress = "percent.mt", verbose = FALSE)
}

find_best_resolution<- function(seurat_obj,level_signature){
  max=0
  best=0
  for(i in seq(0.2,1,0.1)){
    seurat_obj <- FindClusters(seurat_obj, resolution = i)
    states = unique(signatures_mouse[level_signature])[[1]]
    states = states[states != ""]
    clusters = table(seurat_obj[[paste0('integrated_snn_res.',toString(i))]])
    df <- data.frame(matrix(ncol = length(states), nrow = length(clusters)))
    colnames(df)=states
    rownames(df)=as.character(c(0:(length(clusters)-1)))
    print(i)
    for(tissue in as.list(states)){
      a <- DotPlot(object = seurat_obj, features = tissue, group.by = paste0('integrated_snn_res.',toString(i)))
      average = a$data$avg.exp.scaled
      df[tissue]= average
    }
    tmp =sum(sapply(df, function(x) sum(x > 2))> 0 ) +#check number of subset with at least 1 cluster expr >2
      sum(sapply(as.data.frame(t(df)), function(x) sum(x > 2))==1)#check number of cluster with 1 subset
    if (tmp>max){
      max=tmp
      best = paste0('integrated_snn_res.',toString(i))
    }
  }
  return(best)
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
