FindPhenograph = function(object, k=30, cluster.name='pheno_clusters', reduction='pca'){
  mtx = object@reductions[[reduction]]@cell.embeddings
  pheno_out = Rphenograph(mtx, k=k)
  pheno_clusters = factor(membership(pheno_out[[2]]))
  names(pheno_clusters) = rownames(mtx)
  object@meta.data[cluster.name] = pheno_clusters
  return(object)
}