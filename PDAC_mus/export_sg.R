# rds_file list_gene_file file_out_seqgeq header_sg
library(Seurat)
library(dplyr, include.only = c("%>%", "filter"))

args = commandArgs(trailingOnly = TRUE)
rds_file = args[1]
file_out = args[2]
list_gene_file = args[3]
header_sg = args[4]

if(rds_file==file_out){
  stop("input file and output file cannot be the same")
}

print(args)
seurat_obj = readRDS(rds_file)
list_gene = read.table(list_gene_file, header = T)%>%
  unlist(use.names = F)

seqgeq_file <- function(seurat_obj, list_genes, file_out){
  SG_Header <- read.table(header_sg,sep="\t",header=F,row.names=1,check.names = F)
  SG_Header <- unname(SG_Header)
  counts_matrix <- as.data.frame(seurat_obj[["RNA"]]$data) %>%
    filter(row.names(.) %in% list_genes)
  umap = t(seurat_obj[["umap"]]@cell.embeddings)
  total = rbind(counts_matrix,umap)
  
  write.table(SG_Header, file = file_out,sep="\t")
  suppressWarnings(write.table(total, file = file_out,sep="\t",append=T))
}

seqgeq_file(seurat_obj, list_gene, file_out)
