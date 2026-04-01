# subset_rds.R rdsfile_in rdsfile_out SG_file


library(dplyr)
library(stringr)
library(Seurat)

arg = commandArgs(trailingOnly = TRUE)
rds_in = arg[1]
rds_out = arg[2]
SG_file = arg[3]

if(rds_in==rds_out){
  stop("rds input file and rds output file cannot be the same")
}

extract_cellids_from_SG = function(SG_file){
  SG_df = read.csv(SG_file, skip = 5, header = T)
  cellids = colnames(SG_df)[2:length(colnames(SG_df))]
  # SG replaces - with . we have to change it back
  cellids = str_replace_all(cellids, '\\.', '-')
  return(cellids)
}

seurat_obj = readRDS(rds_in)
cellids = extract_cellids_from_SG(SG_file)
seurat_obj_sub = subset(seurat_obj, cells = cellids)
saveRDS(seurat_obj_sub, rds_out)