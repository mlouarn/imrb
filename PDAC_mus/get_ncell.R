library(dplyr)
library(Seurat)

arguments = commandArgs(trailingOnly = TRUE)
rds_files = arguments[1:length(arguments)-1]%>%unlist()
outfile = arguments[length(arguments)]

get_ncell = function(rdsfile){
  print(rdsfile)
  seu = readRDS(rdsfile)
  print(seu)
  print(length(rownames(seu)))
  return(length(rownames(seu)))
}

print(rds_files)
ncell = lapply(rds_files, get_ncell)%>%unlist(use.names = F)
print(ncell)

df = data.frame("rds_file"=rds_files, "ncell"=ncell)
write.table(df, outfile)
