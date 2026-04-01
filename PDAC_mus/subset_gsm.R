library(Seurat)
library(stringr)
library(dplyr)

# RDS_IN = "data/Integrated_10_total_DC.rds"
# RDS_OUT = "data/pymt_mreg.rds"
# GSM_VEC = c('GSM7159211','GSM4623140')

arguments = commandArgs(trailingOnly = TRUE)
RDS_IN = arguments[1]
RDS_OUT = arguments[2]
GSM_VEC = arguments[3:length(arguments)]

if(RDS_IN==RDS_OUT){
  stop("rds input file and rds output file cannot be the same")
}

print(GSM_VEC)
seu = readRDS(RDS_IN)
gsm_tokeep = unique(seu$orig.ident)[!unique(seu$orig.ident) %in% GSM_VEC]
print(gsm_tokeep)
print(seu)
seu = subset(seu, subset = (orig.ident %in% gsm_tokeep))
print(seu)
saveRDS(seu, RDS_OUT)
