# run in the dir with the GSM
library(Seurat)
library(stringr)

for(GSMrds in list.files('.',pattern = '^GSM.*\\.rds')){
  GSMid = str_extract(GSMrds, '^GSM[0-9]{5,9}')
  print(GSMid)
  seu = readRDS(GSMrds)
  seu$orig.ident = GSMid
  saveRDS(seu, GSMrds)
}