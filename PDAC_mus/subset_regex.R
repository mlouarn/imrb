library(Seurat)
library(stringr)
library(dplyr)

# RDS_IN = "data/GSM5914681.rds"
# RDS_OUT = "data/GSM5914681_mreg.rds"
# BARCODE_PATH = "data/GSM5914681_mreg_barcodes.csv"

arguments = commandArgs(trailingOnly = TRUE)
RDS_IN = arguments[1]
RDS_OUT = arguments[2]
BARCODE_PATH = arguments[3]

if(RDS_IN==RDS_OUT){
  stop("rds input file and rds output file cannot be the same")
}

keep_selected_barcodes = function(barcode, list_selected_barcodes){
  detection_count = sum(str_detect(list_selected_barcodes, barcode))
  if(detection_count>1){
    print('===============')
    print(barcode)
    print(detection_count)
    print('===============')
  }
  if(detection_count==1){
    # print(barcode)
    # print(str_subset(list_selected_barcodes, barcode))
    return(barcode)
  }
  else
    return('NA')
}

extract_cellids_from_SG = function(SG_file){
  SG_df = read.csv(SG_file, skip = 5, header = T)
  cellids = colnames(SG_df)[2:length(colnames(SG_df))]
  # SG replaces - with . we have to change it back
  cellids = str_replace_all(cellids, '\\.', '-')
  return(cellids)
}

barcodes_mreg = extract_cellids_from_SG(BARCODE_PATH)
seu = readRDS(RDS_IN)

barcodes_kept = lapply(colnames(seu), keep_selected_barcodes, 
                       barcodes_mreg)%>%
  unlist()
barcodes_kept = barcodes_kept[barcodes_kept!='NA']
print(barcodes_kept[1])
print(barcodes_mreg[1])
print(seu)
seu = subset(seu, cells = barcodes_kept)
print(seu)
saveRDS(seu, RDS_OUT)
