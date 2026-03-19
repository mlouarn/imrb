suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))

setwd('/home/marine-louarn/Documents/PDAC/')

folders = list.dirs(path = "Datasets", full.names = TRUE, recursive = FALSE)
#folders = folders[-c(length(folders)-2,length(folders)-1,length(folders))]
folders <- gsub('./','',folders)

for (folder_GSE in folders) {
  print(folder_GSE)
  int.integrated = readRDS(paste0("Datasets/",folder_GSE,"/",folder_GSE,"_integrated.rds"))
  int.integrated[["CellName"]] <- colnames(int.integrated)
  
  cell_tokeep = suppressMessages(as.data.frame(read_csv(paste0("xtract_myelo/",folder_GSE,"_integrated_Myeloid.csv"), skip = 5)))
  cell_tokeep = colnames(cell_tokeep)
  cell_tokeep = cell_tokeep[-c(1,2)]
    nb_cell = length(colnames(cell.myelo))
  print(paste0('nb cells total: ',nb_cell))
  
  cell.myelo <- subset(cell.myelo, subset = CellName %in% cell_tokeep)
  nb_cell = length(colnames(cell.myelo))
  print(paste0('nb cells myeloid: ',nb_cell))
  saveRDS(cell.myelo, paste0("Datasets/",folder_GSE,"/",folder_GSE,"_myeloid.rds"))
}
