library(Seurat,future)
library(dplyr)
library(data.table)
suppressMessages(library(tidyverse))
library(Azimuth)
source('/home/marinelouarn/imrb/packages.R')
source('/home/marinelouarn/imrb/functions.R')


#Le sample à baseline du patient N6 de NIRADO (N6_bl) n'aura pas de TCR associé dans l'objet. 
#Nous avons un souci avec les fastq TCR de ce sample en particulier. 

#Le patient N5 de NIRADO n'a pas de sample au 3eme timepoint qui correspond au combo, 
#mais il y a un sample à baseline "hors site" (N5_bl_Hs) : normalement c'est un prélèvement péri-tumoral. 

options(future.globals.maxSize = 8000 * 1024^2)
setwd(dir = "/home/marinelouarn/Documents/TCR_scRNA/")


Ari_P1 <- readRDS( "ARIANES_Seurat_obj_TCR/P1_processed_TCR.Rds")
DefaultAssay(Ari_P1) <- "RNA"
Ari_P1 <- JoinLayers(Ari_P1)
DimPlot(Ari_P1, reduction = "SCT_pca_Harmony_21_umap",  group.by = 'CELL_TYPE')

Nira_N1 <- readRDS( "NIRADO_Seurat_obj_TCR/N1_processed_TCR.Rds")
DimPlot(Nira_N1, reduction = "SCT_pca_Harmony_40_umap",  group.by = 'cell_names')

Nira_N6 <- readRDS( "NIRADO_Seurat_obj_TCR/N6_processed_TCR.Rds")
DimPlot(Nira_N6, reduction = "SCT_pca_Harmony_40_umap",  group.by = 'cell_names')





