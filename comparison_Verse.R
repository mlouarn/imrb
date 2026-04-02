source('/home/marine-louarn/imrb/packages.R')
source('/home/marine-louarn/imrb/functions.R')

setwd("/home/marine-louarn/Documents/Alexandre_mouse/Done/")

myelo_int <- readRDS("Integrated_10_total_DC.rds")
ref_mreg <- readRDS('/home/marine-louarn/ref/20211018_MigrDC_Verse_Reference_Prep.rds')
ref_DC <- readRDS('/home/marine-louarn/ref/Reference_mDC_Verse.rds')

avPYMT <- AverageExpression(myelo_int,group.by ="RNA_snn_res.1")
avPYMT <- avPYMT$RNA[rowSums(avPYMT$RNA) != 0, ]

taxIDs <- setNames(c(9606,10090),
                   c("human","mouse"))
build_orthology_table(taxIDs = taxIDs,  primaryTaxID = 10090, 
                      outputFilePrefix="mouse_human_orthologs")

convert <- read.csv("mouse_human_orthologs_20260402.csv")
convert_by_symbol <- convert[!(is.na(convert$human_Symbol)|is.na(convert$mouse_Symbol)),c("human_Symbol","mouse_Symbol")] # Remove NAs from conversion table
convert_by_symbol <- convert_by_symbol[is.element(convert_by_symbol$mouse_Symbol,rownames(avPYMT)),] # Remove genes not in data matrix
dataOut <- avPYMT[match(convert_by_symbol$mouse_Symbol,rownames(avPYMT)),] # Subset data to include only genes with mouse othologs
rownames(dataOut) <- convert_by_symbol$human_Symbol # Convert gene names to mouse

avDC <- AverageExpression(ref_DC,group.by ="DC.Phenograph.Clusters")
avDC <- avDC$RNA[rowSums(avDC$RNA) != 0, ]
intersect_genes = rownames(avDC)[rownames(avDC) %in% rownames(dataOut)]

avDC = avDC[intersect_genes,]
dataOut=dataOut[intersect_genes,]

avDC<- CreateSeuratObject(avDC)
dataOut<- CreateSeuratObject(dataOut)

avDC_a <- NormalizeData(avDC)
avDC_a <- ScaleData(avDC_a,features = rownames(avDC))
avDC_scal = avDC_a@assays$RNA@layers$scale.data
colnames(avDC_scal) <- colnames(avDC_a)
dataOut_a <- NormalizeData(dataOut)
dataOut_a <- ScaleData(dataOut_a,features = rownames(dataOut))
dataOut_scal = dataOut_a@assays$RNA@layers$scale.data
colnames(dataOut_scal) <- as.character(seq(1,20))

cor <- cor(avDC_scal,dataOut_scal, use = 'pairwise.complete.obs', method = 'spearman')
ComplexHeatmap::pheatmap(cor, legend = TRUE, cluster_rows = T,
                         clustering_method = "ward.D", cluster_cols = T, heatmap_legend_param = list(title = "correlation"))
