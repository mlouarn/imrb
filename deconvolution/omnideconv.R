library(omnideconv)
library(Seurat)

sc_mat = read.table("data-sc/GSE132257_GEO_processed_protocol_and_fresh_frozen_raw_UMI_count_matrix.txt.gz",
                     header = T, row.names = 1) # SMC01
metadata = read.table("data-sc/GSE132257_processed_protocol_and_fresh_frozen_cell_annotation.txt",
                      header = T, sep = "\t")
cell_type_annotations = metadata$Cell_type
bk_mat = read.table("bulkrna_counts.tsv")
batch_ids = colnames(bk_mat)

## Pour savoir quel format mettre (brut ou TPM)
## et s'il faut initialiser avec build_model :
## https://github.com/omnideconv/omnideconv?tab=readme-ov-file#signature-matrixmodel-building

deconvolution <- omnideconv::deconvolute(bulk_gene_expression = bk_mat,
                                         method="CPM",
                                         single_cell_object=sc_mat, cell_type_annotations, batch_ids)

saveRDS(deconvolution, "deconvolution.rds")
write.table(deconvolution, "deconvolution.txt", sep = "\t")

# signature = omnideconv::build_model(single_cell_object = sc_mat,
#                                     cell_type_annotations = cell_type_annotations,
#                                     method = "DWLS",
#                                     batch_ids = batch_ids,
#                                     bulk_gene_expression = bk_mat)
# saveRDS(signature, "signature.rds")

# deconvolution <- deconvolute(bulk_gene_expression = bk_mat,
#                              signature = signature,
#                              method='DWLS')
