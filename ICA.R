###ICA
int.combined <- RunICA(int.combined)

genes = signatures_mouse$Signature_gene_mouse.style_symbol

df = as.data.frame(int.combined@reductions$ica@feature.loadings)
df$id = rownames(df)
df_s =df[df$id %in% genes, ]
df_s$id <- NULL
sum = colSums(df_s)
sum = rank(-sum)
sum = sum[sum<6]
rownames(as.matrix(sum))

#2nd option
df$Signature = with(df, ifelse(rownames(df) %in% genes, 'In_signature', 'Not'))
wilcox.test(IC_3~ Signature, 
            data = df,
            exact = FALSE)
calculate_Uscore(int.combined@reductions$ica@feature.loadings, genes)


