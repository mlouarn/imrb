library(AnnotationDbi)
library(dplyr)
library(org.Hs.eg.db)
library(stringr)
library(tibble)


raw_counts_ricky = read.table("ImmunoMSI_genecounts_FCcounts_s1_O.tsv")%>%
  mutate(gene = rownames(.))
raw_counts_nipicol = read.table("merged_raw_counts_Nipicol.tsv")%>%
  mutate(gene = rownames(.))
gene_length = read.table("gene_length.txt", header = T)

raw_counts = raw_counts_ricky%>%
  inner_join(raw_counts_nipicol, by = "gene")%>%
  inner_join(gene_length, by="gene")%>%
  column_to_rownames("gene")

rpk = raw_counts%>%
  mutate_at(str_subset(colnames(.),"RICKI|IMSI"), 
            function(x){x/raw_counts$mean}*10**3)%>%
  mutate(gene = raw_counts$gene)

tpm = rpk%>%
  mutate_at(str_subset(colnames(.),"RICKI|IMSI"), 
            function(x){x/sum(x)*10**6})%>%
  mutate(symbol = mapIds(org.Hs.eg.db, keys=rownames(.), column="SYMBOL", keytype="ENSEMBL"))



write.table(tpm, file = "bulkrna_tpm.tsv")
