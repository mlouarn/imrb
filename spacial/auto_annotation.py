adata=sc.read_h5ad('/home/marine-louarn/Documents/Spacial_Roussy/documents_20260410/SFTP_SHARING_ARC_PGA/region_20EN-1063.h5ad')

#conda activate annotationPy
import celltypist
from celltypist import models
import scanpy as sc

sample1= sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260528_Sample1_novae.h5ad")
liver = sc.read_h5ad("/home/marine-louarn/ref/MoMac_n_DC_Liver_Verse.h5ad")
hcc = sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/refs/HCC_ref.h5ad")
#sample1=sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1.h5ad")
hcc_dutertre=sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/refs/GSE156625_HCCscanpyobj.h5ad")
#gbm_ref=sc.read_h5ad('/home/marine-louarn/Documents/Spacial_Roussy/refs/GSE328755_HK177_02_processed.h5ad')

sc.pp.highly_variable_genes(hcc_dutertre, flavor="seurat", n_top_genes=2000)
new_model = celltypist.train(hcc_dutertre[:, hcc_dutertre.var.highly_variable], labels = hcc_dutertre.obs['louvain'], check_expression = False)
new_model.write('/home/marine-louarn/Documents/Xenium_Calderaro/refs/model_HCC.pkl')

sample1= sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260528_Sample1_novae.h5ad")
sample1.X=sample1.layers['counts']
sc.pp.normalize_total(sample1, target_sum=1e4)
sc.pp.log1p(sample1)

predictions = celltypist.annotate(sample1, model = '/home/marine-louarn/Documents/Xenium_Calderaro/refs/model_HCC.pkl')
predictions.to_table(folder = '/home/marine-louarn/Documents/Xenium_Calderaro/', prefix = '')
adata = predictions.to_adata(insert_labels = True, insert_conf = True)
adata.obs=adata.obs.rename(columns={"predicted_labels": "Annotation_CellTypist"})
adata.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260602_Sample1_SimpleRef.h5ad")

#tangram
adata_simple= sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260602_Sample1_SimpleRef.h5ad")
sdata_v2 = sd.read_zarr("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1.zarr")
sdata_v2["table"] = adata_simple
adata_simple.X = adata_simple.layers["counts"] 

sopa.utils.tangram_annotate(sdata_v2, hcc, "Cell_Type")
adata_simple.obs=adata_simple.obs.rename(columns={"Cell_Type": "Annotation_Tangram"})
adata_scDutertre.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260602_Sample1_tangram.h5ad")


#conda activate annotationPy
import singler
import scanpy as sc
adata_scDutertre= sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260602_Sample1_SimpleRef.h5ad")
hcc_dutertre=sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/refs/GSE156625_HCCscanpyobj.h5ad")

features_obj = adata_scDutertre.var.index
count_matrix_obj = adata_scDutertre.layers['counts'].T
features_ref = hcc.var.index
count_matrix_ref = hcc.layers['counts'].T
label_ref = hcc.obs['Cell_Type']
results = singler.annotate_single(
    test_data = count_matrix_obj,
    test_features = features_obj,
    ref_data = count_matrix_ref,
    ref_features=features_ref,
    ref_labels = label_ref,
)
adata_scDutertre['Annotation_SingleR']=results['best']
adata_scDutertre.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260602_Sample1_singleR.h5ad")

#conda activate scconcept
from concept import scConcept
import scanpy as sc
concept = scConcept(cache_dir='tmp')
concept.load_config_and_model(model_name='corpus40M-model30M') 
result = concept.extract_embeddings(adata=adata, gene_id_column='gene_id')
adata.obsm['X_scConcept'] = result['cls_cell_emb']
adata.write_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1_scconcept.h5ad")


#conda activate spacialPy
import pyucell as uc
from collections import defaultdict
import scanpy as sc
import pandas as pd
signature = pd.read_csv("/home/marine-louarn/Documents/test/20260204_Signatures_Cell_populations_HUMAN.csv")
sample1=sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260528_Sample1_ICPC.h5ad")

filtered_signature = signature[signature['gene'].isin(sample1.var_names.tolist())]
marker_gene_dictionary = defaultdict(list)
for idx, row in filtered_signature.iterrows():
    marker_gene_dictionary[row['Major_Cell_Populations']].append(row['gene'])

marker_gene_dictionary=dict(marker_gene_dictionary)
uc.compute_ucell_scores(sample1, signatures=marker_gene_dictionary, chunk_size=500)
sc.pl.umap(sample1, color="Macro_UCell", cmap="viridis", size=20)
sc.pl.embedding(sample1, color="Macro_UCell", cmap="viridis", size=20,basis='umapICPC')