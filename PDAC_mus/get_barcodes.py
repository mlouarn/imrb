## create a txt for each GSM with the barcodes 
## (without the part that makes them unique in the merge)
## to then extract the mreg from each GSM separately

import scanpy as sc
import pandas as pd
import anndata as ad

CSV_HEADER = "data/csv_header.csv"

adpdac = sc.read_h5ad('data/pdac-mus_mreg.h5ad')
adpymt = sc.read_h5ad('data/pymt-rna_mreg.h5ad')
adata = ad.concat([adpdac, adpymt], label = 'pdac_pymt', keys=['pdac','pymt'], index_unique="-")
pdac_df = pd.read_csv('data/20260317 PDAC scRNAsseq datasets.xlsx - Mouse GSM Infos & Metadata.csv')
pymt_df = pd.read_csv('data/pyMT datasets to integrate - Feuille 2.csv')

gsm_to_remove = list(pdac_df.loc[pdac_df['GSM_In_mregDC-VERSE']=='no','GSM']) + list(pymt_df.loc[pymt_df['GSM_In_mregDC-VERSE']=='no','Sample ID'])

for gsm in set(adata.obs['orig.ident']):
    if (gsm in gsm_to_remove):
        continue
    print(gsm)
    barcodes = adata[adata.obs['orig.ident']==gsm].obs_names
    print(barcodes[0])
    pdac_or_pymt = adata[adata.obs['orig.ident']==gsm].obs['pdac_pymt'].iloc[0]
    with open('data/mreg_csv/'+gsm+'_'+pdac_or_pymt+'_mreg_barcodes.csv', 'w') as file:
        with open(CSV_HEADER, 'r') as header:
            file.write(header.read())
        file.write('genes,' + ','.join(barcodes) + '\n')
        file.write('genes,' + ','.join(barcodes) + '\n')
