# make_h5ad.py GSE_dir

import anndata as ad
import scanpy as sc
import sys
import os
import re

ad.settings.allow_write_nullable_strings=True 
GSE_dir = sys.argv[1]
filenames = os.listdir(GSE_dir)
matrixfiles = [file for file in filenames if bool(re.search("matrix.mtx.gz",file))]
list_prefix = [re.findall(r".*(?=matrix.mtx.gz)", file)[0] for file in matrixfiles]
print(list_prefix)

for prefix in list_prefix:
    adata = sc.read_10x_mtx(path = GSE_dir, prefix=prefix)
    

    gsm_prefix = re.findall(r'GSM[0-9]*(?=_)', prefix)[0]
    adata.write_h5ad(os.path.join(GSE_dir, gsm_prefix + ".h5ad"))


