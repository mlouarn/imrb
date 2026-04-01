## subset_h5ad.py input_file output_file csv_seggeq

import anndata as ad
import scanpy as sc
import sys

ad.settings.allow_write_nullable_strings = True

LEN_HEADER_SG = 5 # the line where the barcodes/columns are
NB_COL_IGNORE = 1  # first col "gene"

h5ad_input = sys.argv[1]
h5ad_output = sys.argv[2]
csv_sg = sys.argv[3]

if h5ad_input==h5ad_output:
    print("input and output can't be the same")
    exit()

with open(csv_sg, 'r') as csv:
    cellids = csv.readlines()[LEN_HEADER_SG].split(',')[NB_COL_IGNORE:]

if cellids[0].isdigit():
    cellids = [int(cellids.strip()) for cellids in cellids] # last one has \n
else:
    cellids = [str(cellids.strip()) for cellids in cellids]
# try: 
#     cellids = [int(cellids.strip()) for cellids in cellids] # last one has \n
# except:
#     cellids = [cellids.strip() for cellids in cellids]

adata = sc.read_h5ad(h5ad_input)
print(f'subset {len(cellids)} cells from {len(adata.obs_names)}, first cell is {cellids[0]}')
adata[cellids,:].write_h5ad(h5ad_output)