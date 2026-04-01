import scanpy as sc
import sys

h5ad_files = sys.argv[1:-1]
outfile = sys.argv[-1]

if outfile.endswith('h5ad'):
    print('output file (last arg) cannot be h5ad')
    exit()

with open(outfile, 'w') as out:
    for h5ad in h5ad_files:
        print(h5ad)
        adata = sc.read_h5ad(h5ad)
        ncell = len(adata.obs_names)
        print(ncell)
        out.write(f'{h5ad}\t{ncell}\n')
