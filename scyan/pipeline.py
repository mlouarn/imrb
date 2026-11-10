#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 16:53:44 2026

@author: marine-louarn
"""

import scyan
import pandas as pd
from pathlib import Path
import anndata
import os


project_name=Path("/home/marine-louarn/Documents/Data_for_Scyan/Laeticia_v3") #path where to save the project finishing by the project name (ex: "/documents/project1")
if not os.path.isdir(project_name):
    os.mkdir(project_name)
#if multiple file (only fcs)
multiple_fcs=True
path_multiple=Path("/home/marine-louarn/Documents/Data_for_Scyan/CD45+_Run1_ReExport") #path to directory containing all .fcs
#if single file
path_to_file =Path("/home/marine-louarn/Documents/Data_for_Scyan/concat_60_SEPTIMET_run1.fcs") #Add here path to cytometry file

batch_correction=True
is_cytof = False #True if CyTOF data, False is flow/spectralflow
path_to_marker_table =Path("/home/marine-louarn/Documents/Data_for_Scyan/Marker Table for Scyan_SIMPLE.v2.csv") #add here path to marker. Use .csv
umap_only_marker = True #for umap computation, True to use only cell-type marker, False for all markers

#TODO : add to the list of markers to the one in your fcs file and their replacement
markers='^cd|^hla|epcam|^ccr|^cx|^lox|mertk|^pdl' #add ^XXX to keep all strating with XXX, xxx for the specific xxx marker

def read_one_fcs(path):
    adata = scyan.read_fcs(path,marker_regex=markers)
    adata.obs["file"] = path.stem
    adata.obs["batch"] = path.stem # If you have batches, add here the batch of the corresponding path
    return adata


if multiple_fcs :
    folder_path = Path(path_multiple) # Replace "." by the path to your folder containing FCS files
    fcs_paths = [path for path in folder_path.iterdir() if path.suffix == ".fcs"]
    adata = anndata.concat([read_one_fcs(p) for p in fcs_paths], index_unique="-")
    #adata.obs['ID_cell']=adata.obs['FSC-A'].index
else:
    if path_to_file.suffix == ".csv":
        adata = scyan.read_csv(path_to_file,marker_regex=markers)
        #adata.obs['ID_cell']=adata.obs['FSC-A'].index
    elif path_to_file.suffix == ".fcs":
        adata = scyan.read_fcs(path_to_file,marker_regex=markers)
        #adata.obs['ID_cell']=adata.obs['FSC-A'].index


print(f"Created anndata object with {adata.n_obs} cells and {adata.n_vars} markers.\n\n-> The markers names are: {', '.join(adata.var_names)}\n-> The non-marker names are: {', '.join(adata.obs.columns)}")
original=adata.copy()

if is_cytof:  # we recommend asinh for CyTOF data
    scyan.preprocess.asinh_transform(adata)
else:  # we recommend auto_logicle for flow or spectral flow
    scyan.preprocess.auto_logicle_transform(adata)
scyan.preprocess.scale(adata)

print('preprocess OK')

table = pd.read_csv(path_to_marker_table, index_col=0)

nb_cells=len(adata.obs)
if umap_only_marker:
    scyan.tools.umap(adata, markers=table.columns, min_dist=0.1, n_cells=nb_cells)
else:
    scyan.tools.umap(adata, min_dist=0.1, n_cells=nb_cells)

scyan.plot.umap
print('UMAP OK')

#model
if batch_correction:
    model = scyan.Scyan(adata, table, batch_key="batch", prior_std=0.35) #change batch key if needed
else:
    model = scyan.Scyan(adata, table)#, prior_std=0.25, lr=0.0001)
model.fit()
model.predict()

print('Model OK')

#adata
original.obs['scyan_pop']=adata.obs['scyan_pop']
original.obsm['X_umap']=adata.obsm['X_umap']

scyan_pop_dic ={'B cells':1,
  'DC1':2,
  'DC2':3,
  'DC3 CD14+':4,
  'DC3 CD14-':5,
  'Eosinophils':6,
  'NK 56-':7,
  'NK 56bright':8,
  'NK 56dim':9,
  'PMN':10,
  'T cells':11,
  'cMo':12,
  'iMo':13,
  'ncMo':14,
  'pDC':15,
  'pre-DC':16}
original.obs['scyan_pop_ID']=adata.obs['scyan_pop'].replace(scyan_pop_dic)

original.obs['batch_ID']=original.obs['batch']
if 'batch' in original.obs.columns :
        original.obs['batch_ID'] = original.obs['batch_ID'].replace({"SEPTIMET_001_":"","_CD45\+":""}, regex = True)
        original.obs['batch_ID'] = original.obs['batch_ID'].astype(int)

scyan.write_fcs(original,os.path.join(project_name, 'scyan_fulloutput_concat.fcs'))
scyan.write_csv(original,os.path.join(project_name, 'scyan_fulloutput_concat.csv'))

for sample in original.obs['batch'].unique():
    adata_sample = original[original.obs['batch']==sample]
    scyan.write_fcs(adata_sample,os.path.join(project_name, 'scyanOUT_'+sample+'.fcs'))

##keep only specific cluster
#filtered_df = original[original.obs['scyan_pop_ID'] == 10]
#scyan.write_fcs(filtered_df,os.path.join(project_name, 'scyan_onlyPopulation10.fcs'))


'''
#population discovery

adata.obs['scyan_selected']='unselected'

#Subpopulation
scyan.tools.subcluster(adata, "PMN", resolution=0.1, markers=table.columns)

scyan.plot.umap(adata, color='scyan_subcluster_PMN')
scyan.plot.pops_expressions(model, key='scyan_subcluster_PMN')
'''
