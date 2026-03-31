#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 11:23:48 2026

@author: marine-louarn
"""

import spatialdata as sd
import sopa
from spatialdata_io import xenium
from pathlib import Path
import shutil
import spatialdata_plot
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns 
 

"""Spacial
"""
#sdata_ar = sopa.io.xenium("/home/marine-louarn/Documents/Xenium_Calderaro/Xenium_exp1/Sample1",cells_boundaries=True,nucleus_labels= True, cells_labels=True, nucleus_boundaries= True)

#sdata_ar.write("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1.zarr")
sdata = sd.read_zarr("/home/marine-louarn/Documents/Xenium_Calderaro/Sample1.zarr")

sopa.segmentation.tissue(sdata)
sopa.make_image_patches(sdata, patch_width=1200, patch_overlap=50)
sopa.utils.get_channel_names(sdata)
sopa.settings.parallelization_backend = "dask"
sopa.segmentation.cellpose(sdata, channels=["DAPI"], diameter=35) #, gpu=True)
sopa.make_transcript_patches(sdata, patch_width=500, prior_shapes_key="cellpose_boundaries")
sopa.segmentation.baysor(sdata, min_area=10)

sopa.aggregate(sdata)

sdata.pl.render_images("morphology_focus").pl.render_shapes(
    "cellpose_boundaries", outline_alpha=1, fill_alpha=0, outline_color="#ffffff"
).pl.show("global")
sdata.tables["table"].obs["region"] = "cell_boundaries"
sdata.set_table_annotates_spatialelement("table", region="cell_boundaries")
sdata.pl.render_shapes("cell_boundaries", color="total_counts").pl.show()

sdata.pl.render_points(
    "transcripts",
    color="feature_name",
    groups=["ERBB2", "A2ML1"],
    palette=["orange", "blue"],
).pl.show()

"""Transcript
"""
adata = sdata["table"]
adata.layers["counts"] = adata.X.copy()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True)
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)