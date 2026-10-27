import sys
import os
import tempfile
import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scvi
import scviva
import spatialdata as sd
import pandas as pd
from scviva.model._resolvi import ResolVI

adata_p1 = sc.read("home/marine-louarn/Documents/Lung_Spatial/Patient1/20261015_Patient1_withAllSubtypes.h5ad")

#model training for cType id
scviva.model.ResolVI.setup_anndata(adata_p1, labels_key="predicted_celltype", layer="counts")
supervised_resolvi = scviva.model.ResolVI(adata_p1, semisupervised=True)
supervised_resolvi.train(max_epochs=100)
adata_p1.obsm["resolvi_celltypes"] = supervised_resolvi.predict(adata_p1, num_samples=3, soft=True)
adata_p1.obs["resolvi_predicted"] = adata_p1.obsm["resolvi_celltypes"].idxmax(axis=1)
adata_p1.obsm["X_resolVI"] = supervised_resolvi.get_latent_representation(adata_p1)

sc.pl.umap(adata_p1, color=["predicted_celltype", "resolvi_predicted"])
sc.pl.spatial(adata_p1, color=["predicted_celltype", "resolvi_predicted"], spot_size=30)


samples_corr = supervised_resolvi.sample_posterior(
    model=supervised_resolvi.module.model_corrected,
    return_sites=["px_rate"],
    summary_fun={"post_sample_q50": np.median},
    num_samples=3,
    summary_frequency=30,
)
samples_corr = pd.DataFrame(samples_corr).T

samples = supervised_resolvi.sample_posterior(
    model=supervised_resolvi.module.model_residuals,
    return_sites=["mixture_proportions"],
    summary_fun={"post_sample_means": np.mean},
    num_samples=3,
    summary_frequency=100,
)
samples = pd.DataFrame(samples).T
adata_p1.obs[["true_proportion", "diffusion_proportion", "background_proportion"]] = samples.loc[ "post_sample_means", "mixture_proportions"]

sc.pl.umap(adata_p1,color=["total_counts", "true_proportion", "diffusion_proportion", "background_proportion"],)
adata_p1.layers["generated_expression"] = samples_corr.loc["post_sample_q50", "px_rate"]

sc.pl.umap(adata_p1, color=["resolvi_predicted",'CD79A'], layer="counts", vmax="p98")
sc.pl.umap(adata_p1, color=["resolvi_predicted", "CD79A"], layer="generated_expression", vmax="p98")

adata_p1.write_h5ad("/home/marine-louarn/Documents/Lung_Spatial/Patient1/20261026_Patient1_resolVI.h5ad")