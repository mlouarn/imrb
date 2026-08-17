#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Random Tests
Created on 2026-08-12

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
import numpy as np
import pandas as pd
import novae
import biomart
from collections import defaultdict
import pyucell as uc
import functions_spatial as func

#test for different conda for novae computation (res1 : novae, spacial: spacialPy)
test_novae = sc.read_h5ad("/home/marine-louarn/Documents/Xenium_Calderaro/20260811_Sample1_novaeV2.h5ad")
df= test_novae.obs[['novae_domains_res1','novae_domains_res1_spacial']]
df["val"] = 1
newdf = pd.pivot_table(df, values="val", index="novae_domains_res1", columns="novae_domains_res1_spacial", aggfunc="sum")

sns.heatmap(newdf,annot=True)
plt.show()

df2=newdf.div(newdf.sum(axis=1), axis=0) #normalisation by row
sns.heatmap(df2,annot=True)
plt.show()

df3=newdf.T.div(newdf.T.sum(axis=1), axis=0) #normalisation by row
sns.heatmap(df3.T,annot=True)
plt.show()