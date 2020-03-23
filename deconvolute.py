#deconvolute.py
#Basil Khuder

import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sci
import seaborn as sns
from sklearn.svm import NuSVR
from autogenes import AutoGenes


def normalize_cells(counts, var_genes = 5000, pca_plot = False, scatter_plot = False):
    """Log normalize cells and focus downstream on 5000 variable genes. Default algorithm is Seurat"""
    counts_norm = sc.pp.normalize_per_cell(counts, copy=True) 
    counts_log = sc.pp.log1p(counts_norm, copy=True) 
    sc.pp.highly_variable_genes(counts_log, flavor='seurat', n_top_genes=var_genes + 1)
    counts_proc = counts_norm[:, counts_log.var[counts_log.var['highly_variable']==True].index]
    sc.pp.pca(counts_log, n_comps=30, use_highly_variable=True, svd_solver='arpack')
    if(pca_plot == True):
        print(sc.pl.pca_variance_ratio(counts_log, log=True))
    if(scatter_plot == True):
        counts_log.obs['cells'] = [x.split('.', 1)[0] for x in counts_log.obs_names]
        counts_log.obsm['X_pca'] *= -1
        print(sc.pl.pca_scatter(counts_log, color='cells'))
    return counts_proc[counts_log.obs_names]

def celltype_mean(counts, clusters):
    """Finds cell-type averages across"""
    sc_mean = pd.DataFrame(index=counts.var_names,columns=clusters)
    for cluster in clusters:
        cells = [x for x in counts.obs_names if x.startswith(cluster)]
        sc_part = counts[cells,:].X.T
        sc_mean[cluster] = pd.DataFrame(np.mean(sc_part,axis=1),index=counts.var_names)
    return sc_mean


def runAg(sc_mean, ngen = 5000, nfeatures = 400, print_plot = False):
    """Runs 5000 generations of optimizations using AutoGeneS and produces 400 marker genes to be used for deconvolution"""
    ag = AutoGenes(sc_mean.T)
    ag.run(ngen=ngen,seed=0,nfeatures=nfeatures,mode='fixed')
    if(print_plot == True):
        print(ag.plot(size='large',weights=(1,-1)))
    pareto = ag.pareto
    return sc_mean[pareto[len(pareto)-1]]


def produceProportions(ag, bulk_data, clusters):
    """Produces normalized cell-type proportion for each cell-type"""
    bulk_data = bulk_data.loc[ag.index,:]
    bulk_data = bulk_data.dropna()
    ag = ag.loc[bulk_data.index]
    proportions_NuSVR = pd.DataFrame(columns=clusters)
    for column in bulk_data.columns:
        regr_NuSVR = NuSVR(nu=0.5,C=0.5,kernel='linear')
        regr_NuSVR.fit(ag, bulk_data[column])
        proportions_NuSVR.loc[column] = regr_NuSVR.coef_[0]
    proportions_NuSVR[proportions_NuSVR < 0] = 0
    for raw in proportions_NuSVR.index:
        data_sum = proportions_NuSVR[raw].sum()
        proportions_NuSVR[raw] = np.divide(proportions_NuSVR.loc[raw],data_sum)
    print(proportions_NuSVR) 
