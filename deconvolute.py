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


def normalize_cells(counts, var_genes = 5000, plot_pca = False, plot_umap = False, n_neighbors=15, n_pcs = None):
    counts_norm = sc.pp.normalize_per_cell(counts, copy=True) 
    counts_log = sc.pp.log1p(counts_norm, copy=True) 
    sc.pp.highly_variable_genes(counts_log, flavor='seurat', n_top_genes=var_genes + 1)
    counts_proc = counts_norm[:, counts_log.var[counts_log.var['highly_variable']==True].index]
    sc.pp.pca(counts_proc, n_comps=30, svd_solver='arpack')
    counts_proc.obs['cells'] = [x.split('-', 1)[0] for x in counts_proc.obs_names]
    if(plot_pca == True):
        sc.pl.pca_variance_ratio(counts_proc, log=True)
    if(plot_umap == True):
        sc.pp.neighbors(counts_proc, n_neighbors = 15, n_pcs = n_pcs)
        sc.tl.umap(counts_proc)
        sc.pl.umap(counts_proc, color = 'cells')
    return counts_proc


def run_ag(counts, clusters, ngen = 5000, nfeatures = 400, print_plot = False):
    sc_mean = pd.DataFrame(index=counts.var_names,columns=clusters)
    for i in clusters:
        cells = [x for x in counts.obs_names if x.startswith(i)]
        sc_part = counts[cells,:].X.T
        sc_mean[i] = pd.DataFrame(np.mean(sc_part,axis=1),index=counts.var_names)
    ag = AutoGenes(sc_mean.T)
    ag.run(ngen=ngen,seed=0,nfeatures=nfeatures,mode='fixed')
    if(print_plot == True):
        print(ag.plot(size='large',weights=(1,-1)))
    return sc_mean[ag.pareto[len(ag.pareto)-1]]

def celltype_correlation_map(ag):
    corr = pd.DataFrame(data = np.corrcoef(ag.T), columns = ag.columns, index = ag.columns)
    mask = np.zeros_like(corr)
    mask[np.triu_indices_from(mask)] = True
    return(sns.clustermap(np.abs(corr),cmap=sns.color_palette("GnBu", 1000), robust=True))
    
def marker_genes_mapped(ag,bulk_data, clusters):
    bulk_data = bulk_data.loc[ag.index,:]
    bulk_data = bulk_data.dropna()
    ag_subset = ag_one.loc[bulk_data.index]
    m_genes = "{} marker genes remaining. {}% of marker genes mapped".format(len(bulk_data.index), (len(bulk_data.index)/len(ag.index)*100))
    print(m_genes)
    

def produce_proportions(ag, bulk_data, clusters):
    bulk_data = bulk_data.loc[ag.index,:]
    bulk_data = bulk_data.dropna()
    ag = ag.loc[bulk_data.index]
    proportions_NuSVR = pd.DataFrame(columns=clusters)
    for i in bulk_data.columns:
        regr_NuSVR = NuSVR(nu=0.5,C=0.5,kernel='linear')
        regr_NuSVR.fit(ag, bulk_data[i])
        proportions_NuSVR.loc[i] = regr_NuSVR.coef_[0]
    proportions_NuSVR[proportions_NuSVR < 0] = 0
    for i in proportions_NuSVR.index:
        data_sum = proportions_NuSVR.loc[i].sum()
        proportions_NuSVR.loc[i] = np.divide(proportions_NuSVR.loc[i],data_sum)
    return(proportions_NuSVR) 
