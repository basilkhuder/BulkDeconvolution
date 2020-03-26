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

class Deconvolute():
  def __init__(self, counts, clusters):
    self.counts = counts
    self.clusters = clusters
    self.ag_pareto = []
    self.proportions = []

  def normalize_cells(self, var_genes = 5000, plot_pca = False, plot_umap = False, n_neighbors=15, n_pcs = None):
    self.counts_norm = sc.pp.normalize_per_cell(self.counts, copy=True) 
    self.counts_log = sc.pp.log1p(self.counts, copy=True) 
    sc.pp.highly_variable_genes(self.counts_log, flavor='seurat', n_top_genes=var_genes)
    self.counts_norm = self.counts_norm[:, self.counts_log.var[self.counts_log.var['highly_variable']==True].index]
    sc.pp.pca(self.counts_norm, n_comps=30, svd_solver='arpack')
    self.counts_norm.obs['cells'] = [x.split('-', 1)[0] for x in self.counts_norm.obs_names]
    if(plot_pca == True):
        sc.pl.pca_variance_ratio(self.counts_norm, log=True)
    if(plot_umap == True):
        sc.pp.neighbors(self.counts_norm, n_neighbors = 15, n_pcs = n_pcs)
        sc.tl.umap(self.counts_norm)
        sc.pl.umap(self.counts_norm, color = 'cells')
      
  def run_ag(self, ngen = 5000, nfeatures = 400, nfeatures_increment = 10, print_plot = False):
    sc_mean = pd.DataFrame(index=self.counts_norm.var_names,columns=self.clusters)
    for i in self.clusters:
        index_cells = [x for x in self.counts_norm.obs_names if x.startswith(i)]
        sc_part = self.counts_norm[index_cells,:].X.T
        sc_mean[i] = pd.DataFrame(np.mean(sc_part,axis=1),index=self.counts_norm.var_names)
    self.ag = AutoGenes(sc_mean.T)
    if(type(nfeatures) == list):
        nfeatures = [i for i in range(nfeatures[0],nfeatures[1] + nfeatures_increment, nfeatures_increment)]
        for i in nfeatures:
            self.ag.run(ngen=ngen,seed=0,nfeatures=i,mode='fixed', verbose = False)
            if(print_plot == True):
                print(self.ag.plot(size='large',weights=(1,-1)))
            self.ag_pareto.append(sc_mean[self.ag.pareto[len(self.ag.pareto)-1]])
    else:
        self.ag.run(ngen=ngen,seed=0,nfeatures=nfeatures,mode='fixed', verbose = False)
        self.ag_pareto.append(sc_mean[self.ag.pareto[len(self.ag.pareto)-1]])
        
  def produce_proportions(self,bulk_data):
    proportions_NuSVR = pd.DataFrame(columns=self.clusters)
    for i in range(len(self.ag_pareto)):
        data_subset = bulk_data.loc[self.ag_pareto[i].index,:]
        data_subset = data_subset.dropna()
        ag = self.ag_pareto[i].loc[data_subset.index]
        for j in data_subset.columns:
            regr_NuSVR = NuSVR(nu=0.5,C=0.5,kernel='linear')
            regr_NuSVR.fit(ag, data_subset[j])
            proportions_NuSVR.loc[j] = regr_NuSVR.coef_[0]
        proportions_NuSVR[proportions_NuSVR < 0] = 0
        for k in proportions_NuSVR.index:
            data_sum = proportions_NuSVR.loc[k].sum()
            proportions_NuSVR.loc[k] = np.divide(proportions_NuSVR.loc[k],data_sum)
        self.proportions.append(proportions_NuSVR)
            
    
  def cell_corr_map(self):
    for i in range(len(self.ag_pareto)):
        corr = pd.DataFrame(data = np.corrcoef(self.ag_pareto[i].T), columns = self.ag_pareto[i].columns, index = self.ag_pareto[i].columns)
        mask = np.zeros_like(corr)
        mask[np.triu_indices_from(mask)] = True
        sns.clustermap(np.abs(corr),cmap=sns.color_palette("GnBu", 1000), robust=True)
        
  def __repr__(self):
    return("Deconvolute object with {} clusters, {} genes and {} cells".format(len(self.clusters),len(self.counts.obs),len(self.counts.var)))
    
