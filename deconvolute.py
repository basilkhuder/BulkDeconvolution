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
    def __init__(self, sc_counts, bulk):
        if sc_counts.endswith('.csv'):
            self.scc = pd.read_csv(sc_counts, index_col = 0)
        if sc_counts.endswith('.txt'):
            self.scc = pd.read_csv(sc_counts, index_col = 0, sep = "\t")
        self.clusters = np.unique(np.array([x.split(".")[0] for x in list(self.scc.columns)]))
        if bulk.endswith('.csv'):
            self.bulk = pd.read_csv(bulk, index_col = 0)
        if bulk.endswith('.txt'):
            self.bulk = pd.read_csv(bulk, index_col = 0, sep = "\t")
        self.scc = sc.AnnData(self.scc.T)
        self.samples = np.array(self.bulk.columns)
        self.ag = []
        self.pro = []


    def normalize_cells(self, var_genes=2000, plot_pca=False):
        sc.pp.normalize_total(self.scc)
        self.log_counts = sc.pp.log1p(self.scc, copy=True)
        sc.pp.highly_variable_genes(self.log_counts, flavor='seurat', n_top_genes=var_genes + 1)
        self.scc = self.scc[:, self.log_counts.var[self.log_counts.var['highly_variable'] == True].index]
        sc.pp.pca(self.log_counts, n_comps=30, use_highly_variable=True, svd_solver='arpack')
        if (plot_pca == True):
            print(sc.pl.pca_variance_ratio(self.log_counts, log=True))
        self.scc = self.scc[self.log_counts.obs_names]


    def run_ag(self, gen=5000, markers=400, incr=0):
        cm = pd.DataFrame(index=self.scc.var_names, columns=self.clusters)
        for i in self.clusters:
            cells = [x for x in self.scc.obs_names if x.startswith(i)]
            cm[i] = pd.DataFrame(np.mean(self.scc[cells, :].X, axis=0), index=self.scc.var_names)
        ag = AutoGenes(cm.T)
        if (type(markers) == list):
            markers = [i for i in range(markers[0], markers[1] + incr, incr)]
            for i in markers:
                ag.run(ngen=gen, seed=0, nfeatures=i, mode='fixed')
                self.ag.append(cm[ag.pareto[len(ag.pareto) - 1]])
        else:
            ag.run(ngen=gen, seed=0, nfeatures=markers, mode='fixed')
            self.ag.append(cm[ag.pareto[len(ag.pareto) - 1]])


    def produce_proportions(self):
        for i in range(len(self.ag)):
            pro_df = pd.DataFrame(columns=self.clusters)
            bulk_sub = self.bulk.loc[self.ag[i].index, :]
            bulk_sub = bulk_sub.dropna()
            self.ag[i] = self.ag[i].loc[bulk_sub.index]
            for j in bulk_sub.columns:
                regr = NuSVR(nu=0.5, C=0.5, kernel='linear')
                regr.fit(self.ag[i], bulk_sub[j])
                pro_df.loc[j] = regr.coef_[0]
            pro_df[pro_df < 0] = 0
            for k in pro_df.index:
                summ = pro_df.loc[k].sum()
                pro_df.loc[k] = np.divide(pro_df.loc[k], summ)
            self.pro.append(pro_df)
