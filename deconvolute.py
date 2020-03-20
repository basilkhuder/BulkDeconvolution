def normalize_cells(counts, var_genes = 5000, pca_plot = False, scatter_plot = False):
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
    
def celltype_mean(clusters, counts):
    sc_mean = pd.DataFrame(index=counts.var_names,columns=clusters)
    for cluster in clusters:
        cells = [x for x in counts.obs_names if x.startswith(cluster)]
        sc_part = counts[cells,:].X.T
        sc_mean[cluster] = pd.DataFrame(np.mean(sc_part,axis=1),index=counts.var_names)
    return sc_mean
    
def runAg(sc_mean, ngen = 5000, nfeatures = 400, print_plot = False):
    ag = AutoGenes(sc_mean.T)
    ag.run(ngen=ngen,seed=0,nfeatures=nfeatures,mode='fixed')
    if(print_plot == True):
        print(ag.plot(size='large',weights=(1,-1)))
    pareto = ag.pareto
    return sc_mean[pareto[len(pareto)-1]]

def normalize_proportions(data):
    data[data < 0] = 0
    for raw in data.index:
        sum = data.loc[raw].sum()
        data.loc[raw] = np.divide(data.loc[raw],sum)
    return data
