# Deconvolution of Bulk RNA-Seq
Basil Khuder

-------
Being able to determine cellular heteorgenity in Bulk RNA-Seq data is possible through deconvolution. By taking a scRNA-Seq dataset made 
up of cells of the same tissue, we can find marker genes that will help significantly differentiate our bulk sample(s). 

*deconvolute.py* contains functions that will normalize single-cell datasets and find X amount of marker genes to be used to
deconvolute the bulk sample. It uses the following packages: 

- scanPy
- AnnData
- Seaborn
- Sklearn
- AutoGeneS
- MatplotLib

See the [Jupyter Notebook](https://github.com/basilkhuder/BulkDeconvolution/blob/master/deconvolute.ipynb) for a full run-down. 

## Usage

Import ```denvolute``` and read in single-cell counts via Scanpy. File should have cells classified by cell-type, and counts unnormalized. If cell-types are not uniquely numbered (for example, if all T Cells are listed as T Cells rather than T Cells 1, T Cells 2, ect) use ```.obs_names_make_unique()```

``` python
data = sc.read("counts.csv")
data.obs_names_make_unique()
```

Use ```deconvolute.normalize_cells()``` to log-transform counts. Set plot_pca and plot_umap to ```True``` to visualize PCA elbow plot and UMAP embedding:

``` python
data_proc = deconvolute.normalize_cells(data, var_genes = 2000, plot_pca = True, pca_umap = True)
```

```run_ag()``` calculates average expression across cells and uses AutoGeneS to generate list of marker genes to be used for deconvolution. ```Clusters``` is a numpy array that specifies names of all the cell-types. ```ngene``` is the number of optimization generations to run. 

``` python
data_ag = deconvolute.run_ag(data_proc, ngen = 2000, clusters =  np.array(['NK Cells', 'T Cells' ,'B Cells','DC Cells']))
```

```produce_proportions()``` uses AutoGeneS-generated marker genes to produce cell-type proportions.

``` python
deconvolute.produce_proportions(data_ag, clusters = np.array(['NK Cells', 'T Cells' ,'B Cells','DC Cells']))
```
