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

## Usage

Import ```denvolute``` and read in single-cell counts via Scanpy (use appropriate method depending on counts file format). File should have cells classified by cell-type, and counts unnormalized. If cell-types are not uniquely numbered (for example, if all T Cells are listed as T Cells rather than T Cells 1, T Cells 2, ect) use ```.obs_names_make_unique()```

``` python
data = sc.read_csv("counts.csv")
data.obs_names_make_unique()
```

Use ```deconvolute.normalize_cells()``` to log-transform counts. Set pca_plot and scatter_plot to ```True``` to visualize PCA elbow plot and clustering:

``` python
data_proc = deconvolute.normalize_cells(data, var_genes = 4000, scatter_plot = True, pca_plot = True)
```

Use ```celltype_mean()``` to calculate average-expression across cells. Clusters is a numpy array that specifies names of all the cell-types. 

``` python

data_proc_mean = deconvolute.celltype_mean(data, clusters = np.array(['NK Cells', 'T Cells' ,'B Cells','DC Cells']))
```

```run_ag()``` uses AutoGeneS to generate list of marker genes to be used for deconvolution. 

``` python
data_ag = runAg(mean_counts)
```

```produce_propotion()``` uses AutoGeneS-generated marker genes to produce cell-type proportions.

``` python
produce_proption(data_ag, clusters = np.array(['NK Cells', 'T Cells' ,'B Cells','DC Cells']))
```
