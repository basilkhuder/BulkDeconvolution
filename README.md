# Deconvolution of Bulk RNA-Seq
Basil Khuder
Updated: 4/2/20

-------
Being able to determine cellular heteorgenity in Bulk RNA-Seq data is possible through deconvolution. By taking a scRNA-Seq dataset made 
up of cells of the same tissue, we can find marker genes that will help significantly differentiate our bulk sample(s). 

*deconvolute.py* allows for normalization of single-cell datasets to find an appropriate and representative amount of  marker genes to be used to deconvolute a bulk RNA-Seq sample. It uses the following packages: 

- scanPy
- AnnData
- Seaborn
- Sklearn
- AutoGeneS

## Usage

Import ```deconvolute``` and create a new Deconvolute object with your single-cell reference and bulk counts. Single-cell reference should have cells classified by cell-type, and counts unnormalized.

``` python
import pandas as pd
import numpy as np
import scanpy as sc
from deconvolute import Deconvolute
dc = Deconvolute(sc_counts = "counts.csv", bulk = "bulk.csv")
```
Use ```normalize_cells()``` to perform library-size correction, log-transformation and variable gene identification. Set ```plot_pca``` and ```plot_umap``` to ```True``` to visualize PCA elbow plot and UMAP embedding:

``` python
dc.normalize_cells(var_genes = 2000, plot_pca = True, plot_umap = True)
```

```run_ag()``` calculates average expression across cells and uses AutoGeneS to generate list of marker genes to be used for deconvolution. 

```gen``` is the number of generations to run using AutoGeneS which aims at maximizing the Euclidean distance, while minimizing the correlation of the set of marker genes. 

```markers``` is the number of marker genes to select. If you want to look at a range of marker genes, put them in a list and use ```incr``` to specficy the increments:

``` python
#For only 200 marker genes
dc.run_ag(gen = 2000, markers = 200)
#For 200 to 1000 marker genes in increments of 200
dc.run_ag(gen = 2000, markers = [200,1000], incr = 200)
```

```produce_proportions(bulk_data)``` uses AutoGeneS-generated marker genes to produce cell-type proportions on bulk RNA-Seq counts. 

``` python
dc.produce_proportions()
```
