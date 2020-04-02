# Deconvolution of Bulk RNA-Seq
Basil Khuder

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

Import ```deconvolute``` and read in single-cell counts and bulk counts via Scanpy. File should have cells classified by cell-type, and counts unnormalized. If cell-types are not uniquely numbered (for example, if all T Cells are listed as T Cells rather than T Cells 1, T Cells 2, ect) use ```.obs_names_make_unique()```

``` python
import pandas as pd
import numpy as np
import scanpy as sc
from deconvolute import Deconvolute
data = sc.read("counts.csv")
bulk = sc.read("bulk.csv")
data.obs_names_make_unique()
```
Create a new Deconvolute object with the imported counts and clusters. Clusters is a numpy array specifying all of the cell-types.

``` python
dc = Deconvolute(sc_counts = data, bulk_counts = bulk, clusters =  np.array(['NK Cells', 'T Cells' ,'B Cells','DC Cells']))
```

Use ```normalize_cells()``` to log-transform counts and set amount of variable genes for downstream analysis. Set plot_pca and plot_umap to ```True``` to visualize PCA elbow plot and UMAP embedding:

``` python
dc.normalize_cells(var_genes = 2000, plot_pca = True, plot_umap = True)
```

```run_ag()``` calculates average expression across cells and uses AutoGeneS to generate list of marker genes to be used for deconvolution. ```ngen``` is the number of optimization generations to run, ```nfeatures``` is the number of marker genes to select. If you want to look at a range of marker genes, put them in a list and use ```nfeatures_increment``` to specficy the increments:

``` python
#For only 200 marker genes
dc.run_ag(ngen = 2000, nfeatures = 200)
#For 200 to 1000 marker genes in increments of 200
dc.run_ag(ngen = 2000, nfeatures = [200,1000], nfeatures_increment = 200)
```

```produce_proportions(bulk_data)``` uses AutoGeneS-generated marker genes to produce cell-type proportions on bulk RNA-Seq counts. 

``` python
bulk_data = pd.read_csv("bulk_counts.txt", index_col = 0, sep = "\t")
dc.produce_proportions()
```
