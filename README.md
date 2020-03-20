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
