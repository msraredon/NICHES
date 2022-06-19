# NICHES
Niche Interactions and Cellular Heterogeneity in Extracellular Signaling

## Overview
NICHES is a simple but powerful computational toolset to analyze cell-cell signaling at the truly single-cell level. Our approach goes beyond techniques which rely on cluster-wise averaging or shared neighbor graph binning to create paired ligand and receptor expression values.

## Brief Methods
NICHES creates unique 1:1 pairs of cells and characterizes each cellular interaction across queried ligand-receptor signaling mechanisms. This allows low-dimensional embedding of cellular interactions in signal-space. Existing meta-data is carried over during transformation, allowing easy and powerful downstream statistical analysis using existing single-cell software packages.  NICHES allows UMAP, tSNE, and pseudotemporal exploration of cell-cell interactions and sensed cellular microenvironmments across celltypes, between experimental conditions, and over developmental trajectories.

When applied to traditional scRNAseq datasets, a uniform sampling of unique barcode pairings is taken from each celltype cross. When applied to spatial datasets, edge construction may be limited to direct neighbor-neighbor interactions, or those within a certain user-defined radius, to provide a histologically-grounded portrait of cell-cell signaling and cellular niches.

Inputs and outputs to NICHES may be in either tabular or Seurat format.

For detailed Methods, please see our preprint on bioRxiv: https://www.biorxiv.org/content/10.1101/2022.01.23.477401v1

## Installation
To install `NICHES` in R, you may run:
```
library(devtools)
install_github('msraredon/NICHES', ref = 'master')
```

## Vignettes
For use vignettes, please visit: https://msraredon.github.io/NICHES/

## Reference
If used for publication, please cite: https://doi.org/10.1101/2022.01.23.477401
