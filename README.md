# NICHES
Niche Interactions and Cellular Heterogeneity in Extracellular Signaling

## Overview
NICHES is a simple but powerful computational toolset to analyze cell-cell signaling at the truly single-cell level. Our approach goes beyond techniques which rely on cluster-wise averaging or shared neighbor graph binning to create paired ligand and receptor expression values.

## Brief Methods
NICHES creates unique 1:1 pairs of cells and characterizes each cellular interaction across all queried ligand-receptor signaling mechanisms. This allows low-dimensional embedding of cellular interactions in signal-space. All existing meta-data is carried over during transformation, allowing easy and powerful downstream statistical analysis using existing single-cell software packages such as Seurat and Scanpy.  Niche-interactions within a system may be characterized by summing all signaling edges, by mechanism, landing on a given receiving cell. These two approaches allow UMAP, tSNE, and pseudotemporal visualization of cell-cell interactions and cellular niches across celltypes, between experimental conditions, and over pseudotemporal trajectories.

When applied to traditional scRNAseq datasets, a uniform sampling of unique barcode pairings is taken from each celltype cross. When applied to spatial datasets, edge construction may be limited to direct neighbor-neighbor interactions, or those within a certain user-defined radius, to provide a histologically-grounded portrait of cell-cell signaling and cellular niches.

For detailed Methods, please see our preprint on bioRxiv:

## Installation
To install `NICHES` in R, you may run:
```
library(devtools)
install_github('msraredon/Connectome', ref = 'master')
```

## Reference
If used for publication, please cite:
