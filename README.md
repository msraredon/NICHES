# NICHES
Niche Interactions and Cellular Heterogeneity in Extracellular Signaling

## Overview
NICHES is a simple but powerful computational toolset to analyze cell-cell signaling at the single-cell level. NICHES creates unique one-to-one pairs of cells and characterizes each cellular interaction across queried ligand-receptor signaling mechanisms. This allows low-dimensional embedding of cellular interactions in signal-space. Existing meta-data is carried over during transformation, allowing easy downstream statistical analysis using existing single-cell software packages including Seurat, Scanpy, Scater, and Monocle3.  NICHES allows UMAP, tSNE, and pseudotemporal exploration of cell-cell interactions and sensed cellular microenvironmments across celltypes, between experimental batches and conditions, and over developmental trajectories.

When applied to scRNAseq datasets, a uniform sampling of unique barcode pairings is taken from each celltype cross. When applied to spatial datasets, users may elect to constrain edge construction to cellular nearest-neighbors, or to neighbors within a certain radius, to provide a histologically-grounded portrait of cell-cell signaling and cellular niches. When system-level outputs are requested, NICHES computes all possible connectivity edges between all input cells.

Inputs and outputs to NICHES may be in either tabular or Seurat object format.

For detailed Methods, please see our preprint on bioRxiv: https://www.biorxiv.org/content/10.1101/2022.01.23.477401v1

## Installation
To install `NICHES` in R, you may run:
```
library(devtools)
install_github('msraredon/NICHES', ref = 'master')
```

## Vignettes
For use-case vignettes, please visit: https://msraredon.github.io/NICHES/

## Reference
If used for publication, please cite: https://doi.org/10.1101/2022.01.23.477401
