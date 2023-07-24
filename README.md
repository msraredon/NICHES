# NICHES
Niche Interactions and Cellular Heterogeneity in Extracellular Signaling

## Overview
NICHES is a powerful computational toolset to analyze cell-cell signaling at the single-cell level. NICHES creates unique one-to-one pairs of cells and characterizes the resulting cellular relationships across all queried ligand-receptor signaling mechanisms. This allows dimensional reduction and subsequent vizualization of cellular interactions relative to one another in signal-space. 

We recommend carrying over all existing meta-data when running NICHES, which facilitates easy downstream data manipulation and statistical analysis. NICHES output data can be analyzed using existing single-cell software packages in both R and python, including Seurat, Scanpy, Scater, and Monocle3.  NICHES allows rapid creation of UMAP, tSNE, pseudotemporal and ComplexHeatmap-based exploration of cell-cell interactions and sensed cellular microenvironmments across celltypes, between experimental batches and conditions, and over developmental trajectories.

When applied to scRNAseq datasets without spatial coordinates, a uniform sampling of unique barcode pairings is taken from each celltype cross so that all possible celltype crosses are sampled. When applied to spatial datasets, users may elect to constrain edge construction to spot-wise nearest-neighbors, or to neighbors within a certain radius, to provide a histologically-grounded portrait of cell-cell signaling and cellular niches. These spatially-informed pairings can then be visualized directly in histologic space. When system-level outputs are requested without spatial coordinates (SystemToCell and CellToSystem), NICHES considers all connectivity edges between all input cells for a given system.

Inputs and outputs to NICHES may be in either tabular or Seurat object format, which facilitates intercompatibility across platforms and packages. 

For detailed Methods, please see our publication in Bioinformatics: https://doi.org/10.1093/bioinformatics/btac775

## Installation
To install `NICHES` in R, you may run:
```
library(devtools)
install_github('msraredon/NICHES', ref = 'master')
```

## Vignettes
For use-case vignettes, please visit: https://msraredon.github.io/NICHES/articles/

## Reference
If used for publication, please cite: https://doi.org/10.1093/bioinformatics/btac775
