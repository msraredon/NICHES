% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunNICHES.R
\name{RunNICHES}
\alias{RunNICHES}
\alias{RunNICHES.default}
\alias{RunNICHES.Seurat}
\alias{RunNICHES.matrix}
\alias{RunNICHES.Matrix}
\title{RunNICHES}
\usage{
RunNICHES(object, ...)

\method{RunNICHES}{default}(
  object,
  assay = "RNA",
  LR.database = "fantom5",
  species,
  min.cells.per.ident = NULL,
  min.cells.per.gene = NULL,
  meta.data.to.map = NULL,
  position.x = NULL,
  position.y = NULL,
  custom_LR_database = NULL,
  k = 4,
  rad.set = NULL,
  blend = "mean",
  CellToCell = T,
  CellToSystem = F,
  SystemToCell = F,
  CellToCellSpatial = F,
  CellToNeighborhood = F,
  NeighborhoodToCell = F,
  output_format = "seurat",
  ...
)

\method{RunNICHES}{Seurat}(
  object,
  assay = "RNA",
  LR.database = "fantom5",
  species,
  min.cells.per.ident = NULL,
  min.cells.per.gene = NULL,
  meta.data.to.map = colnames(object@meta.data),
  position.x = NULL,
  position.y = NULL,
  cell_types = NULL,
  custom_LR_database = NULL,
  k = 4,
  rad.set = NULL,
  blend = "mean",
  CellToCell = T,
  CellToSystem = F,
  SystemToCell = F,
  CellToCellSpatial = F,
  CellToNeighborhood = F,
  NeighborhoodToCell = F,
  output_format = "seurat",
  ...
)

\method{RunNICHES}{matrix}(
  object,
  LR.database = "fantom5",
  species,
  min.cells.per.ident = NULL,
  min.cells.per.gene = NULL,
  meta.data.df = NULL,
  position.x = NULL,
  position.y = NULL,
  cell_types = NULL,
  custom_LR_database = NULL,
  k = 4,
  rad.set = NULL,
  blend = "mean",
  CellToCell = T,
  CellToSystem = F,
  SystemToCell = F,
  CellToCellSpatial = F,
  CellToNeighborhood = F,
  NeighborhoodToCell = F,
  output_format = "seurat",
  ...
)

\method{RunNICHES}{Matrix}(
  object,
  LR.database = "fantom5",
  species,
  min.cells.per.ident = NULL,
  min.cells.per.gene = NULL,
  meta.data.df = NULL,
  position.x = NULL,
  position.y = NULL,
  cell_types = NULL,
  custom_LR_database = NULL,
  k = 4,
  rad.set = NULL,
  blend = "mean",
  CellToCell = T,
  CellToSystem = F,
  SystemToCell = F,
  CellToCellSpatial = F,
  CellToNeighborhood = F,
  NeighborhoodToCell = F,
  output_format = "seurat",
  ...
)
}
\arguments{
\item{object}{A normalized data matrix. Rows are genes and columns are cells.}

\item{...}{Arguments passed for other methods}

\item{assay}{string. Default: "RNA". The assay to run the NICHES transformation on.}

\item{LR.database}{string. Default: "fantom5". Currently accepts "fantom5","omnipath", or "custom".}

\item{species}{string. The species of the object that is being processed. Only required when LR.database = 'fantom5' with species being 'human','mouse','rat', or 'pig', or LR.database = 'omnipath' with species being 'human','mouse', or 'rat'.}

\item{min.cells.per.ident}{integer. Default: NULL. A limit on how small (how many cells) a single population can be to participate in connectomic crossings.}

\item{min.cells.per.gene}{integer. Default: NULL. Limits analysis to interactions involving genes expressed above minimum threshold number of cells in the system.}

\item{meta.data.to.map}{character vector. Optional. Default: NULL. A vector of metadata names present in the original object which will be carried to the NICHES objects}

\item{position.x}{string. Optional. Default: NULL. Only required for spatial omics data. The name that specifies location on the spatial x-axis in the input meta.data.df.}

\item{position.y}{string. Optional. Default: NULL. Only required for spatial omics data. The name that specifies location on the spatial y-axis in the input meta.data.df.}

\item{custom_LR_database}{data.frame. Optional. Default: NULL. Only required when LR.database = "custom". Each row is a ligand-receptor mechanism where the first column corresponds to the source genes that express the ligands subunits (separated by '_') and the second column corresponds to the receptor genes that express the receptor subunits (separated by '_').}

\item{k}{integer. Optional. Default: 4. Number of neighbors in a knn graph. Used to compute a mutual nearest neighbor graph based on the spatial coordinates of the spatial transcriptomic datasets.}

\item{rad.set}{numeric. Optional. Default: NULL. The radius threshold to define neighbors based on the spatial coordinates of the spatial transcriptomic datasets. Ignored when 'k' is provided.}

\item{blend}{string. Default: "mean". Choice of linear operator to combine edges in single-cell niche investigations. Defaults to "mean", also accepts "sum".}

\item{CellToCell}{logical. Default: TRUE. Whether to analyze cell-cell interactions without considering spatial coordinates.}

\item{CellToSystem}{logical. Default: FALSE. Whether to analyze summed signaling output to total system coming from each cell. Does not consider Euclidean coordinates.}

\item{SystemToCell}{logical. Default: FALSE. Whether to analyze summed signaling input from total system landing on each cell (cellular microenvironment/niche). Does not consider Euclidean coordinates.}

\item{CellToCellSpatial}{logical. Default: FALSE. Whether to analyze cell-cell interactions between Euclidean neighbors. Only applicable in spatial datasets.}

\item{CellToNeighborhood}{logical. Default: FALSE. Whether to analyze summed signaling output to Euclidean neighbors. Only applicable in spatial datasets.}

\item{NeighborhoodToCell}{logical. Default: FALSE. Whether to analyze summed signaling input from Euclidean neighbors (cellular microenvironment/niche). Only applicable in spatial datasets.}

\item{output_format}{string. Default: "seurat". Choice of the output format. "seurat" will output a list of seurat objects, "raw" will output a list of lists with raw interaction matrix and compiled metadata}

\item{cell_types}{string. Default: NULL. The name of the meta.data.df column specifying the cell types of the cells.}

\item{meta.data.df}{A dataframe. Optional. Default: NULL. A dataframe of the metadata (columns) associated with the cells (rows).}
}
\description{
Performs user-selected NICHES transformations on a Seurat object or a normalized data matrix. 
By default, references RunCellToCell, RunCellToSystem, and RunSystemToCell functions. 
If positional data is provided, similar analyses can be performed which are limited exclusively to cells that directly neighbor each other.
Output is a set of specialized NICHES objects which allow fine analysis of cell-cell interactions.
}
