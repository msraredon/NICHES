% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{compute_edgelist}
\alias{compute_edgelist}
\title{Compute an edgelist based on the spatial coordinates}
\usage{
compute_edgelist(filtered.obj, position.x, position.y, k = 4, rad.set = NULL)
}
\arguments{
\item{filtered.obj}{A filtered Seurat object. The active identity will be used to define populations for connectomic sampling and crossings.}

\item{position.x}{string. Optional. Default: NULL. The name of the meta.data column specifying location on the spatial x-axis. Only required for spatial omics data.}

\item{position.y}{string. Optional. Default: NULL. The name of the meta.data column specifying location on the spatial y-axis. Only required for spatial omics data.}

\item{k}{integer. Optional. Default: 4. Number of neighbors in a knn graph. Used to compute a mutual nearest neighbor graph based on the spatial coordinates of the spatial transcriptomic datasets.}

\item{rad.set}{numeric. Optional. Default: NULL. The radius threshold to define neighbors based on the spatial coordinates of the spatial transcriptomic datasets. Ignored when 'k' is provided.}
}
\description{
Compute an edgelist based on the spatial coordinates
}
