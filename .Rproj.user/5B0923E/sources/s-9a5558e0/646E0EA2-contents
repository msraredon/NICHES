#### Define the functions from Seurat to SCC
#### remember to export
#setGeneric("Idents<-",function(object,organization,...){standardGeneric("Idents<-")})
#setMethod("Idents<-", signature(object = "SCC"), function(object,organization,...) {
#  Seurat::`Idents<-`(object@organizations[[organization]],...)
#})


#' FindVariableFeatures for SCC
#'
#' @param object an SCC object
#' @param organization the organization to run on
#' @param ... parameters passed to base function in Seurat
#'
#' @return an SCC object
#' @export
#'
setGeneric("FindVariableFeatures",function(object,organization,...){standardGeneric("FindVariableFeatures")})
setMethod("FindVariableFeatures", signature(object = "SCC"), function(object,organization,...) {
  object@organizations[[organization]] <-  Seurat::FindVariableFeatures(object@organizations[[organization]],...)
  return(object)
})

# vars <- VariableFeatures(object = scc_tmp,organization = "pair")
#' VariableFeatures
#'
#' @param object an SCC object
#' @param organization the organization to run on
#' @param ... parameters passed to base function in Seurat
#'
#' @return an SCC object
#' @export
#'

setGeneric("VariableFeatures",function(object,organization,...){standardGeneric("VariableFeatures")})
setMethod("VariableFeatures", signature(object = "SCC"), function(object,organization,...) {
  return(Seurat::VariableFeatures(object@organizations[[organization]],...))
})

# Usage: scc_tmp <- ScaleData(scc_tmp,organization = "pair")
#' ScaleData
#'
#' @param object an SCC object
#' @param organization the organization to run on
#' @param ... parameters passed to base function in Seurat
#'
#' @return SCC object
#' @export
#'

setGeneric("ScaleData",function(object,organization,...){standardGeneric("ScaleData")})
setMethod("ScaleData", signature(object = "SCC"), function(object,organization,...) {
  object@organizations[[organization]] <-  Seurat::ScaleData(object@organizations[[organization]],...)
  return(object)
})

# scc_tmp <- RunPCA(scc_tmp,organization = "pair",features = VariableFeatures(object = scc_tmp,organization = "pair"),seed.use = 34)
#' RunPCA
#'
#' @param object an SCC object
#' @param organization the organization to run on
#' @param ... parameters passed to base function in Seurat
#'
#' @return SCC object
#' @export
#'

setGeneric("RunPCA",function(object,organization,...){standardGeneric("RunPCA")})
setMethod("RunPCA", signature(object = "SCC"), function(object,organization,...) {
  object@organizations[[organization]] <-  Seurat::RunPCA(object@organizations[[organization]],...)
  return(object)
})

#' RunTSNE
#'
#' @param object an SCC object
#' @param organization the organization to run on
#' @param ... parameters passed to base function in Seurat
#'
#' @return an SCC object
#' @export
#'
setGeneric("RunTSNE",function(object,organization,...){standardGeneric("RunTSNE")})
setMethod("RunTSNE", signature(object = "SCC"), function(object,organization,...) {
  object@organizations[[organization]] <-  Seurat::RunTSNE(object@organizations[[organization]],...)
  return(object)
})

# Usage: DimPlot(SCC_object,organization = "pair_spatial",...)

#' DimPlot
#'
#' @param object an SCC object
#' @param organization the organization to run on
#' @param ... parameters passed to base function in Seurat
#'
#' @return an SCC object
#' @export
#'

setGeneric("DimPlot",function(object,organization,...){standardGeneric("DimPlot")})
setMethod("DimPlot", signature(object = "SCC"), function(object,organization,...) {
  Seurat::DimPlot(object@organizations[[organization]],...)
})


#' FindNeighbors
#'
#' @param object an SCC object
#' @param organization the organization to run on
#' @param ... parameters passed to base function in Seurat
#'
#' @return an SCC object
#' @export
#'

setGeneric("FindNeighbors",function(object,organization,...){standardGeneric("FindNeighbors")})
setMethod("FindNeighbors", signature(object = "SCC"), function(object,organization,...) {
  object@organizations[[organization]] <-  Seurat::FindNeighbors(object@organizations[[organization]],...)
  return(object)
})


#' FindClusters
#'
#' @param object an SCC object
#' @param organization the organization to run on
#' @param ... parameters passed to base function in Seurat
#'
#' @return an SCC object
#' @export
#'

setGeneric("FindClusters",function(object,organization,...){standardGeneric("FindClusters")})
setMethod("FindClusters", signature(object = "SCC"), function(object,organization,...) {
  object@organizations[[organization]] <-  Seurat::FindClusters(object@organizations[[organization]],...)
  return(object)
})


#' FindAllMarkers
#'
#' @param object an SCC object
#' @param organization the organization to run on
#' @param ... parameters passed to base function in Seurat
#'
#' @return an SCC object
#' @export
#'

setGeneric("FindAllMarkers",function(object,organization,...){standardGeneric("FindAllMarkers")})
setMethod("FindAllMarkers", signature(object = "SCC"), function(object,organization,...) {
  return(Seurat::FindAllMarkers(object@organizations[[organization]],...))
})


#' FeaturePlot
#'
#' @param object an SCC object
#' @param organization the organization to run on
#' @param ... parameters passed to base function in Seurat
#'
#' @return an SCC object
#' @export
#'

setGeneric("FeaturePlot",function(object,organization,...){standardGeneric("FeaturePlot")})
setMethod("FeaturePlot", signature(object = "SCC"), function(object,organization,...) {
  return(Seurat::FeaturePlot(object@organizations[[organization]],...))
})


#' DotPlot
#'
#' @param object an SCC object
#' @param organization the organization to run on
#' @param ... parameters passed to base function in Seurat
#'
#' @return an SCC object
#' @export
#'
setGeneric("DotPlot",function(object,organization,...){standardGeneric("DotPlot")})
setMethod("DotPlot", signature(object = "SCC"), function(object,organization,...) {
  return(Seurat::DotPlot(object@organizations[[organization]],...))
})

