#' RunSCC
#' 
#' Performs Single-Cell Connectivity (SCC) transformations on a Seurat object. 
#' By default, references RunCellCell, RunCellSystem, and RunSystemCell functions. 
#' If positional data is provided, similar analyses can be performed which are limited exclusively to cells that directly neighbor each other.
#' Output is a set of specialized Seurat objects (SCC objects) corresponding to input arguments.
#' 
#' @param object A Seurat 4.0 object. The active identity meta.data will be used to define populations for connectomic sampling and crossings.
#' @param LR.database Accepts either 'fantom5' or a custom data.frame with the first column equal to ligands, second column equal to associated receptors.
#' @param species The species of the object that is being processed. Only required if LR.database = 'fantom5', and allows 'human','mouse','rat', or 'pig'
#' @param assay The assay to run the SCC transformation on. Defaults to "RNA."
#' @param min.cells.per.ident Default 10. A limit on how small (how many cells) a single population can be to participate in connectomic crossings.
#' @param meta.data.to.map A character vector of metadata names present in the original object which will be carried to the SCC objects
#' @param position.info The name of the meta.data column specifying spatial coordinates. Only relevant for spatial omics data.
#' @param ... Additional parameters to pass to RunCellCell, RunSystemCell, RunCellSystem, or spatial equivalents
#'
#' @export


RunSCC <- function(object,
                        LR.database = 'fantom5',
                        species,
                        assay = 'RNA',
                        min.cells.per.ident = 10,
                        meta.data.to.map = NULL,
                        position.info = NULL,
                        CellCell = T,
                        CellSystem = T,
                        SystemCell = T,
                        CellCellSpatial = T,
                        CellNeighborhood = T,
                        NeighborhoodCell = T,
                        ...){
  
  require(Seurat)
  require(dplyr)
  
  # Calculate SCC objects without spatial restrictions
  
  if (CellCell == T){obj1 <- RunCellCell(object)}
  if (CellSystem == T){obj2 <- RunCellSystem(object)}
  if (SystemCell == T){obj3 <- RunSystemCell(object)}
  
  # If requested, additionally calculate spatially-limited analogues
  
  if (CellCellSpatial == T | CellNeighborhood == T | NeighborhoodCell == T){
    
    if (is.null(position.info)){stop("\n Position information not provided. Please specify metadata column containing spatial coordinates.")}
    
    # Define distance between cells here, to use for radius calculations. Pass to the below three functions.
  
  }
  
  if (CellCellSpatial == T){obj4 <- RunCellCellSpatial(object)} #Spatially-limited Cell-Cell vectors
  if (CellNeighborhood == T){obj5 <- RunCellNeighborhood(object)} #Spatially-limited Cell-Neighborhood vectors
  if (NeighborhoodCell == T){obj6 <- RunNeighborhoodCell(object)} #Spatially-limited Neighborhood-Cell vectors (niches)

  # Compile objects for output

}

