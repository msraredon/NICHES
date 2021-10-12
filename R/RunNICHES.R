#' RunNICHES
#' 
#' Performs user-selected NICHES transformations on a Seurat object. 
#' By default, references RunCellToCell, RunCellToSystem, and RunSystemToCell functions. 
#' If positional data is provided, similar analyses can be performed which are limited exclusively to cells that directly neighbor each other.
#' Output is a set of specialized NICHES objects which allow fine analysis of cell-cell interactions.
#' 
#' @param object A Seurat 4.0 object. The active identity will be used to define populations for connectomic sampling and crossings.
#' @param LR.database Currently accepts 'fantom5' or 'omnipath' 
#' @param species The species of the object that is being processed. Only required if LR.database = 'fantom5', and allows 'human','mouse','rat', or 'pig'
#' @param assay The assay to run the NICHES transformation on. Defaults to "RNA."
#' @param min.cells.per.ident Default 1. A limit on how small (how many cells) a single population can be to participate in connectomic crossings.
#' @param min.cells.per.gene Default 10. Limits analysis to interactions involving genes expressed above minimum threshold number of cells in the system. 
#' @param meta.data.to.map A character vector of metadata names present in the original object which will be carried to the NICHES objects
#' @param position.x The name of the meta.data column specifying location on the spatial x-axis. Only relevant for spatial omics data.
#' @param position.y The name of the meta.data column specifying location on the spatial y-axis. Only relevant for spatial omics data.
#' @param ... Additional parameters to pass to RunCellToCell, RunSystemToCell, RunCellToSystem, or spatial equivalents
#' @param rad.set Default 1. The radius in Euclidean space to consider local neighbors.
#' @param blend Choice of linear operator to combine edges in single-cell niche investigations. Defaults to "sum", also accepts "mean".
#' @param CellToCell Default TRUE. Whether to analyze cell-cell interactions without considering spatial coordinates.
#' @param CellToSystem Default FALSE. Whether to analyze summed signaling output to total system coming from each cell. Does not consider Euclidean coordinates.
#' @param SystemToCell Default FALSE. Whether to analyze summed signaling input from total system landing on each cell (cellular microenvironment/niche). Does not consider Euclidean coordinates. 
#' @param CellToCellSpatial Default FALSE. Whether to analyze cell-cell interactions between Euclidean neighbors. Only applicable in spatial datasets.
#' @param CellToNeighborhood Default FALSE. Whether to analyze summed signaling output to Euclidean neighbors. Only applicable in spatial datasets.
#' @param NeighborhoodToCell Default FALSE. Whether to analyze summed signaling input from Euclidean neighbors (cellular microenvironment/niche). Only applicable in spatial datasets.
#'
#' @export


RunNICHES <- function(object,
                        LR.database,
                        species,
                        assay,
                        min.cells.per.ident = 1,
                        min.cells.per.gene = 10,
                        meta.data.to.map = NULL,
                        position.x = NULL,
                        position.y = NULL,
                        rad.set = 1,
                        blend = 'sum',
                        CellToCell = T,
                        CellToSystem = F,
                        SystemToCell = F,
                        CellToCellSpatial = F,
                        CellToNeighborhood = F,
                        NeighborhoodToCell = F,
                        ...){
   # Initialize output structure
  output <- list()
  
  # Calculate NICHES organizations without spatial restrictions
  
  if (CellToCell == T){output[[length(output)+1]] <- RunCellToCell(object,LR.database,assay = assay,species = species,meta.data.to.map = meta.data.to.map,min.cells.per.ident = min.cells.per.ident,min.cells.per.gene = min.cells.per.gene,...)}
  if (CellToSystem == T){output[[length(output)+1]] <- RunCellToSystem(object,LR.database,assay = assay,species = species,meta.data.to.map = meta.data.to.map,min.cells.per.ident = min.cells.per.ident,min.cells.per.gene = min.cells.per.gene,blend = blend,...)}
  if (SystemToCell == T){output[[length(output)+1]] <- RunSystemToCell(object,LR.database,assay = assay,species = species,meta.data.to.map = meta.data.to.map,min.cells.per.ident = min.cells.per.ident,min.cells.per.gene = min.cells.per.gene,blend = blend,...)}
  
  # If requested, additionally calculate spatially-limited NICHES organizations
  
  if (CellToCellSpatial == T | CellToNeighborhood == T | NeighborhoodToCell == T){
    
    if (is.null(position.x) | is.null(position.y)){stop("\n Position information not provided. Please specify metadata columns containing x- and y-axis spatial coordinates.")}
    
  ## Define distance between cells here (?), to use for radius calculations. Pass to the below three functions. This will generalize the function to any spatial dataset (?)
  
  }
  
  if (CellToCellSpatial == T){output[[length(output)+1]] <- RunCellToCellSpatial(object,LR.database,assay = assay,species = species,position.x = position.x,position.y = position.y,meta.data.to.map = meta.data.to.map,rad.set = rad.set,min.cells.per.ident = min.cells.per.ident,min.cells.per.gene = min.cells.per.gene,...)} #Spatially-limited Cell-Cell vectors
  if (CellToNeighborhood == T){output[[length(output)+1]] <- RunCellToNeighborhood(object,LR.database,assay = assay,species = species,position.x = position.x,position.y = position.y,meta.data.to.map = meta.data.to.map,rad.set = rad.set,min.cells.per.ident = min.cells.per.ident,min.cells.per.gene = min.cells.per.gene,...)} #Spatially-limited Cell-Neighborhood vectors
  if (NeighborhoodToCell == T){output[[length(output)+1]] <- RunNeighborhoodToCell(object,LR.database,assay = assay,species = species,position.x = position.x,position.y = position.y,meta.data.to.map = meta.data.to.map,rad.set = rad.set,min.cells.per.ident = min.cells.per.ident,min.cells.per.gene = min.cells.per.gene,...)} #Spatially-limited Neighborhood-Cell vectors (niches)

  # Compile objects for output
  return(output)
}

