#' RunSCC
#' 
#' Performs Single-Cell Connectivity (SCC) transformations on a Seurat object. 
#' By default, references RunCellToCell, RunCellToSystem, and RunSystemToCell functions. 
#' If positional data is provided, similar analyses can be performed which are limited exclusively to cells that directly neighbor each other.
#' Output is a set of specialized Seurat objects (SCC objects) containing cell signaling information.
#' 
#' @param object A Seurat 4.0 object. The active identity will be used to define populations for connectomic sampling and crossings.
#' @param LR.database Accepts either 'fantom5' or a custom data.frame with the first column equal to ligands, second column equal to associated receptors.
#' @param species The species of the object that is being processed. Only required if LR.database = 'fantom5', and allows 'human','mouse','rat', or 'pig'
#' @param assay The assay to run the SCC transformation on. Defaults to "RNA."
#' @param min.cells.per.ident Default 10. A limit on how small (how many cells) a single population can be to participate in connectomic crossings.
#' @param meta.data.to.map A character vector of metadata names present in the original object which will be carried to the SCC objects
#' @param position.x The name of the meta.data column specifying location on the spatial x-axis. Only relevant for spatial omics data.
#' @param position.y The name of the meta.data column specifying location on the spatial y-axis. Only relevant for spatial omics data.
#' @param ... Additional parameters to pass to RunCellToCell, RunSystemToCell, RunCellToSystem, or spatial equivalents
#'
#' @export


RunSCC <- function(object,
                        LR.database = 'fantom5',
                        species,
                        assay = 'RNA',
                        min.cells.per.ident = 10,
                        meta.data.to.map = NULL,
                        position.x = NULL,
                        position.y = NULL,
                        CellToCell = T,
                        CellToSystem = F,
                        SystemToCell = T,
                        CellToCellSpatial = T,
                        CellToNeighborhood = F,
                        NeighborhoodToCell = T,
                        ...){
   # Initialize output structure
  output <- list()
  
  # Calculate SCC organizations without spatial restrictions
  
  if (CellToCell == T){output[[length(output)+1]] <- RunCellToCell(object,species = species,meta.data.to.map = meta.data.to.map,...)}
  if (CellToSystem == T){output[[length(output)+1]] <- RunCellToSystem(object,species = species,meta.data.to.map = meta.data.to.map,...)}
  if (SystemToCell == T){output[[length(output)+1]] <- RunSystemToCell(object,species = species,meta.data.to.map = meta.data.to.map,...)}
  
  # If requested, additionally calculate spatially-limited SCC organizations
  
  if (CellToCellSpatial == T | CellToNeighborhood == T | NeighborhoodToCell == T){
    
    if (is.null(position.x) | is.null(position.y)){stop("\n Position information not provided. Please specify metadata columns containing x- and y-axis spatial coordinates.")}
    
  ## Define distance between cells here (?), to use for radius calculations. Pass to the below three functions. This will generalize the function to any spatial dataset (?)
  
  }
  
  if (CellToCellSpatial == T){output[[length(output)+1]] <- RunCellToCellSpatial(object,species = species,position.x = position.x,position.y = position.y,meta.data.to.map = meta.data.to.map,...)} #Spatially-limited Cell-Cell vectors
  if (CellToNeighborhood == T){output[[length(output)+1]] <- RunCellToNeighborhood(object,species = species,position.x = position.x,position.y = position.y,meta.data.to.map = meta.data.to.map,...)} #Spatially-limited Cell-Neighborhood vectors
  if (NeighborhoodToCell == T){output[[length(output)+1]] <- RunNeighborhoodToCell(object,species = species,position.x = position.x,position.y = position.y,meta.data.to.map = meta.data.to.map,...)} #Spatially-limited Neighborhood-Cell vectors (niches)

  # Compile objects for output
  return(output)
}

