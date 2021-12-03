#' RunNICHES
#' 
#' Performs user-selected NICHES transformations on a Seurat object. 
#' By default, references RunCellToCell, RunCellToSystem, and RunSystemToCell functions. 
#' If positional data is provided, similar analyses can be performed which are limited exclusively to cells that directly neighbor each other.
#' Output is a set of specialized NICHES objects which allow fine analysis of cell-cell interactions.
#' 
#' @param object Seurat object. The active identity will be used to define populations for connectomic sampling and crossings.
#' @param assay string. Default: "RNA". The assay to run the NICHES transformation on. 
#' @param LR.database string. Default: "fantom5". Currently accepts "fantom5","omnipath", or "custom"
#' @param species string. The species of the object that is being processed. Only required when LR.database = 'fantom5' with species being 'human','mouse','rat', or 'pig', or LR.database = 'omnipath' with species being 'human','mouse', or 'rat'
#' @param min.cells.per.ident integer. Default: 1. A limit on how small (how many cells) a single population can be to participate in connectomic crossings.
#' @param min.cells.per.gene integer. Default: 10. Limits analysis to interactions involving genes expressed above minimum threshold number of cells in the system. 
#' @param meta.data.to.map character vector. Optional. Default: NULL. A vector of metadata names present in the original object which will be carried to the NICHES objects
#' @param position.x string. Optional. Default: NULL. The name of the meta.data column specifying location on the spatial x-axis. Only required for spatial omics data.
#' @param position.y string. Optional. Default: NULL. The name of the meta.data column specifying location on the spatial y-axis. Only required for spatial omics data.
#' @param custom_LR_database data.frame. Optional. Default: NULL. Only required when LR.database = "custom". Each row is a ligand-receptor mechanism where the first column corresponds to the source genes that express the ligands subunits (separated by '_') and the second column corresponds to the receptor genes that express the receptor subunits (separated by '_').
#' @param k integer. Optional. Default: 4. Number of neighbors in a knn graph. Used to compute a mutual nearest neighbor graph based on the spatial coordinates of the spatial transcriptomic datasets.  
#' @param rad.set numeric. Optional. Default: NULL. The radius threshold to define neighbors based on the spatial coordinates of the spatial transcriptomic datasets. Ignored when 'k' is provided.
#' @param blend string. Default: "sum". Choice of linear operator to combine edges in single-cell niche investigations. Defaults to "sum", also accepts "mean".
#' @param CellToCell logical. Default: TRUE. Whether to analyze cell-cell interactions without considering spatial coordinates.
#' @param CellToSystem logical. Default: FALSE. Whether to analyze summed signaling output to total system coming from each cell. Does not consider Euclidean coordinates.
#' @param SystemToCell logical. Default: FALSE. Whether to analyze summed signaling input from total system landing on each cell (cellular microenvironment/niche). Does not consider Euclidean coordinates. 
#' @param CellToCellSpatial logical. Default: FALSE. Whether to analyze cell-cell interactions between Euclidean neighbors. Only applicable in spatial datasets.
#' @param CellToNeighborhood logical. Default: FALSE. Whether to analyze summed signaling output to Euclidean neighbors. Only applicable in spatial datasets.
#' @param NeighborhoodToCell logical. Default: FALSE. Whether to analyze summed signaling input from Euclidean neighbors (cellular microenvironment/niche). Only applicable in spatial datasets.
#'
#' @export


RunNICHES <- function(object,
                        assay="RNA",
                        LR.database="fantom5",
                        species,
                        min.cells.per.ident = NULL,
                        min.cells.per.gene = NULL,
                        meta.data.to.map = NULL,
                        position.x = NULL,
                        position.y = NULL,
                        custom_LR_database = NULL,
                        k = 4,
                        rad.set = NULL,
                        blend = 'sum',
                        CellToCell = T,
                        CellToSystem = F,
                        SystemToCell = F,
                        CellToCellSpatial = F,
                        CellToNeighborhood = F,
                        NeighborhoodToCell = F
                        ){
  # TODO: check the parameter validity here, then register the parameters
  # 1: check the data format of the required parameters and the optional parameters (if(!is.null())): int, character, etc.
  # 2. check some of the dependencies of the parameters
  #   (1) If LR.database is 'custom', then 'custom_LR_database' can't be NULL
  #   (2) If CellToCellSpatial,CellToNeighborhood, or NeighborhoodToCell is T, then rad.set, 'position.x' and 'position.y' can't be null
  # 3. check whether other functions have done similar checks
  # 4. think whether to put stop or warning
  
  # check assay:
  if(is.null(object@assays[[assay]])) stop(paste0("Assay ",assay," is NULL"))
  
  # check lr database, species, custom_LR_database
  if(LR.database %in% c("fantom5","omnipath","custom")){
    if(LR.database == "fantom5"){
      if(!(species %in% c('human','mouse','rat','pig'))) 
        stop(paste0("Unsupported species ",species, " for fantom 5 database, only 'human','mouse','rat',or 'pig' supported."))
    }
    if(LR.database == "omnipath"){
      if(!(species %in% c('human','mouse','rat')))
        stop(paste0("Unsupported species ",species, " for omnipath database, only 'human','mouse',or 'rat' supported."))
    }
    if(LR.database == "custom"){
      # check the format of custom_LR_database
      if(is.null(custom_LR_database)) stop("custom_LR_database is NULL")
      # TODO: check gene names
      message("Custom Ligand Receptor database enabled...")
      message("Checking the format of the custom database...")
      if(!is.data.frame(custom_LR_database)){
        warning("Custom database provided is not in dataframe format.")
        warning("Converting to dataframe format...")
        custom_LR_database <- as.data.frame(custom_LR_database)
      }
      if(ncol(custom_LR_database) < 2) stop("Custom database provided contains less than 2 columns.")
    }
  }else stop('\n LR.receptor argument not recognized. Only accepts "omnipath","fantom5" or "custom".')
  
  if(!is.null(custom_LR_database) & LR.database!="custom"){ 
    warning("custom_LR_database is provided but LR.databse is not specified as 'custom'")
}
  # TODO: check min.cells.per.ident and min.cells.per.gene (only integers)?
  ## MSBR:: THIS CAUSES A BUG -- needs to accept NULL as well for if-logic in prepSeurat
  #min.cells.per.ident <- as.integer(min.cells.per.ident)
  #min.cells.per.gene <- as.integer(min.cells.per.gene)
  
  # check meta.data.to.map
  if(!is.null(meta.data.to.map)){
    # check each meta data is present in the provided seurat object
    for(each_meta in meta.data.to.map){
      if(!each_meta %in% colnames(object@meta.data)) 
        warning(paste0("Metadata: "),each_meta," is not present in the provided Seurat object")
    }
  }
  # check indicators
  # jc: Add organization names to the list
  org_names_indicator <- c(CellToCell,CellToSystem,SystemToCell,CellToCellSpatial,CellToNeighborhood,NeighborhoodToCell)
  org_names_indicator <- sapply(org_names_indicator,function(org){
    org_out <- as.logical(org)
    if(is.na(org_out)) warning("Organization indicator ",org," is not set to TRUE or FALSE.")
    return(org_out)
  })
  names(org_names_indicator) <- c("CellToCell","CellToSystem","SystemToCell","CellToCellSpatial","CellToNeighborhood","NeighborhoodToCell")
  
  
  # If requested, additionally calculate spatially-limited NICHES organizations
  if (org_names_indicator["CellToCellSpatial"] == T | org_names_indicator["CellToNeighborhood"] == T | org_names_indicator["NeighborhoodToCell"] == T){
    
    if (is.null(position.x) | is.null(position.y)){stop("\n Position information not provided. Please specify metadata columns containing x- and y-axis spatial coordinates.")}
    if(!is.null(k)) k <- as.integer(k)
    if(!is.null(rad.set)) rad.set <- as.numeric(rad.set)
  }
  
  if((!is.null(position.x) | !is.null(position.y)) & org_names_indicator["CellToCellSpatial"] == F & org_names_indicator["CellToNeighborhood"] == F & org_names_indicator["NeighborhoodToCell"] == T)
    warning("Spatial positions are provided but the spatial organization functions: 'CellToCellSpatial','CellToNeighborhood', and 'NeighborhoodToCell' are set to FALSE.")
    
  # TODO: check blend with the system indicators, blend is not in other functions?
  if(org_names_indicator["CellToSystem"] == T | org_names_indicator["SystemToCell"] == T)
    if(!blend %in% c("sum","mean")) stop("blend paramter is not recognized: need to be 'sum' or 'mean" )
  
  # Initialize output structure
  output <- list()
  
  
  # jc: move the shared preprocessing steps here to avoid redundancy and reduce the number of parameters to be passed to other functions
  sys.small <- prepSeurat(object,assay,min.cells.per.ident,min.cells.per.gene)
  ground.truth <- lr_load(LR.database,custom_LR_database,species,rownames(sys.small@assays[[assay]]))
  if (org_names_indicator["CellToCellSpatial"] == T | org_names_indicator["CellToNeighborhood"] == T | org_names_indicator["NeighborhoodToCell"] == T){
    ## 1. Move the neighbor graph construction here
    ## 2. Enable a k-nearest-neighbor parameter as an alternative
    edgelist <- compute_edgelist(sys.small,position.x,position.y,k,rad.set)
  }
  
  
  # Calculate NICHES organizations without spatial restrictions
  # jc: only pass the processed data to each function
  if (CellToCell == T){output[[length(output)+1]] <- RunCellToCell(sys.small=sys.small,
                                                                   ground.truth=ground.truth,
                                                                   assay = assay,
                                                                   meta.data.to.map = meta.data.to.map
                                                                   )}
  if (CellToSystem == T){output[[length(output)+1]] <- RunCellToSystem(sys.small=sys.small,
                                                                       ground.truth=ground.truth,
                                                                       assay = assay,
                                                                       meta.data.to.map = meta.data.to.map,
                                                                       blend = blend
                                                                       )}
  if (SystemToCell == T){output[[length(output)+1]] <- RunSystemToCell(sys.small=sys.small,
                                                                       ground.truth=ground.truth,
                                                                       assay = assay,
                                                                       meta.data.to.map = meta.data.to.map,
                                                                       blend = blend
                                                                       )}
  
  
  if (CellToCellSpatial == T){output[[length(output)+1]] <- RunCellToCellSpatial(sys.small=sys.small,
                                                                                 ground.truth=ground.truth,
                                                                                 assay = assay,
                                                                                 meta.data.to.map = meta.data.to.map,
                                                                                 edgelist = edgelist
                                                                                 )} #Spatially-limited Cell-Cell vectors
  if (CellToNeighborhood == T){output[[length(output)+1]] <- RunCellToNeighborhood(sys.small=sys.small,
                                                                                   ground.truth=ground.truth,
                                                                                   assay = assay,
                                                                                   meta.data.to.map = meta.data.to.map,
                                                                                   edgelist = edgelist
                                                                                   )} #Spatially-limited Cell-Neighborhood vectors
  if (NeighborhoodToCell == T){output[[length(output)+1]] <- RunNeighborhoodToCell(sys.small=sys.small,
                                                                                   ground.truth=ground.truth,
                                                                                   assay = assay,
                                                                                   meta.data.to.map = meta.data.to.map,
                                                                                   edgelist = edgelist
                                                                                   )} #Spatially-limited Neighborhood-Cell vectors (niches)

  # jc: Add organization names to the list
  names(output) <- names(org_names_indicator)[org_names_indicator]
  
  # Compile objects for output
  return(output)
}

