
#' load the Ligands and Receptors that corresponds to the row names of the input seurat
#'
#' @param LR.database database to use
#' @param species species in the database (fantom5 or omnipath) to use
#' @param input_rownames the genes names to query
#' @param custom_LR_database data.frame. Each row is a ligand-receptor mechanism where the first column corresponds to the source genes that express the ligands subunits (separated by '_') and the second column corresponds to the receptor genes that express the receptor subunits (separated by '_'). 
#'
#' @return a list of ligands and receptors
#' @export

LoadGroundTruth <- function(LR.database,custom_LR_database=NULL,species,input_rownames){
  if (LR.database == 'omnipath'){
    ground.truth <- LoadOmniPath(species = species)
  }else if (LR.database == 'fantom5'){
    ground.truth <- LoadFantom5(species = species)
  }else if (LR.database == 'custom'){
    ground.truth <- LoadCustom(custom_LR_database)
  }else {
    stop('\n LR.receptor argument not recognized. Only accepts "omnipath","fantom5" or "custom".')
  }
  
  # Subset to only map against mechanisms present in the base Seurat object
  ground.truth <- FilterGroundTruth(ground.truth,
                                    input_rownames = input_rownames)
  
  return(ground.truth)
}