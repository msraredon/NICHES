
#' Seurat input preprocessing
#'
#' @param object input Seurat object 
#' @param assay which assay to use in object
#' @param min.cells.per.ident minimum cells per ident
#' @param min.cells.per.gene minimum cells which must express a gene for it to be considered in connectomics
#'
#' @return a processed seurat object
#' @export

prepSeurat <- function(object,assay,min.cells.per.ident,min.cells.per.gene){
  # Set default assay (not necessary, but just in case)
  Seurat::DefaultAssay(object) <- assay
  
  # Stash object
  sys.small <- object
  
  # Limit object to cell populations larger than requested minimum
  if (!is.null(min.cells.per.ident)){
    message(paste("\n",'Subsetting to populations with greater than',min.cells.per.ident,'cells'))
    idents.include <- names(table(Seurat::Idents(sys.small)))[table(Seurat::Idents(sys.small)) > min.cells.per.ident]
    sys.small <- subset(sys.small,idents = idents.include)
  }
  
  # Limit analysis to interaction involving genes expressed above minimum threshold number of cells in the system
  if (!is.null(min.cells.per.gene)){
    message(paste("\n",'Subsetting to genes expressed in greater than',min.cells.per.gene,'cells'))
    cells.per.gene <- data.frame(non.zero.cells = Matrix::rowSums(sys.small@assays[[assay]]@counts>0))
    GOI <- subset(cells.per.gene,non.zero.cells > min.cells.per.gene)
    sys.small <- sys.small[rownames(GOI),]
  }
  
  num.cells <- ncol(sys.small)
  message(paste("\n",num.cells,'distinct cells from',
                length(names(table(Seurat::Idents(sys.small)))),'celltypes to be analyzed, across',nrow(GOI),'signaling genes'))
  
  return(sys.small)
}


#' load the Ligands and Receptors that corresponds to the row names of the input seurat
#'
#' @param LR.database database to use
#' @param species 
#' @param input_rownames the genes names to query
#'
#' @return a list of ligands and receptors
#' @export

lr_load <- function(LR.database,species,input_rownames){
  if (LR.database == 'omnipath'){
    ground.truth <- LoadOmniPath(species = species)
  }else if (LR.database == 'fantom5'){
    ground.truth <- LoadFantom5(species = species)
  }else if (LR.database == 'custom'){
    
  }else {
    stop('\n LR.receptor argument not recognized. Only accepts "omnipath","fantom5" or "custom".')
  }
  
  # Subset to only map against mechanisms present in the base Seurat object
  ground.truth <- FilterGroundTruth(ground.truth,
                                    input_rownames = input_rownames)
  
  return(ground.truth)
}

# better: check_celltypes: check whether the idents are cell types, yes to return the unique cell types, no to return an error
#' return cell types
#'
#' @param seurat_object 
#'
#' @return unique cell types
#' @export

return_celltypes <- function(seurat_object){
  message('\n For sampling purposes, please make sure that the active Identity of the input seurat object corresponds to cell types')
  return(names(table(Seurat::Idents(seurat_object)))) # ms: other steps are depdent on the exact format of this output
}
