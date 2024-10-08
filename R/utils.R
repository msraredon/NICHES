#' Wrapper function for GetAssayData based on the version of the SeuratObject
#'
#' @param object input Seurat object 
#' @param assay which assay to use in object
#' @param slot "data","counts",or "scale.data"
#'
#' @return the specified assay data
#' @export
getSeuratAssay <- function(object,assay,slot){
  if(SeuratObject::Version(object) >= "5.0.0") return(Seurat::GetAssayData(object,assay=assay,layer=slot))
  else return(Seurat::GetAssayData(object,assay=assay,slot=slot))
}

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
  filtered.obj <- object
  
  
  # Limit object to cell populations larger than requested minimum
  if (!is.null(min.cells.per.ident)){
    message(paste("\n",'Subsetting to populations with greater than',min.cells.per.ident,'cells'))
    idents.include <- names(table(Seurat::Idents(filtered.obj)))[table(Seurat::Idents(filtered.obj)) > min.cells.per.ident]
    filtered.obj <- subset(filtered.obj,idents = idents.include)
  }
  
  # Limit analysis to interaction involving genes expressed above minimum threshold number of cells in the system
  if (!is.null(min.cells.per.gene)){
    message(paste("\n",'Subsetting to genes expressed in greater than',min.cells.per.gene,'cells'))
    
    # jc: return an error if the count matrix in the given assay is empty
    if(!length(getSeuratAssay(filtered.obj,assay,"counts"))) stop("Unable to subset: the count matrix in the given assay is empty.")
    
    cells.per.gene <- data.frame(non.zero.cells = Matrix::rowSums(getSeuratAssay(filtered.obj,assay,"counts")>0))
    GOI <- subset(cells.per.gene,non.zero.cells > min.cells.per.gene)
    filtered.obj <- filtered.obj[rownames(GOI),]
  }
  
  num.cells <- ncol(filtered.obj)
  message(paste("\n",num.cells,'distinct cells from',
                length(names(table(Seurat::Idents(filtered.obj)))),'celltypes to be analyzed'))
  
  return(filtered.obj)
}


#' load the Ligands and Receptors that corresponds to the row names of the input seurat
#'
#' @param LR.database database to use
#' @param species species in the database (fantom5 or omnipath) to use
#' @param input_rownames the genes names to query
#' @param custom_LR_database data.frame. Each row is a ligand-receptor mechanism where the first column corresponds to the source genes that express the ligands subunits (separated by '_') and the second column corresponds to the receptor genes that express the receptor subunits (separated by '_'). 
#'
#' @return a list of ligands and receptors
#' @export

lr_load <- function(LR.database,custom_LR_database=NULL,species,input_rownames){
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

# better: check_celltypes: check whether the idents are cell types, yes to return the unique cell types, no to return an error
#' return cell types
#'
#' @param seurat_object Input seurat object
#'
#' @return unique cell types
#' @export

return_celltypes <- function(seurat_object){
  message('\n For sampling purposes, please make sure that the active Identity of the input seurat object corresponds to cell types')
  return(names(table(Seurat::Idents(seurat_object)))) # ms: other steps are depdent on the exact format of this output
}


#' Compute an edgelist based on the spatial coordinates
#'
#' @param filtered.obj A filtered Seurat object. The active identity will be used to define populations for connectomic sampling and crossings.
#' @param position.x string. Optional. Default: NULL. The name of the meta.data column specifying location on the spatial x-axis. Only required for spatial omics data.
#' @param position.y string. Optional. Default: NULL. The name of the meta.data column specifying location on the spatial y-axis. Only required for spatial omics data.
#' @param k integer. Optional. Default: 4. Number of neighbors in a knn graph. Used to compute a mutual nearest neighbor graph based on the spatial coordinates of the spatial transcriptomic datasets.  
#' @param rad.set numeric. Optional. Default: NULL. The radius threshold to define neighbors based on the spatial coordinates of the spatial transcriptomic datasets. Ignored when 'k' is provided.
#'
#' @export
#'
compute_edgelist <- function(filtered.obj,
                             position.x,
                             position.y,
                             k=4,
                             rad.set=NULL
                             ){
  
  ### CREATE MAPPING ###
  
  # Create adjacency matrix
  # Adapted from :: https://stackoverflow.com/questions/16075232/how-to-create-adjacency-matrix-from-grid-coordinates-in-r
  # Setup numbering and labeling
  # jc: possible bug, change object to filtered.obj
  df <- data.frame(x = filtered.obj[[position.x]], y = filtered.obj[[position.y]])
  df$barcode <- rownames(df)
  df$x <- as.character(df$x)
  df$y <- as.character(df$y)
  df$x <- as.numeric(df$x)
  df$y <- as.numeric(df$y)
  df <- df[,c('x','y')]
  
  # Compute the euclidean distance matrix
  distance_mat <- apply(df, 1, function(pt)
    (sqrt(abs(pt["x"] - df$x)^2 + abs(pt["y"] - df$y)^2))
  )
  # keep the cell names
  rownames(distance_mat) <- colnames(distance_mat)
  
  if(!is.null(k)){
    message("Compute edgelist based on mutual nearest neighbors.")
    if(!is.null(rad.set))warning("'k' is not NULL. Parameter 'rad.set' will be ignored.")
    # neighbor_mat: the indices of k nearest neighbors  
    neighbor_mat <- apply(distance_mat,1,function(dis_vec){
      order(dis_vec)[1:(k+1)]
    })
    
    # k nearest neighbor adjacency matrix
    adj_mat <- matrix(data=0,nrow=nrow(distance_mat),ncol = ncol(distance_mat))
    rownames(adj_mat) <- rownames(distance_mat)
    colnames(adj_mat) <- colnames(distance_mat)
    for(cell in colnames(neighbor_mat)) adj_mat[cell,neighbor_mat[,cell]] <- 1
    # mutual nearest neighbor adjacency matrix
    adj_mat_final <- 1*(adj_mat & t(adj_mat))
    
  }
  else if(!is.null(rad.set)){
    message("\n Compute edgelist based on spatial radius threshold.")
    # Alternatively the adjacency matrix is defined by the absolute radius
    adj_mat_final <- 1*(distance_mat <= rad.set)
  }
  else{
    stop("Both k and rad.set are NULL.")
  }
  
  # Convert adj matrix to edgelist
  edgelist <- igraph::graph.adjacency(adj_mat_final)
  edgelist <- igraph::get.data.frame(edgelist)
  
  return(edgelist)
  
}