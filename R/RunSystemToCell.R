#' RunSystemToCell
#' 
#' Condenses signaling edges landing on each cell within a Seurat object, from all other cells in the system treated as a single source. 
#' Outputs another Seurat object, but where the rows of the matrix are ligand-receptor mechanisms
#' and the columns are each a single receiving cell barcode. The information in the matrix is a sum (or an average, depending on user preference) of
#' all signaling edges landing on that particular cell, from all cells in the system (including from itself.)
#' This transformation allows rapid manipulation and dimensional reduction of how a cell is connected within the system.
#' The default assay of this object is called "SystemToCell" to distinguish it from other Seurat objects.
#' Meta.data slots by default contain "ReceivingType" information, which is the celltypes for each point, 
#' and "ReceivingCell" which is the exact cell barcode present in the original Seurat object.
#' 
#' @param object A Seurat 3.0 object.  The active identity meta.data will be used to define populations for connectomic sampling and crossings.
#' @param LR.database Accepts either 'fantom5' or a custom data.frame with the first column equal to ligands, second column equal to associated receptors.
#' @param species The species of the object that is being processed.  Only required if LR.database = 'fantom5', and allows 'human','mouse','rat', or 'pig'
#' @param assay The assay to run the SystemCell transformation on. Defaults to "RNA."
#' @param min.cells.per.ident Default 1. A limit on how small (how many cells) a single population can be to participate in connectomic crossings.
#' @param blend Choice of linear operator to combine edges. Defaults to "sum", also accepts "mean"
#'
#' @export


RunSystemToCell <- function(object,
                   LR.database = 'fantom5',
                   species,
                   assay = 'RNA',
                   min.cells.per.ident = 1,
                   blend = 'sum',
                   meta.data.to.map = NULL){
  
  # jc: wrapped the preprocessing steps
  sys.small <- prepSeurat(object,assay,min.cells.per.ident)
  
  # jc: Load corresponding ligands and receptors
  ground.truth <- lr_load(LR.database,species,rownames(sys.small@assays[[assay]]))

  ### CREATE MAPPING ###
  
  # Ligand data
  subunit.list <- list() # Builds sending (ligand) data for any number of ligand subunits
  for (s in 1:ncol(ground.truth$source.subunits)){ #For each subunit column...
    subunit.list[[s]] <- matrix(data = NA_real_,nrow = nrow(ground.truth$source.subunits),ncol = ncol(sys.small)) #initialize a mechanism x barcode matrix of all NAs
    colnames(subunit.list[[s]]) <- colnames(sys.small)
    rownames(subunit.list[[s]]) <- rownames(ground.truth$source.subunits)
    non.na.indices <- !is.na(ground.truth$source.subunits[,s]) #Identify rows in the s-th column of the ground truth which are not NA
    subunit.list[[s]][non.na.indices,] <- as.matrix(sys.small@assays[[assay]]@data[ground.truth$source.subunits[non.na.indices,s],])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
  }
  lig.map <- Reduce('*',subunit.list)
  rm(subunit.list)
  
  # Make COMBINED-ACROSS-SYSTEM LIGAND INFO
    if (blend == 'sum'){
    lig.map2 <- Matrix::rowSums(lig.map,dims = 1)
  }
  if (blend == 'mean'){
    lig.map2 <- Matrix::rowMeans(lig.map,dims = 1)
  }
  lig.map2 <- do.call(cbind, replicate(ncol(lig.map), lig.map2, simplify=FALSE))
  
  # Receptor data
  subunit.list <- list() # Builds receiving (receptor) data for any number of receptor subunits
  for (t in 1:ncol(ground.truth$target.subunits)){
    subunit.list[[s]] <- matrix(data = NA_real_,nrow = nrow(ground.truth$target.subunits),ncol = ncol(sys.small)) #initialize a mechanism x barcode matrix of all NAs
    colnames(subunit.list[[s]]) <- colnames(sys.small)
    rownames(subunit.list[[s]]) <- rownames(ground.truth$target.subunits)
    non.na.indices <- !is.na(ground.truth$target.subunits[,s]) #Identify rows in the s-th column of the ground truth which are not NA
    subunit.list[[s]][non.na.indices,] <- as.matrix(sys.small@assays[[assay]]@data[ground.truth$target.subunits[non.na.indices,s],])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
  }
  rec.map <- Reduce('*',subunit.list)
  rm(subunit.list)
  
  # Merged map (can be done with any operator, here is multiplication (RECOMMENDED: preserves zeroes and is quantitative))
  sc.connectome <- lig.map2*rec.map


  # Create the rownames (directed ligands and receptors)
  rownames(sc.connectome) <- paste(rownames(lig.map),rownames(rec.map),sep = '-')
  # Create the column names (directed System-cell)
  colnames(sc.connectome) <- paste("System",colnames(rec.map),sep = '-')
  
  #Use this matrix to create a Seurat object:
  demo <- Seurat::CreateSeuratObject(counts = as.matrix(sc.connectome),assay = 'SystemToCell')
  
  # Add metadata to the Seurat object
  meta.data.to.add <- data.frame(as.character(colnames(rec.map)))
  rownames(meta.data.to.add) <- paste("System",colnames(rec.map),sep = '-')
  demo <- Seurat::AddMetaData(demo,metadata = meta.data.to.add,col.name = 'ReceivingCell')
  demo <- Seurat::AddMetaData(demo,metadata = Seurat::Idents(sys.small),col.name = 'ReceivingType')
  
  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
    # Identify sending and receiving barcodes
    #sending.barcodes <- colnames(do.call(cbind,lig.data)) # Only receiving cell metadata applies for this function
    receiving.barcodes <- colnames(rec.map) 
    # Pull and format sending and receiving metadata
    #sending.metadata <- as.matrix(object@meta.data[,meta.data.to.map][sending.barcodes,])
    receiving.metadata <- as.matrix(object@meta.data[,meta.data.to.map][receiving.barcodes,])
    # Make joint metadata
    #datArray <- abind(sending.metadata,receiving.metadata,along=3)
    #joint.metadata <- as.matrix(apply(datArray,1:2,function(x)paste(x[1],"-",x[2])))
    # Define column names
    #colnames(joint.metadata) <- paste(colnames(sending.metadata),'Joint',sep = '.')
    #colnames(sending.metadata) <- paste(colnames(sending.metadata),'Sending',sep='.')
    #colnames(receiving.metadata) <- paste(colnames(receiving.metadata),'Receiving',sep='.')
    # Compile
    meta.data.to.add.also <- receiving.metadata
    rownames(meta.data.to.add.also) <- paste('System',receiving.barcodes,sep='-')
    # Add additional metadata
    demo <- Seurat::AddMetaData(demo,metadata = as.data.frame(meta.data.to.add.also))
  }
  
  # How many vectors were captured by this sampling?
  
  message(paste("\n",length(unique(demo$ReceivingCell)),'System-To-Cell edges were computed, across',length(unique(demo$ReceivingType)),'cell types'))
  return(demo)
}

