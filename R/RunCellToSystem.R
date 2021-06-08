#' RunCellToSystem
#' 
#' Condenses all signaling edges coming from each cell within a Seurat object connecting to every other cell in the system, treated as a single sink.
#' Outputs another Seurat object, but where the rows of the matrix are ligand-receptor mechanisms
#' and the columns are each a single sending cell barcode. The information in the matrix is a sum (or an average, depending on user preference) of
#' all signaling edges coming from that particular cell, to every other cell in the system (including to itself.)
#' This transformation allows rapid manipulation and dimensional reduction of how a cell is connected within the system.
#' The default assay of this object is called "CellToSystem" to distinguish it from other Seurat objects.
#' Meta.data slots by default contain "SendingType" information, which is the celltypes for each point, 
#' and "SendingCell" which is the exact cell barcode present in the original Seurat object.
#' 
#' @param object A Seurat 3.0 object.  The active identity meta.data will be used to define populations for connectomic sampling and crossings.
#' @param LR.database Accepts either 'fantom5' or a custom data.frame with the first column equal to ligands, second column equal to associated receptors.
#' @param species The species of the object that is being processed.  Only required if LR.database = 'fantom5', and allows 'human','mouse','rat', or 'pig'
#' @param assay The assay to run the CellToSystem transformation on. Defaults to "RNA."
#' @param min.cells.per.ident Default 1. A limit on how small (how many cells) a single population can be to participate in connectomic crossings.
#' @param blend Choice of linear operator to combine edges. Defaults to "sum", also accepts "mean"
#'
#' @export


RunCellToSystem <- function(object,
                      LR.database,
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
  
  # Receptor data
  subunit.list <- list() # Builds receiving (receptor) data for any number of receptor subunits
  for (t in 1:ncol(ground.truth$target.subunits)){
    subunit.list[[t]] <- matrix(data = NA_real_,nrow = nrow(ground.truth$target.subunits),ncol = ncol(sys.small)) #initialize a mechanism x barcode matrix of all NAs
    colnames(subunit.list[[t]]) <- colnames(sys.small)
    rownames(subunit.list[[t]]) <- rownames(ground.truth$target.subunits)
    non.na.indices <- !is.na(ground.truth$target.subunits[,t]) #Identify rows in the s-th column of the ground truth which are not NA
    subunit.list[[t]][non.na.indices,] <- as.matrix(sys.small@assays[[assay]]@data[ground.truth$target.subunits[non.na.indices,t],])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
  }
  rec.map <- Reduce('*',subunit.list)
  rm(subunit.list)
  
  
  # Make COMBINED-ACROSS-SYSTEM RECEPTOR INFO
  if (blend == 'sum'){
    rec.map2 <- Matrix::rowSums(rec.map,dims = 1)
  }
  if (blend == 'mean'){
    rec.map2 <- Matrix::rowMeans(rec.map,dims = 1)
  }
  rec.map2 <- do.call(cbind, replicate(ncol(rec.map), rec.map2, simplify=FALSE))
  
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
  
  
  # Merged map (can be done with any operator, here is multiplication (RECOMMENDED: preserves zeroes and is quantitative))
  sc.connectome <- lig.map*rec.map2
  
  
  # Create the rownames (directed ligands and receptors)
  rownames(sc.connectome) <- paste(rownames(lig.map),rownames(rec.map),sep = '-')
  # Create the column names (directed cell-System)
  colnames(sc.connectome) <- paste(colnames(lig.map),"System",sep = '-')
  
  # Use this matrix to create a Seurat object:
  demo <- Seurat::CreateSeuratObject(counts = as.matrix(sc.connectome),assay = 'CellToSystem')
  
  # Add metadata to the Seurat object
  meta.data.to.add <- data.frame(as.character(colnames(lig.map)))
  rownames(meta.data.to.add) <- paste(colnames(lig.map),"System",sep = '-')
  demo <- Seurat::AddMetaData(demo,metadata = meta.data.to.add,col.name = 'SendingCell')
  demo <- Seurat::AddMetaData(demo,metadata = Seurat::Idents(sys.small),col.name = 'SendingType')
  
  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
    # Identify sending and receiving barcodes
    sending.barcodes <- colnames(lig.map) # Only sending cell metadata applies for this function
    #receiving.barcodes <- colnames(rec.map) 
    # Pull and format sending and receiving metadata
    sending.metadata <- as.matrix(object@meta.data[,meta.data.to.map][sending.barcodes,])
    #receiving.metadata <- as.matrix(object@meta.data[,meta.data.to.map][receiving.barcodes,])
    # Make joint metadata
    #datArray <- abind(sending.metadata,receiving.metadata,along=3)
    #joint.metadata <- as.matrix(apply(datArray,1:2,function(x)paste(x[1],"-",x[2])))
    # Define column names
    #colnames(joint.metadata) <- paste(colnames(sending.metadata),'Joint',sep = '.')
    #colnames(sending.metadata) <- paste(colnames(sending.metadata),'Sending',sep='.')
    #colnames(receiving.metadata) <- paste(colnames(receiving.metadata),'Receiving',sep='.')
    # Compile
    meta.data.to.add.also <- sending.metadata
    rownames(meta.data.to.add.also) <- paste(sending.barcodes,'System',sep='-')
    # Add additional metadata
    demo <- Seurat::AddMetaData(demo,metadata = as.data.frame(meta.data.to.add.also))
  }
  

  # How many vectors were captured by this sampling?
  
  message(paste("\n",length(unique(demo$SendingCell)),'Cell-To-System edges were computed, across',length(unique(demo$SendingType)),'cell types'))
  return(demo)
}

