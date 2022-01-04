#' RunCellToCellSpatial
#'
#' @param sys.small A filtered Seurat object. The active identity will be used to define populations for connectomic sampling and crossings.
#' @param ground.truth Ground truth signaling mechanisms present in sys.small.
#' @param assay The assay to run the SCC transformation on. Defaults to "RNA."
#' @param meta.data.to.map A character vector of metadata names present in the original object which will be carried to the NICHES objects
#' @param edgelist data.frame. Each row is an directional edge between two spatially connected cells
#' 
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export

RunCellToCellSpatial <- function(sys.small,
                                 ground.truth,
                                 assay,
                                 meta.data.to.map,
                                 edgelist
                                 ){

  
    # Make ligand matrix

    #lig.data <- sys.small@assays[[assay]]@data[ligands,edgelist$from]

    subunit.list <- list() # Builds sending (ligand) data for any number of ligand subunits
    for (s in 1:ncol(ground.truth$source.subunits)){ #For each subunit column...
      subunit.list[[s]] <- matrix(data = 1,nrow = nrow(ground.truth$source.subunits),ncol = ncol(sys.small@assays[[assay]]@data[,edgelist$from])) #initialize a mechanism x barcode matrix of all NAs
      colnames(subunit.list[[s]]) <- colnames(sys.small@assays[[assay]]@data[,edgelist$from])
      rownames(subunit.list[[s]]) <- rownames(ground.truth$source.subunits)
      non.na.indices <- !is.na(ground.truth$source.subunits[,s]) #Identify rows in the s-th column of the ground truth which are not NA
      subunit.list[[s]][non.na.indices,] <- as.matrix(sys.small@assays[[assay]]@data[ground.truth$source.subunits[non.na.indices,s],edgelist$from])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
    }
    lig.data <- Reduce('*',subunit.list)
    rm(subunit.list)
    
  # Make receptor matrix

    #rec.data <- sys.small@assays[[assay]]@data[receptors,edgelist$to]
    
    subunit.list <- list() # Builds receiving (receptor) data for any number of receptor subunits
    for (t in 1:ncol(ground.truth$target.subunits)){
      subunit.list[[t]] <- matrix(data = 1,nrow = nrow(ground.truth$target.subunits),ncol = ncol(sys.small@assays[[assay]]@data[,edgelist$to])) #initialize a mechanism x barcode matrix of all NAs
      colnames(subunit.list[[t]]) <- colnames(sys.small@assays[[assay]]@data[,edgelist$to])
      rownames(subunit.list[[t]]) <- rownames(ground.truth$target.subunits)
      non.na.indices <- !is.na(ground.truth$target.subunits[,t]) #Identify rows in the t-th column of the ground truth which are not NA
      subunit.list[[t]][non.na.indices,] <- as.matrix(sys.small@assays[[assay]]@data[ground.truth$target.subunits[non.na.indices,t],edgelist$to])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
    }
    rec.data <- Reduce('*',subunit.list)
    rm(subunit.list)

  # Make SCC matrix
  scc <- lig.data*rec.data
  rownames(scc) <- paste(rownames(lig.data),rownames(rec.data),sep = '—')
  colnames(scc) <- paste(colnames(lig.data),colnames(rec.data),sep = '—')
  sending.cell.idents <- as.character(Seurat::Idents(sys.small)[colnames(lig.data)])
  receiving.cell.idents <- as.character(Seurat::Idents(sys.small)[colnames(rec.data)])
  dim(scc)

  # Use this matrix to create a Seurat object:
  demo <- Seurat::CreateSeuratObject(counts = as.matrix(scc),assay = 'CellToCellSpatial')

  # Add key metadata

  meta.data.to.add <- data.frame(SendingType = sending.cell.idents,
                                 ReceivingType = receiving.cell.idents)

  rownames(meta.data.to.add) <- colnames(scc)
  meta.data.to.add$VectorType <- paste(meta.data.to.add$SendingType,
                                       meta.data.to.add$ReceivingType,
                                       sep = '—')

  #Add metadata to the Seurat object
  demo <- Seurat::AddMetaData(demo,metadata = meta.data.to.add)

  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
    # Identify sending and receiving barcodes
    sending.barcodes <- colnames(lig.data)
    receiving.barcodes <- colnames(rec.data)
    # Pull and format sending and receiving metadata
    # jc: possible bug, change object to sys.small
    sending.metadata <- as.matrix(sys.small@meta.data[,meta.data.to.map][sending.barcodes,])
    receiving.metadata <- as.matrix(sys.small@meta.data[,meta.data.to.map][receiving.barcodes,])
    # Make joint metadata
    datArray <- abind::abind(sending.metadata,receiving.metadata,along=3)
    joint.metadata <- as.matrix(apply(datArray,1:2,function(x)paste(x[1],"-",x[2])))
    # Define column names
    colnames(joint.metadata) <- paste(colnames(sending.metadata),'Joint',sep = '.')
    colnames(sending.metadata) <- paste(colnames(sending.metadata),'Sending',sep='.')
    colnames(receiving.metadata) <- paste(colnames(receiving.metadata),'Receiving',sep='.')
    # Compile
    meta.data.to.add.also <- cbind(sending.metadata,receiving.metadata,joint.metadata)
    rownames(meta.data.to.add.also) <- paste(sending.barcodes,receiving.barcodes,sep='—')
    # Add additional metadata
    demo <- Seurat::AddMetaData(demo,metadata = as.data.frame(meta.data.to.add.also))
  }
  # Set initial identity
  Seurat::Idents(demo) <- demo$VectorType
  
  # How many vectors were captured by this sampling?
  message(paste("\n",length(unique(demo$VectorType)),'distinct VectorTypes were computed, out of',length(table(Seurat::Idents(sys.small)))^2,'total possible'))

  return(demo)
}
