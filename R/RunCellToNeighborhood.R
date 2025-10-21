#' RunCellToNeighborhood
#'

#' @param node.object A Seurat object containing cells-as-barcode data. The active identity will be used to define populations for connectomic sampling and crossings.
#' @param ground.truth Ground truth signaling mechanisms to be queried.
#' @param assay The assay to run the SCC transformation on. Defaults to "RNA."
#' @param meta.data.to.map A character vector of metadata names present in the original object which will be carried to the NICHES objects
#' @param blend Choice of linear operator to combine edges. Defaults to "mean", also accepts "sum"
#' @param edgelist data.frame. Each row is an directional edge between two spatially connected cells 
#' @param output_format string. Choice of the output format. "seurat" will output a list of seurat objects, "raw" will output a list of lists with raw interaction matrix and compiled metadata
#'
#' @export


RunCellToNeighborhood <- function(node.object,
                                  ground.truth,
                                  assay,
                                  meta.data.to.map,
                                  blend="mean",
                                  edgelist,
                                  output_format
                                  ){

  # Make ligand matrix
  

  #lig.data <- node.object@assays[[assay]]@data[ligands,edgelist$from]
  
  subunit.list <- list() # Builds sending (ligand) data for any number of ligand subunits
  for (s in 1:ncol(ground.truth$source.subunits)){ #For each subunit column...
    subunit.list[[s]] <- matrix(data = 1,nrow = nrow(ground.truth$source.subunits),ncol = ncol(GetSeuratAssay(node.object,assay,"data")[,edgelist$from])) #initialize a mechanism x barcode matrix of all NAs
    colnames(subunit.list[[s]]) <- colnames(GetSeuratAssay(node.object,assay,"data")[,edgelist$from])
    rownames(subunit.list[[s]]) <- rownames(ground.truth$source.subunits)
    non.na.indices <- !is.na(ground.truth$source.subunits[,s]) #Identify rows in the s-th column of the ground truth which are not NA
    subunit.list[[s]][non.na.indices,] <- as.matrix(GetSeuratAssay(node.object,assay,"data")[ground.truth$source.subunits[non.na.indices,s],edgelist$from])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices

  }
  lig.data <- Reduce('*',subunit.list)
  rm(subunit.list)
  
  # Make receptor matrix
  

  #rec.data <- node.object@assays[[assay]]@data[receptors,edgelist$to]
  
  subunit.list <- list() # Builds receiving (receptor) data for any number of receptor subunits
  for (t in 1:ncol(ground.truth$target.subunits)){
    subunit.list[[t]] <- matrix(data = 1,nrow = nrow(ground.truth$target.subunits),ncol = ncol(GetSeuratAssay(node.object,assay,"data")[,edgelist$to])) #initialize a mechanism x barcode matrix of all NAs
    colnames(subunit.list[[t]]) <- colnames(GetSeuratAssay(node.object,assay,"data")[,edgelist$to])
    rownames(subunit.list[[t]]) <- rownames(ground.truth$target.subunits)
    non.na.indices <- !is.na(ground.truth$target.subunits[,t]) #Identify rows in the t-th column of the ground truth which are not NA
    subunit.list[[t]][non.na.indices,] <- as.matrix(GetSeuratAssay(node.object,assay,"data")[ground.truth$target.subunits[non.na.indices,t],edgelist$to])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices

  }
  rec.data <- Reduce('*',subunit.list)
  rm(subunit.list)
  
  
  # Make SCC matrix
  scc <- lig.data*rec.data
  rownames(scc) <- paste(rownames(lig.data),rownames(rec.data),sep = '—')

  # Condense by column name
  colnames(scc) <- colnames(lig.data) # Make colnames equal to sending cell
  scc <- as.matrix(scc)
  if(blend == "sum") scc <- t(rowsum(t(scc), colnames(scc)))
  else if(blend == "mean")  scc <- sapply(unique(colnames(scc)), function(sending_name) 
    rowMeans(scc[,colnames(scc)== sending_name,drop=FALSE], na.rm=TRUE) )
  
  # Label columns properly
  barcodes <- colnames(scc)
  colnames(scc) <- paste(barcodes,'Neighborhood',sep = '—')

  # Use this matrix to create a Seurat object:
  demo <- Seurat::CreateSeuratObject(counts = as.matrix(scc),assay = 'CellToNeighborhood')
  # JC: Seurat V5 will not create data slot automatically, the following step is to manually add this slot
  if(SeuratObject::Version(demo) >= "5.0.0"){
    demo <- Seurat::NormalizeData(demo,assay = "CellToNeighborhood")  # Seura Object need to be >= 5.0.1
    demo@assays$CellToNeighborhood@layers$data <- demo@assays$CellToNeighborhood@layers$counts # Seura Object need to be >= 5.0.1
    
  }
  
  # Add metadata based on ident slot
  # bug fix: add the Neighborhood - prefix

  sending_type.meta <- data.frame(Seurat::Idents(node.object)[barcodes])
  rownames(sending_type.meta) <- paste(rownames(sending_type.meta),"Neighborhood",sep = '—')

  
  demo <- Seurat::AddMetaData(demo,metadata = sending_type.meta,col.name = c("SendingCell","SendingType"))

  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
    # Identify sending and receiving barcodes
    sending.barcodes <- barcodes # Only sending cell metadata applies for this function
    #receiving.barcodes <- colnames(rec.map)
    # Pull and format sending and receiving metadata

    # jc: possible bug, change object to node.object
    sending.metadata <- as.matrix(node.object@meta.data[,meta.data.to.map,drop=FALSE][sending.barcodes,])

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
    rownames(meta.data.to.add.also) <- paste(sending.barcodes,'Neighborhood',sep='—')
    # Add additional metadata
    demo <- Seurat::AddMetaData(demo,metadata = as.data.frame(meta.data.to.add.also))
  }
  # Set initial identity
  Seurat::Idents(demo) <- demo$SendingType
  # How many vectors were captured by this sampling?
  message(paste("\n",length(unique(demo$SendingCell)),'Cell-To-Neighborhood edges were computed, across',length(unique(demo$SendingType)),'cell types'))

  if(output_format == "seurat") return(demo)
  else{
    output_list <- vector(mode = "list",length=2)
    names(output_list) <- c("CellToNeighborhoodMatrix","metadata")
    output_list[["CellToNeighborhoodMatrix"]] <- GetSeuratAssay(demo,"CellToNeighborhood","counts")
    output_list[["metadata"]] <- demo@meta.data
    return(output_list)
  }
  
}
