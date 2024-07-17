#' RunNeighborhoodToCell
#' 
#' @param sys.small A filtered Seurat object. The active identity will be used to define populations for connectomic sampling and crossings.
#' @param ground.truth Ground truth signaling mechanisms present in sys.small.
#' @param assay The assay to run the SCC transformation on. Defaults to "RNA."
#' @param meta.data.to.map A character vector of metadata names present in the original object which will be carried to the NICHES objects
#' @param blend Choice of linear operator to combine edges. Defaults to "mean", also accepts "sum"
#' @param edgelist data.frame. Each row is an directional edge between two spatially connected cells
#' @param output_format string. Choice of the output format. "seurat" will output a list of seurat objects, "raw" will output a list of lists with raw interaction matrix and compiled metadata
#'
#' @export

RunNeighborhoodToCell <- function(sys.small,
                                  ground.truth,
                                  assay,
                                  meta.data.to.map,
                                  blend="mean",
                                  edgelist,
                                  output_format
                                  ){
  
  # Make ligand matrix
  
  #lig.data <- sys.small@assays[[assay]]@data[ligands,edgelist$from]
  
  subunit.list <- list() # Builds sending (ligand) data for any number of ligand subunits
  for (s in 1:ncol(ground.truth$source.subunits)){ #For each subunit column...
    subunit.list[[s]] <- matrix(data = 1,nrow = nrow(ground.truth$source.subunits),ncol = ncol(getSeuratAssay(sys.small,assay,"data")[,edgelist$from])) #initialize a mechanism x barcode matrix of all NAs
    colnames(subunit.list[[s]]) <- colnames(getSeuratAssay(sys.small,assay,"data")[,edgelist$from])
    rownames(subunit.list[[s]]) <- rownames(ground.truth$source.subunits)
    non.na.indices <- !is.na(ground.truth$source.subunits[,s]) #Identify rows in the s-th column of the ground truth which are not NA
    subunit.list[[s]][non.na.indices,] <- as.matrix(getSeuratAssay(sys.small,assay,"data")[ground.truth$source.subunits[non.na.indices,s],edgelist$from])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
  }
  lig.data <- Reduce('*',subunit.list)
  rm(subunit.list)
  
  # Make receptor matrix
  
  #rec.data <- sys.small@assays[[assay]]@data[receptors,edgelist$to]
  
  subunit.list <- list() # Builds receiving (receptor) data for any number of receptor subunits
  for (t in 1:ncol(ground.truth$target.subunits)){
    subunit.list[[t]] <- matrix(data = 1,nrow = nrow(ground.truth$target.subunits),ncol = ncol(getSeuratAssay(sys.small,assay,"data")[,edgelist$to])) #initialize a mechanism x barcode matrix of all NAs
    colnames(subunit.list[[t]]) <- colnames(getSeuratAssay(sys.small,assay,"data")[,edgelist$to])
    rownames(subunit.list[[t]]) <- rownames(ground.truth$target.subunits)
    non.na.indices <- !is.na(ground.truth$target.subunits[,t]) #Identify rows in the t-th column of the ground truth which are not NA
    subunit.list[[t]][non.na.indices,] <- as.matrix(getSeuratAssay(sys.small,assay,"data")[ground.truth$target.subunits[non.na.indices,t],edgelist$to])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
  }
  rec.data <- Reduce('*',subunit.list)
  rm(subunit.list)
  
  
  # Make SCC matrix
  scc <- lig.data*rec.data
  rownames(scc) <- paste(rownames(lig.data),rownames(rec.data),sep = '—')
  
  # Condense by column name
  colnames(scc) <- colnames(rec.data) # Make colnames equal to receiving cell
  scc <- as.matrix(scc)
  if(blend == "sum") scc <- t(rowsum(t(scc), colnames(scc)))
  else if(blend == "mean") scc <- sapply(unique(colnames(scc)), function(receiving_name) 
    rowMeans(scc[,colnames(scc)== receiving_name,drop=FALSE], na.rm=TRUE) )
    
  # Label columns properly
  barcodes <- colnames(scc)
  colnames(scc) <- paste('Neighborhood',colnames(scc),sep = '—')
  
  # Use this matrix to create a Seurat object:
  demo <- Seurat::CreateSeuratObject(counts = as.matrix(scc),assay = 'NeighborhoodToCell')
  # JC: Seurat V5 will not create data slot automatically, the following step is to manually add this slot
  if(SeuratObject::Version(demo) >= "5.0.0"){
    demo <- NormalizeData(demo,assay = "NeighborhoodToCell")  # Seura Object need to be >= 5.0.1
    demo@assays$NeighborhoodToCell@layers$data <- demo@assays$NeighborhoodToCell@layers$counts # Seura Object need to be >= 5.0.1
    
  }
  
  # Add metadata based on ident slot
  demo <- Seurat::AddMetaData(demo,metadata = barcodes,col.name = 'ReceivingCell')
  # bug fix: add the Neighborhood - prefix
  receiving_type.meta <- data.frame(Seurat::Idents(sys.small)[barcodes])
  rownames(receiving_type.meta) <- paste("Neighborhood",rownames(receiving_type.meta),sep = '—')
  
  demo <- Seurat::AddMetaData(demo,metadata = receiving_type.meta,col.name = 'ReceivingType')
  
  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
    # Identify sending and receiving barcodes
    #sending.barcodes <- barcodes 
    receiving.barcodes <- barcodes # Only receiving cell metadata applies for this function
    # Pull and format sending and receiving metadata
    #sending.metadata <- as.matrix(object@meta.data[,meta.data.to.map][sending.barcodes,])
    # jc: possible bug, change object to sys.small
    receiving.metadata <- as.matrix(sys.small@meta.data[,meta.data.to.map,drop=FALSE][receiving.barcodes,])
    # Make joint metadata
    #datArray <- abind(sending.metadata,receiving.metadata,along=3)
    #joint.metadata <- as.matrix(apply(datArray,1:2,function(x)paste(x[1],"-",x[2])))
    # Define column names
    #colnames(joint.metadata) <- paste(colnames(sending.metadata),'Joint',sep = '.')
    #colnames(sending.metadata) <- paste(colnames(sending.metadata),'Sending',sep='.')
    #colnames(receiving.metadata) <- paste(colnames(receiving.metadata),'Receiving',sep='.')
    # Compile
    meta.data.to.add.also <- receiving.metadata
    rownames(meta.data.to.add.also) <- paste('Neighborhood',receiving.barcodes,sep='—')
    # Add additional metadata
    demo <- Seurat::AddMetaData(demo,metadata = as.data.frame(meta.data.to.add.also))
  }
  # Set initial identity
  Seurat::Idents(demo) <- demo$ReceivingType
  # How many vectors were captured by this sampling?
  message(paste("\n",length(unique(demo$ReceivingCell)),'Neighborhood-To-Cell edges were computed, across',length(unique(demo$ReceivingType)),'cell types'))
  
  if(output_format == "seurat") return(demo)
  else{
    output_list <- vector(mode = "list",length=2)
    names(output_list) <- c("NeighborhoodToCellMatrix","metadata")
    output_list[["NeighborhoodToCellMatrix"]] <- getSeuratAssay(demo,"NeighborhoodToCell","counts")
    output_list[["metadata"]] <- demo@meta.data
    return(output_list)
  }
  
}

