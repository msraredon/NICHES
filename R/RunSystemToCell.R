#' RunSystemToCell
#' 
#' Condenses signaling edges landing on each cell within a Seurat object, from all other cells in the system treated as a single source. 
#' Outputs another Seurat object, but where the rows of the matrix are ligand-receptor mechanisms
#' and the columns are each a single receiving cell barcode. The information in the matrix is an average (or a sum, depending on user preference) of
#' all signaling edges landing on that particular cell, from all cells in the system (including from itself.)
#' This transformation allows rapid manipulation and dimensional reduction of how a cell is connected within the system.
#' The default assay of this object is called "SystemToCell" to distinguish it from other Seurat objects.
#' Meta.data slots by default contain "ReceivingType" information, which is the celltypes for each point, 
#' and "ReceivingCell" which is the exact cell barcode present in the original Seurat object.
#' 

#' @param node.object A filtered Seurat object. The active identity will be used to define populations for connectomic sampling and crossings. 
#' @param ground.truth Ground truth signaling mechanisms present in node.object.
#' @param assay The assay to run the SystemToCell transformation on. Defaults to "RNA."
#' @param blend Choice of linear operator to combine edges. Defaults to "mean", also accepts "sum", "mean.adj"
#' @param meta.data.to.map A character vector of metadata names present in the original object which will be carried to the NICHES objects
#' @param output_format string. Choice of the output format. "seurat" will output a list of seurat objects, "raw" will output a list of lists with raw interaction matrix and compiled metadata
#'
#' @export



RunSystemToCell <- function(node.object,
                            ground.truth,
                            assay,
                            blend = 'mean',
                            meta.data.to.map,
                            output_format
                            ){
  

  ### CREATE MAPPING ###
  
  # Ligand data
  subunit.list <- list() # Builds sending (ligand) data for any number of ligand subunits
  for (s in 1:ncol(ground.truth$source.subunits)){ #For each subunit column...

    subunit.list[[s]] <- matrix(data = 1,nrow = nrow(ground.truth$source.subunits),ncol = ncol(node.object)) #initialize a mechanism x barcode matrix of all NAs
    colnames(subunit.list[[s]]) <- colnames(node.object)
    rownames(subunit.list[[s]]) <- rownames(ground.truth$source.subunits)
    non.na.indices <- !is.na(ground.truth$source.subunits[,s]) #Identify rows in the s-th column of the ground truth which are not NA
    subunit.list[[s]][non.na.indices,] <- as.matrix(GetSeuratAssay(node.object,assay,"data")[ground.truth$source.subunits[non.na.indices,s],])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices

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
  if (blend == 'mean.adj'){
    binary <- lig.map > 0
    percent.expression <- Matrix::rowSums(binary)/ncol(binary)
    lig.map2 <- Matrix::rowMeans(lig.map,dims=1)
    lig.map2 <- lig.map2/percent.expression
  }
  lig.map2 <- do.call(cbind, replicate(ncol(lig.map), lig.map2, simplify=FALSE))
  
  # Receptor data
  subunit.list <- list() # Builds receiving (receptor) data for any number of receptor subunits
  for (t in 1:ncol(ground.truth$target.subunits)){

    subunit.list[[t]] <- matrix(data = 1,nrow = nrow(ground.truth$target.subunits),ncol = ncol(node.object)) #initialize a mechanism x barcode matrix of all NAs
    colnames(subunit.list[[t]]) <- colnames(node.object)
    rownames(subunit.list[[t]]) <- rownames(ground.truth$target.subunits)
    non.na.indices <- !is.na(ground.truth$target.subunits[,t]) #Identify rows in the s-th column of the ground truth which are not NA
    subunit.list[[t]][non.na.indices,] <- as.matrix(GetSeuratAssay(node.object,assay,"data")[ground.truth$target.subunits[non.na.indices,t],])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices

  }
  rec.map <- Reduce('*',subunit.list)
  rm(subunit.list)
  
  # Merged map (can be done with any operator, here is multiplication (RECOMMENDED: preserves zeroes and is quantitative))
  sc.connectome <- lig.map2*rec.map


  # Create the rownames (directed ligands and receptors)
  rownames(sc.connectome) <- paste(rownames(lig.map),rownames(rec.map),sep = '—')
  # Create the column names (directed System-cell)
  colnames(sc.connectome) <- paste("System",colnames(rec.map),sep = '—')
  
  #Use this matrix to create a Seurat object:
  demo <- Seurat::CreateSeuratObject(counts = as.matrix(sc.connectome),assay = 'SystemToCell')
  # JC: Seurat V5 will not create data slot automatically, the following step is to manually add this slot
  if(SeuratObject::Version(demo) >= "5.0.0"){
    demo <- Seurat::NormalizeData(demo,assay = "SystemToCell")  # Seurat Object need to be >= 5.0.1
    demo@assays$SystemToCell@layers$data <- demo@assays$SystemToCell@layers$counts # Seurat Object need to be >= 5.0.1
    
  }
  
  # Add metadata to the Seurat object
  meta.data.to.add <- data.frame(as.character(colnames(rec.map)))
  rownames(meta.data.to.add) <- paste("System",colnames(rec.map),sep = '—')
  demo <- Seurat::AddMetaData(demo,metadata = meta.data.to.add,col.name = 'ReceivingCell')
  
  # bug fix: add the system - prefix

  receiving_type.meta <- data.frame(Seurat::Idents(node.object))

  rownames(receiving_type.meta) <- paste("System",rownames(receiving_type.meta),sep = '—')
  
  demo <- Seurat::AddMetaData(demo,metadata = receiving_type.meta,col.name = 'ReceivingType')
  
  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
    # Identify sending and receiving barcodes
    #sending.barcodes <- colnames(do.call(cbind,lig.data)) # Only receiving cell metadata applies for this function
    receiving.barcodes <- colnames(rec.map) 
    # Pull and format sending and receiving metadata
    #sending.metadata <- as.matrix(object@meta.data[,meta.data.to.map][sending.barcodes,])

    # jc: possible bug, change object to node.object
    receiving.metadata <- as.matrix(node.object@meta.data[,meta.data.to.map,drop=FALSE][receiving.barcodes,])

    # Make joint metadata
    #datArray <- abind(sending.metadata,receiving.metadata,along=3)
    #joint.metadata <- as.matrix(apply(datArray,1:2,function(x)paste(x[1],"-",x[2])))
    # Define column names
    #colnames(joint.metadata) <- paste(colnames(sending.metadata),'Joint',sep = '.')
    #colnames(sending.metadata) <- paste(colnames(sending.metadata),'Sending',sep='.')
    #colnames(receiving.metadata) <- paste(colnames(receiving.metadata),'Receiving',sep='.')
    # Compile
    meta.data.to.add.also <- receiving.metadata
    rownames(meta.data.to.add.also) <- paste('System',receiving.barcodes,sep='—')
    # Add additional metadata
    demo <- Seurat::AddMetaData(demo,metadata = as.data.frame(meta.data.to.add.also))
  }
  
  # Define initial identity
  Seurat::Idents(demo) <- demo$ReceivingType
  
  # How many vectors were captured by this sampling?
  
  message(paste("\n",length(unique(demo$ReceivingCell)),'System-To-Cell edges were computed, across',length(unique(demo$ReceivingType)),'cell types'))
  
  if(output_format == "seurat") return(demo)
  else{
    output_list <- vector(mode = "list",length=2)
    names(output_list) <- c("SystemToCellMatrix","metadata")
    output_list[["SystemToCellMatrix"]] <- GetSeuratAssay(demo,"SystemToCell","counts")
    output_list[["metadata"]] <- demo@meta.data
    return(output_list)
  }
  
}

