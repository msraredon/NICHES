#' RunCellToCell
#' 
#' Performs cell-cell transformation on a Seurat object, with structural downsampling to avoid data-inflation. Outputs another Seurat object, but where the columns of the matrix are 
#' barcode-barcode pairs, and the rows of the matrix are ligand-receptor mechanisms. This allows rapid manipulation and dimensional reduction of cell-cell connectivity data.
#' The default assay of this object is called "CellToCell" to distinguish it from normal Seurat objects.
#' Meta.data slots by default contain "SendingType" "ReceivingType" and "VectorType" information.
#' 
#' @param sys.small A filtered Seurat object. The active identity will be used to define populations for connectomic sampling and crossings. 
#' @param ground.truth Ground truth signaling mechanisms present in sys.small.
#' @param assay The assay to run the SCC transformation on. Defaults to "RNA."
#' @param meta.data.to.map A character vector of metadata names present in the original object which will be carried to the SCC objects
#' @param output_format string. Choice of the output format. "seurat" will output a list of seurat objects, "raw" will output a list of lists with raw interaction matrix and compiled metadata
#'
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export


RunCellToCell <- function(sys.small,
                          ground.truth,
                          assay,
                          meta.data.to.map,
                          output_format){

  
  ### CREATE MAPPING ###

  # jc: Identify celltypes:names(table(Idents(sys.small))). Better to run check_celltypes, but harder to check
  celltypes <- return_celltypes(sys.small)
  
  # Ligand dataset (listwise, for each celltype)
  lig.list <- list()
  for (i in 1:length(celltypes)){
    temp <- subset(sys.small,idents = celltypes[i])
    subunit.list <- list() # Builds sending (ligand) data for any number of ligand subunits
    for (s in 1:ncol(ground.truth$source.subunits)){ #For each subunit column...
      subunit.list[[s]] <- matrix(data = 1,nrow = nrow(ground.truth$source.subunits),ncol = ncol(temp)) #initialize a mechanism x barcode matrix of all NAs
      colnames(subunit.list[[s]]) <- colnames(temp)
      rownames(subunit.list[[s]]) <- rownames(ground.truth$source.subunits)
      non.na.indices <- !is.na(ground.truth$source.subunits[,s]) #Identify rows in the s-th column of the ground truth which are not NA
      subunit.list[[s]][non.na.indices,] <- as.matrix(getSeuratAssay(temp,assay,"data")[ground.truth$source.subunits[non.na.indices,s],])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
    }
    lig.list[[i]] <- Reduce('*',subunit.list)
    rm(subunit.list)
  }
  
  # Receptor dataset (listwise, for each celltype)
  rec.list <- list()
  for (i in 1:length(celltypes)){
    temp <- subset(sys.small,idents = celltypes[i])
    subunit.list <- list() # Builds receiving (receptor) data for any number of receptor subunits
    for (t in 1:ncol(ground.truth$target.subunits)){
      subunit.list[[t]] <- matrix(data = 1,nrow = nrow(ground.truth$target.subunits),ncol = ncol(temp)) #initialize a mechanism x barcode matrix of all NAs
      colnames(subunit.list[[t]]) <- colnames(temp)
      rownames(subunit.list[[t]]) <- rownames(ground.truth$target.subunits)
      non.na.indices <- !is.na(ground.truth$target.subunits[,t]) #Identify rows in the t-th column of the ground truth which are not NA
      subunit.list[[t]][non.na.indices,] <- as.matrix(getSeuratAssay(temp,assay,"data")[ground.truth$target.subunits[non.na.indices,t],])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
    }
    rec.list[[i]] <- Reduce('*',subunit.list)
    rm(subunit.list)
  }
  
  # For each sending celltype, create all the outgoing edges 
  # to all receiving celltypes -- this covers autocrine AND bi-directional signaling -- 
  # sampled to be representative of all crosses and manageable in size
  lig.data <- list()
  rec.data <- list()
  scc.data <- list()
  sending.cell.idents <- list()
  receiving.cell.idents <- list()
  
  for (i in 1:length(celltypes)){
    
    # Define maximum number of comparisons for each pairing
    num <- as.data.frame(table(Seurat::Idents(sys.small)))
    num$sender.freq <- ncol(lig.list[[i]])
    rownames(num) <- num$Var1
    num <- num[,-1]
    num <- num %>% dplyr::rowwise() %>% dplyr::mutate(max.possible = min(.data$Freq, .data$sender.freq))
    
    # Craft the ligand side for a single sender to all other types
    lig.temp <- list()
    for (j in 1:length(celltypes)){ # here 'j' is every receiving cell type
      lig.temp[[j]] <- lig.list[[i]][,sample(ncol(lig.list[[i]]), size = num[j,]$max.possible), drop = FALSE]
    }
    lig.data[[i]] <- do.call(cbind,lig.temp)
    
    # Craft the receptor side for a single sender to all other types
    rec.temp <- list()
    for (j in 1:length(celltypes)){
      rec.temp[[j]] <- rec.list[[j]][,sample(ncol(rec.list[[j]]), size = num[j,]$max.possible), drop = FALSE]
    }
    rec.data[[i]] <- do.call(cbind,rec.temp)
    
    # Combine into partial SCC matrix  
    scc.data[[i]] <- lig.data[[i]]*rec.data[[i]]
    
    rownames(scc.data[[i]]) <- paste(rownames(lig.data[[i]]),rownames(rec.data[[i]]),sep = '—')
    colnames(scc.data[[i]]) <- paste(colnames(lig.data[[i]]),colnames(rec.data[[i]]),sep = '—')
    
    sending.cell.idents[[i]] <- as.character(Seurat::Idents(sys.small)[colnames(lig.data[[i]])])
    receiving.cell.idents[[i]] <- as.character(Seurat::Idents(sys.small)[colnames(rec.data[[i]])])
    
  }
  
  # Combine all of these to make the full SCC matrix
  scc <- do.call(cbind,scc.data)
  
  #Use this matrix to create a Seurat object:
  demo <- Seurat::CreateSeuratObject(counts = as.matrix(scc),assay = 'CellToCell')
  
  # Gather and assemble metadata based on "ident" slot
  sending.cell.idents.2 <- do.call(c,sending.cell.idents)
  receiving.cell.idents.2 <- do.call(c,receiving.cell.idents)
  meta.data.to.add <- data.frame(SendingType = sending.cell.idents.2,
                                 ReceivingType = receiving.cell.idents.2)
  rownames(meta.data.to.add) <- colnames(scc)
  meta.data.to.add$VectorType <- paste(meta.data.to.add$SendingType,
                                       meta.data.to.add$ReceivingType,
                                       sep = '—')
  #Add ident metadata
  demo <- Seurat::AddMetaData(demo,metadata = meta.data.to.add)
  
  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
  # Identify sending and receiving barcodes
  sending.barcodes <- colnames(do.call(cbind,lig.data)) #This can be simplified if the above SCC construction is simplified
  receiving.barcodes <- colnames(do.call(cbind,rec.data)) #This can be simplified if the above SCC construction is simplified
  # Pull and format sending and receiving metadata    jc: possible bug, change object to sys.small
  sending.metadata <- as.matrix(sys.small@meta.data[,meta.data.to.map,drop=FALSE][sending.barcodes,])
  receiving.metadata <- as.matrix(sys.small@meta.data[,meta.data.to.map,drop=FALSE][receiving.barcodes,])
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
  # Define initial identity
  Seurat::Idents(demo) <- demo$VectorType
  
  # How many vectors were captured by this sampling?
  message(paste("\n",ncol(demo),'Cell-To-Cell edges computed, sampling',length(unique(demo$VectorType)),'distinct VectorTypes, out of',length(table(Seurat::Idents(sys.small)))^2,'total possible'))
  
  if(output_format == "seurat") return(demo)
  else{
    output_list <- vector(mode = "list",length=2)
    names(output_list) <- c("CellToCellMatrix","metadata")
    output_list[["CellToCellMatrix"]] <- getSeuratAssay(demo,"CellToCell","counts")
    output_list[["metadata"]] <- demo@meta.data
    return(output_list)
  }

}

