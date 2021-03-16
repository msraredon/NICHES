#' RunCellToCell
#' 
#' Performs cell-cell transformation on a Seurat object, with structural downsampling to avoid data-inflation. Outputs another Seurat object, but where the columns of the matrix are 
#' barcode-barcode pairs, and the rows of the matrix are ligand-receptor mechanisms. This allows rapid manipulation and dimensional reduction of cell-cell connectivity data.
#' The default assay of this object is called "CellToCell" to distinguish it from normal Seurat objects.
#' Meta.data slots by default contain "SendingType" "ReceivingType" and "VectorType" information.
#' 
#' @param object A Seurat 3.0 object.  The active identity meta.data will be used to define populations for connectomic sampling and crossings.
#' @param LR.database Accepts either 'fantom5' or a custom data.frame with the first column equal to ligands, second column equal to associated receptors.
#' @param species The species of the object that is being processed.  Only required if LR.database = 'fantom5', and allows 'human','mouse','rat', or 'pig'
#' @param assay The assay to run the SCC transformation on. Defaults to "RNA."
#' @param min.cells.per.ident Default 1. A limit on how small (how many cells) a single population can be to participate in connectomic crossings.
#' @param meta.data.to.map A character vector of metadata names present in the original object which will be carried to the SCC objects
#'
#' @export


RunCellToCell <- function(object,
                   LR.database = 'fantom5',
                   species,
                   assay = 'RNA',
                   min.cells.per.ident = 1,
                   meta.data.to.map = NULL){

  require(Seurat)
  require(dplyr)
  require(abind)
  
  # jc: wrapped the preprocessing steps
  sys.small <- prepSeurat(object,assay,min.cells.per.ident)
  
  # jc: Load corresponding ligands and receptors
  lrs <- lr_load(LR.database,species,rownames(sys.small@assays[[assay]]))
  ligands <- lrs[['ligands']]
  receptors <- lrs[['receptors']]
  
  ### CREATE MAPPING ###

  # jc: Identify celltypes:names(table(Idents(sys.small))). Better to run check_celltypes, but harder to check
  celltypes <- return_celltypes(sys.small)
  
  # Ligand dataset
  lig.list <- list()
  for (i in 1:length(celltypes)){
    temp <- subset(sys.small,idents = celltypes[i])
    lig.list[[i]] <- temp@assays[[assay]]@data[ligands,]
  }
  
  # Receptor dataset
  rec.list <- list()
  for (i in 1:length(celltypes)){
    temp <- subset(sys.small,idents = celltypes[i])
    rec.list[[i]] <- temp@assays[[assay]]@data[receptors,]
  }
  
  # For each celltype, create all the outgoing edges 
  # (to all celltypes -- this covers autocrine AND bi-directional signaling)
  lig.data <- list()
  rec.data <- list()
  scc.data <- list()
  sending.cell.idents <- list()
  receiving.cell.idents <- list()
  
  for (i in 1:length(celltypes)){
    
    # Define maximum number of comparisons for each pairing
    num <- as.data.frame(table(Idents(sys.small)))
    num$sender.freq <- ncol(lig.list[[i]])
    rownames(num) <- num$Var1
    num <- num[,-1]
    num <- num %>% rowwise() %>% mutate(max.possible = min(Freq, sender.freq))
    
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
    
    rownames(scc.data[[i]]) <- paste(rownames(lig.data[[i]]),rownames(rec.data[[i]]),sep = '-')
    colnames(scc.data[[i]]) <- paste(colnames(lig.data[[i]]),colnames(rec.data[[i]]),sep = '-')
    
    sending.cell.idents[[i]] <- as.character(Idents(sys.small)[colnames(lig.data[[i]])])
    receiving.cell.idents[[i]] <- as.character(Idents(sys.small)[colnames(rec.data[[i]])])
    
  }
  
  # Combine all of these to make the full SCC matrix
  scc <- do.call(cbind,scc.data)
  
  #Use this matrix to create a Seurat object:
  demo <- CreateSeuratObject(counts = as.matrix(scc),assay = 'CellToCell')
  
  # Gather and assemble metadata based on "ident" slot
  sending.cell.idents.2 <- do.call(c,sending.cell.idents)
  receiving.cell.idents.2 <- do.call(c,receiving.cell.idents)
  meta.data.to.add <- data.frame(SendingType = sending.cell.idents.2,
                                 ReceivingType = receiving.cell.idents.2)
  rownames(meta.data.to.add) <- colnames(scc)
  meta.data.to.add$VectorType <- paste(meta.data.to.add$SendingType,
                                       meta.data.to.add$ReceivingType,
                                       sep = '-')
  #Add ident metadata
  demo <- AddMetaData(demo,metadata = meta.data.to.add)
  
  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
  # Identify sending and receiving barcodes
  sending.barcodes <- colnames(do.call(cbind,lig.data)) #This can be simplified if the above SCC construction is simplified
  receiving.barcodes <- colnames(do.call(cbind,rec.data)) #This can be simplified if the above SCC construction is simplified
  # Pull and format sending and receiving metadata
  sending.metadata <- as.matrix(object@meta.data[,meta.data.to.map][sending.barcodes,])
  receiving.metadata <- as.matrix(object@meta.data[,meta.data.to.map][receiving.barcodes,])
  # Make joint metadata
  datArray <- abind(sending.metadata,receiving.metadata,along=3)
  joint.metadata <- as.matrix(apply(datArray,1:2,function(x)paste(x[1],"-",x[2])))
  # Define column names
  colnames(joint.metadata) <- paste(colnames(sending.metadata),'Joint',sep = '.')
  colnames(sending.metadata) <- paste(colnames(sending.metadata),'Sending',sep='.')
  colnames(receiving.metadata) <- paste(colnames(receiving.metadata),'Receiving',sep='.')
  # Compile
  meta.data.to.add.also <- cbind(sending.metadata,receiving.metadata,joint.metadata)
  rownames(meta.data.to.add.also) <- paste(sending.barcodes,receiving.barcodes,sep='-')
  # Add additional metadata
  demo <- AddMetaData(demo,metadata = as.data.frame(meta.data.to.add.also))
  }
  
  # How many vectors were captured by this sampling?
  message(paste("\n",ncol(demo),'Cell-To-Cell edges computed, sampling',length(unique(demo$VectorType)),'distinct VectorTypes, out of',length(table(Idents(sys.small)))^2,'total possible'))
  return(demo)
}

