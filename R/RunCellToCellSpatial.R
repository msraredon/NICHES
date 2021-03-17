#' RunCellToCellSpatial
#' 
#' @param object A Seurat 4.0 object. The active identity will be used to define populations for connectomic sampling and crossings.
#' @param LR.database Accepts either 'fantom5' or a custom data.frame with the first column equal to ligands, second column equal to associated receptors.
#' @param species The species of the object that is being processed.  Only required if LR.database = 'fantom5', and allows 'human','mouse','rat', or 'pig'
#' @param assay The assay to run the SCC transformation on. Defaults to "RNA."
#' @param min.cells.per.ident Default 1. A limit on how small (how many cells) a single population can be to participate in connectomic crossings.
#' @param position.x The name of the meta.data column specifying location on the spatial x-axis. Only relevant for spatial omics data.
#' @param position.y The name of the meta.data column specifying location on the spatial y-axis. Only relevant for spatial omics data.
#'
#' @export

RunCellToCellSpatial <- function(object,
                          LR.database = 'fantom5',
                          species,
                          assay = 'RNA',
                          min.cells.per.ident = 1,
                          position.x,
                          position.y,
                          meta.data.to.map = NULL){
  
  require(Seurat)
  require(dplyr)
  
  # Check if setup is correct
  if (class(LR.database) == 'character'){
    if (LR.database == 'fantom5' & is.null(species)){stop("\nPlease select species for FANTOM5 mapping. Allows 'human','mouse','rat', or 'pig' ")}
  }else{}
  
  # Set default assay (not necessary, but just in case)
  DefaultAssay(object) <- assay
  
  # Stash object
  sys.small <- object
  
  # Limit object to cell populations larger than requested minimum
  if (!is.null(min.cells.per.ident)){
    message(paste("\n",'Subsetting to populations with greater than',min.cells.per.ident,'cells'))
    idents.include <- names(table(Idents(sys.small)))[table(Idents(sys.small)) > min.cells.per.ident]
    sys.small <- subset(sys.small,idents = idents.include)
  }
  
  num.cells <- ncol(sys.small)
  message(paste("\n",num.cells,'distinct cells from',length(names(table(Idents(sys.small)))),'celltypes to be analyzed'))
  
  # Identify paired ligands and receptors in the dataset to map against
  if(class(LR.database) == 'character'){
    if (LR.database == 'fantom5'){
      # Load ground-truth database (FANTOM5, species-converted as appropriate, per methodlogy in Raredon et al 2019, DOI: 10.1126/sciadv.aaw3851)
      if (species == 'human'){
        fantom <- Connectome::ncomms8866_human
      }
      if (species == 'mouse'){
        fantom <- Connectome::ncomms8866_mouse
      }
      if (species == 'rat'){
        fantom <- Connectome::ncomms8866_rat
      }
      if (species == 'pig'){
        fantom <- Connectome::ncomms8866_pig
      }}}else{
        num.mechs <- nrow(LR.database)
        message(paste("\n","Custom mapping requested. Mapping cells against",num.mechs,"mechanisms provided via LR.database argument"))
        fantom <- data.frame(Ligand.ApprovedSymbol = as.character(LR.database[,1]),
                             Receptor.ApprovedSymbol = as.character(LR.database[,2]))
      }
  
  # Subset to only mechanisms present in the object
  fantom.specific <- subset(fantom,
                            Ligand.ApprovedSymbol %in% rownames(sys.small@assays[[assay]]) & Receptor.ApprovedSymbol %in% rownames(sys.small@assays[[assay]]))
  ligands <- fantom.specific$Ligand.ApprovedSymbol
  receptors <- fantom.specific$Receptor.ApprovedSymbol
  
  
  ### CREATE MAPPING ###
  
  # Create adjacency matrix
  # Adapted from :: https://stackoverflow.com/questions/16075232/how-to-create-adjacency-matrix-from-grid-coordinates-in-r
  # Setup numbering and labeling
  df <- data.frame(x = object[[position.x]], y = object[[position.y]])
  df$barcode <- rownames(df)
  df$x <- as.character(df$x)
  df$y <- as.character(df$y)
  df$x <- as.numeric(df$x)
  df$y <- as.numeric(df$y)
  df <- df[,c('x','y')]
  
  # Make adj matrix
  # Within a circle of radius "rad" around each coordinate (Set rad = 1 for only direct neighbors)
  rad = 1
  result <- apply(df, 1, function(pt) 
    (sqrt(abs(pt["x"] - df$x)^2 + abs(pt["y"] - df$y)^2) <= rad) 
  )
  
  diag(result) <- 1
  rownames(result) <- colnames(result)
  
  # Convert adj matrix to edgelist
  edgelist <- igraph::graph.adjacency(result)
  edgelist <- igraph::get.data.frame(edgelist)
  
  # Make ligand matrix

    lig.data <- sys.small@assays[[assay]]@data[ligands,edgelist$from]
  
  # Make receptor matrix

    rec.data <- sys.small@assays[[assay]]@data[receptors,edgelist$to]
  
  
  # Make SCC matrix
  scc <- lig.data*rec.data
  rownames(scc) <- paste(rownames(lig.data),rownames(rec.data),sep = '-')
  colnames(scc) <- paste(colnames(lig.data),colnames(rec.data),sep = '-')
  sending.cell.idents <- as.character(Idents(sys.small)[colnames(lig.data)])
  receiving.cell.idents <- as.character(Idents(sys.small)[colnames(rec.data)])
  
  # Use this matrix to create a Seurat object:
  demo <- CreateSeuratObject(counts = as.matrix(scc),assay = 'CellToCellSpatial')
  
  # Add key metadata
  
  meta.data.to.add <- data.frame(SendingType = sending.cell.idents,
                                 ReceivingType = receiving.cell.idents)
  
  rownames(meta.data.to.add) <- colnames(scc)
  meta.data.to.add$VectorType <- paste(meta.data.to.add$SendingType,
                                       meta.data.to.add$ReceivingType,
                                       sep = '-')
  
  #Add metadata to the Seurat object
  demo <- AddMetaData(demo,metadata = meta.data.to.add)
  
  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
    # Identify sending and receiving barcodes
    sending.barcodes <- colnames(lig.data) 
    receiving.barcodes <- colnames(rec.data) 
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
  message(paste("\n",length(unique(demo$VectorType)),'distinct VectorTypes were computed, out of',length(table(Idents(sys.small)))^2,'total possible'))
  
  return(demo)
}

