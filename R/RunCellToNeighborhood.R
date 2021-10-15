#' RunCellToNeighborhood
#'
#' @param object A Seurat 4.0 object. The active identity will be used to define populations for connectomic sampling and crossings.
#' @param LR.database Currently accepts 'fantom5' or 'omnipath' 
#' @param species The species of the object that is being processed.  Only required if LR.database = 'fantom5', and allows 'human','mouse','rat', or 'pig'
#' @param assay The assay to run the SCC transformation on. Defaults to "RNA."
#' @param min.cells.per.ident A limit on how small (how many cells) a single population can be to participate in connectomic crossings.
#' @param min.cells.per.gene Limits analysis to interactions involving genes expressed above minimum threshold number of cells in the system. 
#' @param position.x The name of the meta.data column specifying location on the spatial x-axis. Only relevant for spatial omics data.
#' @param position.y The name of the meta.data column specifying location on the spatial y-axis. Only relevant for spatial omics data.
#' @param rad.set The radius in Euclidean space to consider local neighbors.
#' @param meta.data.to.map A character vector of metadata names present in the original object which will be carried to the NICHES objects
#' @param ... 
#'
#' @export

RunCellToNeighborhood <- function(object,
                               LR.database,
                               species,
                               assay,
                               min.cells.per.ident,
                               min.cells.per.gene,
                               position.x,
                               position.y,
                               rad.set,
                               meta.data.to.map,...){

  # jc: wrapped the preprocessing steps
  sys.small <- prepSeurat(object,assay,min.cells.per.ident,min.cells.per.gene)

  # jc: Load corresponding ligands and receptors
  ground.truth <- lr_load(LR.database,species,rownames(sys.small@assays[[assay]]))

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
  rad = rad.set
  result <- apply(df, 1, function(pt)
    (sqrt(abs(pt["x"] - df$x)^2 + abs(pt["y"] - df$y)^2) <= rad)
  )

  diag(result) <- 1
  rownames(result) <- colnames(result)

  # Convert adj matrix to edgelist
  edgelist <- igraph::graph.adjacency(result)
  edgelist <- igraph::get.data.frame(edgelist)

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
  rownames(scc) <- paste(rownames(lig.data),rownames(rec.data),sep = '-')

  # Condense by column name
  colnames(scc) <- colnames(lig.data) # Make colnames equal to sending cell
  scc <- as.matrix(scc)
  scc <- t(rowsum(t(scc), colnames(scc)))

  # Label columns properly
  barcodes <- colnames(scc)
  colnames(scc) <- paste(colnames(scc),'Neighborhood',sep = '-')

  # Use this matrix to create a Seurat object:
  demo <- Seurat::CreateSeuratObject(counts = as.matrix(scc),assay = 'CellToNeighborhood')

  # Add metadata based on ident slot
  demo <- Seurat::AddMetaData(demo,metadata = barcodes,col.name = 'SendingCell')
  demo <- Seurat::AddMetaData(demo,metadata = Seurat::Idents(sys.small)[barcodes],col.name = 'SendingType')

  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
    # Identify sending and receiving barcodes
    sending.barcodes <- barcodes # Only sending cell metadata applies for this function
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
    rownames(meta.data.to.add.also) <- paste(sending.barcodes,'Neighborhood',sep='-')
    # Add additional metadata
    demo <- Seurat::AddMetaData(demo,metadata = as.data.frame(meta.data.to.add.also))
  }
  # Set initial identity
  Seurat::Idents(demo) <- demo$SendingType
  # How many vectors were captured by this sampling?
  message(paste("\n",length(unique(demo$SendingCell)),'Cell-To-Neighborhood edges were computed, across',length(unique(demo$SendingType)),'cell types'))

  return(demo)
}
