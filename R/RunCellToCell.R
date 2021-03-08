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
#' @param min.cells.per.ident Default 10. A limit on how small (how many cells) a single population can be to participate in connectomic crossings.
#'
#' @export


RunCellToCell <- function(object,
                   LR.database = 'fantom5',
                   species,
                   assay = 'RNA',
                   min.cells.per.ident = 10){

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

  # Identify celltypes
  celltypes <- names(table(Idents(sys.small)))
  
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
  
  # Cool, but in order to interpret, we need the additional metadata of cell types so that we can color it 
  # by sending cell type, receiving cell type, and overall celltype-to-celltype vector. 
  sending.cell.idents.2 <- do.call(c,sending.cell.idents)
  receiving.cell.idents.2 <- do.call(c,receiving.cell.idents)
  meta.data.to.add <- data.frame(SendingType = sending.cell.idents.2,
                                 ReceivingType = receiving.cell.idents.2)
  rownames(meta.data.to.add) <- colnames(scc)
  meta.data.to.add$VectorType <- paste(meta.data.to.add$SendingType,
                                       meta.data.to.add$ReceivingType,
                                       sep = '-')
  
  #Add metadata to the Seurat object
  demo <- AddMetaData(demo,metadata = meta.data.to.add)
  
  # How many vectors were captured by this sampling?
 
  message(paste("\n",length(unique(demo$VectorType)),'distinct VectorTypes were computed, out of',length(table(Idents(sys.small)))^2,'total possible'))
  return(demo)
}

