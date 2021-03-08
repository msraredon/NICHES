#' RunCellToSystem
#' 
#' Condenses all signaling edges coming from each cell within a Seurat object connecting to every other cell in the system, including itself. Outputs another Seurat object, but where the rows of the matrix are ligand-receptor mechanisms
#' and the columns are each a single cell barcode. The information in the matrix is a sum (or an average, depending on user preference) of
#' all signaling edges coming from that particular cell, to every other cell in the system (including to itself.)
#' This transformation allows rapid manipulation and dimensional reduction of how a cell is connected within the system.
#' The default assay of this object is called "CellToSystem" to distinguish it from other Seurat objects.
#' Meta.data slots by default contain "SendingType" information, which is the celltypes for each point, 
#' and "SendingCell" which is the exact cell barcode present in the original Seurat object
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
                      LR.database = 'fantom5',
                      species,
                      assay = 'RNA',
                      min.cells.per.ident = 1,
                      blend = 'sum'){
  
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
  
  # Make SUMMED RECEPTOR INFO
  rec.map <- sys.small@assays[[assay]]@data[receptors,]
 
  if (blend == 'sum'){
    rec.map2 <- Matrix::rowSums(rec.map,dims = 1)
  }
  if (blend == 'mean'){
    rec.map2 <- Matrix::rowMeans(rec.map,dims = 1)
  }
  rec.map2 <- do.call(cbind, replicate(ncol(rec.map), rec.map2, simplify=FALSE))
  
  # Ligand Map from imputed slot
  lig.map <- sys.small@assays[[assay]]@data[ligands,]

  
  # Merged map (can be done with any operator, here is multiplication (RECOMMENDED: preserves zeroes and is quantitative))
  sc.connectome <- lig.map*rec.map2
  
  
  # Create the rownames (directed ligands and receptors)
  rownames(sc.connectome) <- paste(rownames(lig.map),rownames(rec.map),sep = '-')
  # Create the column names (directed cell-System)
  colnames(sc.connectome) <- paste(colnames(lig.map),"System",sep = '-')
  
  # Use this matrix to create a Seurat object:
  demo <- CreateSeuratObject(counts = as.matrix(sc.connectome),assay = 'CellToSystem')
  
  # Add metadata to the Seurat object
  meta.data.to.add <- data.frame(as.character(colnames(lig.map)))
  rownames(meta.data.to.add) <- paste(colnames(lig.map),"System",sep = '-')
  demo <- AddMetaData(demo,metadata = meta.data.to.add,col.name = 'SendingCell')
  demo <- AddMetaData(demo,metadata = Idents(sys.small),col.name = 'SendingType')
  
  # How many vectors were captured by this sampling?
  
  message(paste("\n",length(unique(demo$SendingCell)),'Cell-System edges were computed, across',length(unique(demo$SendingType)),'cell types'))
  return(demo)
}

