

#' Seurat input preprocessing
#'
#' @param object input Seurat object 
#' @param assay which assay to use in object
#' @param min.cells.per.ident minimum cells per ident
#' @param min.cells.per.gene minimum cells which must express a gene for it to be considered in connectomics
#'
#' @return a processed seurat object
#' @export

prepSeurat <- function(object,assay,min.cells.per.ident,min.cells.per.gene){
  # Set default assay (not necessary, but just in case)
  Seurat::DefaultAssay(object) <- assay
  
  # Stash object
  sys.small <- object
  
  
  # Limit object to cell populations larger than requested minimum
  if (!is.null(min.cells.per.ident)){
    message(paste("\n",'Subsetting to populations with greater than',min.cells.per.ident,'cells'))
    idents.include <- names(table(Seurat::Idents(sys.small)))[table(Seurat::Idents(sys.small)) > min.cells.per.ident]
    sys.small <- subset(sys.small,idents = idents.include)
  }
  
  # Limit analysis to interaction involving genes expressed above minimum threshold number of cells in the system
  if (!is.null(min.cells.per.gene)){
    message(paste("\n",'Subsetting to genes expressed in greater than',min.cells.per.gene,'cells'))
    
    # jc: return an error if the count matrix in the given assay is empty
    if(!length(getSeuratAssay(sys.small,assay,"counts"))) stop("Unable to subset: the count matrix in the given assay is empty.")
    
    cells.per.gene <- data.frame(non.zero.cells = Matrix::rowSums(getSeuratAssay(sys.small,assay,"counts")>0))
    GOI <- subset(cells.per.gene,non.zero.cells > min.cells.per.gene)
    sys.small <- sys.small[rownames(GOI),]
  }
  
  num.cells <- ncol(sys.small)
  message(paste("\n",num.cells,'distinct cells from',
                length(names(table(Seurat::Idents(sys.small)))),'celltypes to be analyzed'))
  
  return(sys.small)
}