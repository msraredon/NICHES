
prepSeurat <- function(object,assay,min.cells.per.ident){
  # Set default assay (not necessary, but just in case)
  Seurat::DefaultAssay(object) <- assay
  
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
  
  return(sys.small)
}


# load the L-R that corresponds to the row names of the input seurat
lr_load <- function(LR.database,species,input_rownames){
  
  if(class(LR.database) == 'character'){
    if (LR.database == 'fantom5'){
      if(is.null(species)){
        stop("\nPlease select species for FANTOM5 mapping. Allows 'human','mouse','rat', or 'pig' ")}
      # Load ground-truth database (FANTOM5, species-converted as appropriate, per methodlogy in Raredon et al 2019, DOI: 10.1126/sciadv.aaw3851)
      else{
        switch(species,
               human = {fantom <- get(data(ncomms8866_human))},
               mouse = {fantom <- get(data(ncomms8866_mouse))},
               rat = {fantom <- get(data(ncomms8866_rat))},
               pig = {fantom <- get(data(ncomms8866_pig))},
               stop("input species not recognized"))
      } 
    }
  }
  else{
    num.mechs <- nrow(LR.database)
    message(paste("\n","Custom mapping requested. Mapping cells against",
                  num.mechs,"mechanisms provided via LR.database argument"))
    fantom <- data.frame(Ligand.ApprovedSymbol = as.character(LR.database[,1]),
                         Receptor.ApprovedSymbol = as.character(LR.database[,2]))
  }
  
  
  # Subset to only mechanisms present in the object
  fantom.specific <- subset(fantom,
                            fantom$Ligand.ApprovedSymbol %in% input_rownames 
                            & fantom$Receptor.ApprovedSymbol %in% input_rownames)
  ligands <- fantom.specific$Ligand.ApprovedSymbol
  receptors <- fantom.specific$Receptor.ApprovedSymbol
  
  return(list("ligands" = ligands, "receptors" = receptors))
}

# better: check_celltypes: check whether the idents are cell types, yes to return the unique cell types, no to return an error
return_celltypes <- function(seurat_object){
  warning("Please make sure that Identity of the input seurat object corresponds to cell types")
  return(unique(Idents(sys.small)))
}
