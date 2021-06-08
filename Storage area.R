# Temp Holding Area

# Ligand dataset (listwise, for each celltype)
lig.list <- list()
for (i in 1:length(celltypes)){
  temp <- subset(sys.small,idents = celltypes[i])
  subunit.list <- list() # Builds sending (ligand) data for any number of ligand subunits
  for (s in 1:ncol(ground.truth$source.subunits)){ #For each subunit column...
    subunit.list[[s]] <- matrix(data = NA_real_,nrow = nrow(ground.truth$source.subunits),ncol = ncol(temp)) #initialize a mechanism x barcode matrix of all NAs
    colnames(subunit.list[[s]]) <- colnames(temp)
    rownames(subunit.list[[s]]) <- rownames(ground.truth$source.subunits)
    non.na.indices <- !is.na(ground.truth$source.subunits[,s]) #Identify rows in the s-th column of the ground truth which are not NA
    subunit.list[[s]][non.na.indices,] <- as.matrix(temp@assays[[assay]]@data[ground.truth$source.subunits[non.na.indices,s],])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
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
    subunit.list[[t]] <- matrix(data = NA_real_,nrow = nrow(ground.truth$target.subunits),ncol = ncol(temp)) #initialize a mechanism x barcode matrix of all NAs
    colnames(subunit.list[[t]]) <- colnames(temp)
    rownames(subunit.list[[t]]) <- rownames(ground.truth$target.subunits)
    non.na.indices <- !is.na(ground.truth$target.subunits[,t]) #Identify rows in the t-th column of the ground truth which are not NA
    subunit.list[[t]][non.na.indices,] <- as.matrix(temp@assays[[assay]]@data[ground.truth$target.subunits[non.na.indices,t],])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
  }
  rec.list[[i]] <- Reduce('*',subunit.list)
  rm(subunit.list)
}