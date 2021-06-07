#' Filter the ground truth to only include mechanisms in the base object
#'
#' @export

FilterGroundTruth <- function(ground.truth,input_rownames){
  
  source.indices <- which(apply(ground.truth$source.subunits, 1, function(r) any(r %in% input_rownames)))
  target.indices <- which(apply(ground.truth$target.subunits, 1, function(r) any(r %in% input_rownames)))
  indices.keep <- intersect(source.indices,target.indices)
  
  ground.truth$source.subunits <- ground.truth$source.subunits[indices.keep,,drop = F]
  ground.truth$target.subunits <- ground.truth$target.subunits[indices.keep,,drop = F]
  
  # Only keep columns which have non-NA values remaining
  ground.truth$source.subunits <- ground.truth$source.subunits[,colSums(!is.na(ground.truth$source.subunits))>0,drop=F]
  ground.truth$target.subunits <- ground.truth$target.subunits[,colSums(!is.na(ground.truth$target.subunits))>0,drop=F]

  return(ground.truth)
}
