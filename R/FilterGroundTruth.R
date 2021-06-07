FilterGroundTruth <- function(ground.truth,input_rownames){
  
  source.indices <- which(apply(ground.truth$source.subunits, 1, function(r) any(r %in% input_rownames)))
  target.indices <- which(apply(ground.truth$target.subunits, 1, function(r) any(r %in% input_rownames)))
  indices.keep <- intersect(source.indices,target.indices)
  
  ground.truth$source.subunits <- ground.truth$source.subunits[indices.keep,]
  ground.truth$target.subunits <- ground.truth$target.subunits[indices.keep,]
  
  return(ground.truth)
}