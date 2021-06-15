#' Filter the ground truth to only include mechanisms in the base object
#'
#' @export

FilterGroundTruth <- function(ground.truth,input_rownames){

  GOI <- c(input_rownames,NA,'NA')
  
  message('\n Limiting ground truth to genes within dataset')
  
  source.indices <- which(apply(ground.truth$source.subunits, 1, function(r) all(r %in% GOI))) #All sending subunits must be in dataset or be NA
  target.indices <- which(apply(ground.truth$target.subunits, 1, function(r) all(r %in% GOI))) #All receiving subunits must be in dataset or be NA
  indices.keep <- intersect(source.indices,target.indices)
  
  ground.truth$source.subunits <- ground.truth$source.subunits[indices.keep,,drop = F]
  ground.truth$target.subunits <- ground.truth$target.subunits[indices.keep,,drop = F]
  
  # Only keep columns which have non-NA values remaining
  ground.truth$source.subunits <- ground.truth$source.subunits[,colSums(!is.na(ground.truth$source.subunits))>0,drop=F]
  ground.truth$target.subunits <- ground.truth$target.subunits[,colSums(!is.na(ground.truth$target.subunits))>0,drop=F]
  
  message(paste('\n Mapping against',nrow(ground.truth$source.subunits),'ground truth signaling mechanisms'))
  
  return(ground.truth)
}
