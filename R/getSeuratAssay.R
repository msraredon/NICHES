#' Wrapper function for GetAssayData based on the version of the SeuratObject
#'
#' @param object input Seurat object 
#' @param assay which assay to use in object
#' @param slot "data","counts",or "scale.data"
#'
#' @return the specified assay data
#' @export
getSeuratAssay <- function(object,assay,slot){
  if(SeuratObject::Version(object) >= 5) return(Seurat::GetAssayData(object,assay=assay,layer=slot))
  else return(Seurat::GetAssayData(object,assay=assay,slot=slot))
}