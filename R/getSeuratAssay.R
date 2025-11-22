#' Wrapper function for GetAssayData based on the version of the SeuratObject
#'
#' @param object input Seurat object 
#' @param assay which assay to use in object
#' @param slot "data","counts",or "scale.data"
#'
#' @return the specified assay data
#' @export
# GetSeuratAssay <- function(object,assay,slot){
#   if(SeuratObject::Version(object) >= 5) return(Seurat::GetAssayData(object,assay=assay,layer=slot))
#   else return(Seurat::GetAssayData(object,assay=assay,slot=slot))
# }
# per https://github.com/msraredon/NICHES/issues/43#issuecomment-2377502362
getSeuratAssay <- function (object, assay, slot) 
{
  if (SeuratObject::Version(object) >= "5.0.0") 
    Seurat::GetAssayData(object, assay = assay, layer = slot)
  else Seurat::GetAssayData(object, assay = assay, slot = slot)
}