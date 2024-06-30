
# better: check_celltypes: check whether the idents are cell types, yes to return the unique cell types, no to return an error
#' return cell types
#'
#' @param seurat_object Input seurat object
#'
#' @return unique cell types
#' @export

return_celltypes <- function(seurat_object){
  message('\n For sampling purposes, please make sure that the active Identity of the input seurat object corresponds to cell types')
  return(names(table(Seurat::Idents(seurat_object)))) # ms: other steps are depdent on the exact format of this output
}