#' Load the custom Ligand Receptor database from the user
#'
#' @param custom_LR_database data.frame. Each row is a ligand-receptor mechanism where the first column corresponds to the source genes that express the ligands subunits (separated by '_') and the second column corresponds to the receptor genes that express the receptor subunits (separated by '_').
#'
#' @export
#'
LoadCustom <- function(custom_LR_database){
  
  if(is.null(custom_LR_database)) stop("Custom Ligand Receptor database is not provided.")
  #TODO: check the validity of the custom database
  message("Using custom Ligand Receptor Database...")
  message("'species' parameter ignored")
  #TODO: check the format of custom_LR_database: 2 columns, each entry is '_' separated
  #TODO: check gene symbols
  
  # custom_LR_database is required to be formatted like OminiPath database 
  # m * 2 dataframe
  # each row is a mechanism, the first column contains ligands with subunits separated by '_', the second column contains receptors with subunits separated by '_'
  colnames(custom_LR_database) <- c("source_genesymbol","target_genesymbol")
  # combine each LR to a mechanism 
  custom_LR_database$mechanism <- paste(custom_LR_database$source_genesymbol,custom_LR_database$target_genesymbol,sep = '-')
  # count the maximum number of subunits for Ligands and Receptors
  source_sub_max <- max(sapply(custom_LR_database$source_genesymbol,function(x) length(strsplit(x,split="_")[[1]])))
  target_sub_max <- max(sapply(custom_LR_database$target_genesymbol,function(x) length(strsplit(x,split="_")[[1]])))
  # initialize the col names
  source_col_names <- paste0("source_",c(1:source_sub_max))
  target_col_names <- paste0("target_",c(1:target_sub_max))
  
  # Split into individual columns
  temp <- tidyr::separate(data = custom_LR_database,
                          col = source_genesymbol, # Split Source genes
                          into = source_col_names, # Uses initialized column names
                          sep = '_',
                          remove = F)
  temp <- tidyr::separate(data = temp,
                          col = target_genesymbol, # Split Target genes
                          into = target_col_names, # Uses initialized column names
                          sep = '_',
                          remove = F)
  # Export subunit dataframe
  source.subunits <- as.matrix(temp[,source_col_names]) #allows duplicate rownames
  rownames(source.subunits) <- temp$source_genesymbol
  target.subunits <- as.matrix(temp[,target_col_names]) #allows duplicate rownames
  rownames(target.subunits) <- temp$target_genesymbol
  
  ground.truth <- list('source.subunits' = source.subunits,
                       'target.subunits' = target.subunits)
  return(ground.truth)
  
  
}