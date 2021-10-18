#' Load OmniPath as ground truth
#'
#' @param species Species reference in the Ominipath database
#'
#' @export

LoadOmniPath <- function(species){
  
# Setup species call
  if (species == 'human'){
    organism = 9606
    }else if (species == 'mouse'){
      organism = 10090
      }else if(species == 'rat'){
        organism = 10116
        }else{
          stop("\nPlease select species for OmniPath mapping. Allows 'human','mouse',or 'rat' ")
        }

#Ligand-Receptor Network ----
lr_Interactions_Omnipath <- OmnipathR::import_ligrecextra_interactions(organism = organism) %>%
  dplyr::select(source_genesymbol,target_genesymbol) %>%
  dplyr::distinct()

#Tag with mechanism name
lr_Interactions_Omnipath$mechanism <- paste(lr_Interactions_Omnipath$source_genesymbol,lr_Interactions_Omnipath$target_genesymbol,sep = '-')

# Identify max number of ligand subunits and max number of receptor subunits (based on "_" as a separator, used in current OmniPath iteration as of 2021-06-07)
source_sub_max <- max(stringr::str_count(lr_Interactions_Omnipath$source_genesymbol, "_"))
target_sub_max <- max(stringr::str_count(lr_Interactions_Omnipath$target_genesymbol, "_"))

# Initialize column names based on how many subunits are in initial database
source_col_names <- c()
for (i in 1:length(c(1:source_sub_max))){
  source_col_names[i] <- paste('source_',c(1:source_sub_max)[i],sep='')
}
target_col_names <- c()
for (i in 1:length(c(1:target_sub_max))){
  target_col_names[i] <- paste('target_',c(1:target_sub_max)[i],sep='')
}

# Split into individual columns
temp <- tidyr::separate(data = lr_Interactions_Omnipath,
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
