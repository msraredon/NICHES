LoadOmniPath <- function(){
  
#OmniPath
library(OmnipathR)
library(dplyr)
library(stringr)

#Ligand-Receptor Network ----
lr_Interactions_Omnipath <- import_ligrecextra_interactions()

# Identify max number of ligand subunits and max number of receptor subunits
source_sub_max <- max(str_count(lr_Interactions_Omnipath$source_genesymbol, "_"))
target_sub_max <- max(str_count(lr_Interactions_Omnipath$target_genesymbol, "_"))

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

# Export ligand subunit dataframe
source_subunits <<- temp[,source_col_names]
target_subunits <<- temp[,target_col_names]


}
