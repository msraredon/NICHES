
#OmniPath
library(OmnipathR)
library(dplyr)
library(stringr)

#Ligand-Receptor Network ----
lr_Interactions_Omnipath <- import_ligrecextra_interactions() # %>%
  #dplyr::select(source_genesymbol,target_genesymbol,sources) %>%
  #dplyr::rename(from=source_genesymbol, to=target_genesymbol) %>% 
  #dplyr::filter(from != to) %>% 
  #dplyr::distinct()

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
                        col = source_genesymbol, # Split Target genes
                        into = target_col_names, # Uses initialized column names
                        sep = '_',
                        remove = F)

