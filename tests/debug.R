library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(NICHES)

library(SeuratData)
#InstallData("ifnb")
data("ifnb")

data_list <- SplitObject(ifnb, split.by="stim")

niches_obs <- lapply(data_list, function(x){
  RunNICHES(x,
            assay = 'RNA',
            species = 'human',
            LR.database = 'fantom5',
            SystemToCell = T,
            CellToCell = F) 
}) 

niche_merge <- merge(niches_obs[[1]]$SystemToCell,niches_obs[[2]]$SystemToCell,add.cell.ids = c("CTRL","STIM"))
