require(SCC)
require(Seurat)
# Load Data
emb <- readRDS("~/Box Sync/Kluger_Lab/Spatial SCC/spatial_seurat_processed_main.rds")

#Pull out one slide only
FFPE2 <- subset(emb,idents=c("E10.5 tail"))
DefaultAssay(FFPE2) <- "RNA"
Idents(FFPE2) <- FFPE2$cell_label_dom

# Register position
FFPE2@meta.data$position <- gsub("_2","",colnames(FFPE2))

# Index position to x and y
temp <- data.frame(stringr::str_split_fixed(FFPE2$position,pattern = 'x',n = 2))
FFPE2$x <- temp$X1
FFPE2$y <- temp$X2

# Testing runs
test <- RunSCC(FFPE2,species = 'mouse',position.x = 'x',position.y = 'y',
               CellToCell = T,CellToSystem = T,SystemToCell = T,
               CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F,meta.data.to.map = c('nCount_SCT','orig.ident')) #works

test <- RunSCC(FFPE2,species = 'mouse',position.x = 'x',position.y = 'y',
               CellToCell = F,CellToSystem = T,SystemToCell = F,
               CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F) #works
test <- RunSCC(FFPE2,species = 'mouse',position.x = 'x',position.y = 'y',
               CellCell = F,CellSystem = F,SystemCell = T,
               CellCellSpatial = F,CellNeighborhood = F,NeighborhoodCell = F) #works

test <- RunSCC(FFPE2,species = 'mouse',position.x = 'x',position.y = 'y',
               CellToCell = T,CellToSystem = T,SystemToCell = T,
               CellToCellSpatial = T,CellToNeighborhood = F,NeighborhoodToCell = F) #works

test <- RunSCC(FFPE2,species = 'mouse',position.x = 'x',position.y = 'y',
               CellCell = F,CellSystem = F,SystemCell = F,
               CellCellSpatial = F,CellNeighborhood = T,NeighborhoodCell = F)
test <- RunSCC(FFPE2,species = 'mouse',position.x = 'x',position.y = 'y',
               CellCell = F,CellSystem = F,SystemCell = F,
               CellCellSpatial = F,CellNeighborhood = F,NeighborhoodCell = T)
test <- RunSCC(FFPE2,species = 'mouse',position.x = 'x',position.y = 'y',
               CellCell = T,CellSystem = T,SystemCell = T,
               CellCellSpatial = T,CellNeighborhood = F,NeighborhoodCell = F) #works
test
