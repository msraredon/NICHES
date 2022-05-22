#devtools::document()
#devtools::check()
#devtools::install()
library(NICHES)
library(Seurat)
# Load Data
emb <- readRDS("/data/jyc/github_proj/SCC_paper/code_test/data/spatial_seurat_processed_main.RDS")

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
test <- RunNICHES(object = FFPE2,assay = "SCT",LR.database = 'fantom5',species = 'mouse',position.x = 'x',position.y = 'y',rad.set = 1,k=NULL,
               CellToCell = T,CellToSystem = T,SystemToCell = T,
               CellToCellSpatial = T,CellToNeighborhood = T,NeighborhoodToCell = T,meta.data.to.map = c('nCount_SCT','orig.ident')) #works

test1 <- RunNICHES(object = FFPE2,assay = "SCT",LR.database = 'fantom5',species = 'mouse',position.x = 'x',position.y = 'y',rad.set = 1,k=NULL,
                  CellToCell = T,CellToSystem = T,SystemToCell = T,blend = "sum",
                  CellToCellSpatial = T,CellToNeighborhood = T,NeighborhoodToCell = T,meta.data.to.map = c('nCount_SCT','orig.ident')) #works

# Cluster
Idents(test[[1]]) <- 'VectorType'
Idents(test[[2]]) <- 'SendingType'
Idents(test[[3]]) <- 'ReceivingType'
Idents(test[[4]]) <- 'VectorType'

test[[1]] <- FindVariableFeatures(test[[1]],selection.method = "disp")
test[[1]] <- ScaleData(test[[1]])
test[[1]] <- RunPCA(test[[1]])
PCHeatmap(test[[1]],dims = 1:9,balanced = T,cells = 100)

# Experiment to see if neighbors have enrichment

table(Idents(test[[1]]))
table(Idents(test[[2]]))

samp1 <- subset(test[[1]],idents = 'Neural tube and notochord-Mesenchymal')
samp1$Sample <- 'NonSpatial'
samp2 <- subset(test[[2]],idents = 'Neural tube and notochord')
samp2$Sample <- 'Spatial'

merge <- merge(samp1,samp2)
Idents(merge) <- merge[['Sample']]
merge <- NormalizeData(merge)
merge <- ScaleData(merge)
markers <- FindAllMarkers(merge,min.pct = 0,logfc.threshold = 0)

# Other tests

test <- RunNICHES(FFPE2,assay = "SCT",species = 'mouse',position.x = 'x',position.y = 'y',
               CellToCell = F,CellToSystem = T,SystemToCell = F,LR.database = 'fantom5',
               CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F) 

test <- RunNICHES(FFPE2,assay = "SCT",species = 'mouse',position.x = 'x',position.y = 'y',
                  CellToCell = F,CellToSystem = F,SystemToCell = T,LR.database = 'fantom5',
                  CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F) 

test <- RunNICHES(FFPE2,assay = "SCT",species = 'mouse',position.x = 'x',position.y = 'y',
               CellToCell = T,CellToSystem = T,SystemToCell = T,LR.database = 'fantom5',
               CellToCellSpatial = T,CellToNeighborhood = F,NeighborhoodToCell = F)

test <- RunNICHES(FFPE2,assay = "SCT",species = 'mouse',position.x = 'x',position.y = 'y',
                  CellToCell = F,CellToSystem = F,SystemToCell = F,LR.database = 'fantom5',
                  CellToCellSpatial = F,CellToNeighborhood = T,NeighborhoodToCell = F)

test <- RunNICHES(FFPE2,assay = "SCT",species = 'mouse',position.x = 'x',position.y = 'y',
                  CellToCell = F,CellToSystem = F,SystemToCell = F,LR.database = 'fantom5',
                  CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = T)

test <- RunNICHES(FFPE2,assay = "SCT",species = 'mouse',position.x = 'x',position.y = 'y',
                  CellToCell = T,CellToSystem = T,SystemToCell = T,LR.database = 'fantom5',
                  CellToCellSpatial = T,CellToNeighborhood = F,NeighborhoodToCell = F) 

