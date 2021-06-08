require(SCC)
require(Seurat)
require(ggplot2)

# Load Data
load("~/Box Sync/Science Response/Base objects for computational work (Use Human GEN2 for revision work)/objects/rat_cca 2018-08-04 .Robj")

# Upate
rat <- UpdateSeuratObject(rat_cca)
rat <- subset(rat,idents = c('Sox9+','Fib_Col13a1+','ATI','ATII','EC_cap','EC_vasc','Mac_alv'))

# Testing runs
test <- RunSCC(object = rat,LR.database = 'fantom5',species = 'rat',
               CellToCell = T,CellToSystem = T,SystemToCell = T,
               CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F,meta.data.to.map = c('nCount_RNA','exp')) #works
# Cluster
Idents(test[[1]]) <- 'VectorType'
Idents(test[[2]]) <- 'SendingType'
Idents(test[[3]]) <- 'ReceivingType'

test[[1]] <- FindVariableFeatures(test[[1]])
test[[1]] <- ScaleData(test[[1]])
test[[1]] <- RunPCA(test[[1]])
PCHeatmap(test[[1]],dims = 1:9,balanced = T,cells = 100)
test[[1]] <- RunUMAP(test[[1]],dims = 1:20)
DimPlot(test[[1]])+NoLegend()
Idents(test[[1]]) <- 'SendingType'
table(Idents(test[[1]]))
sub <- subset(test[[1]],idents = c('Fib_Col13a1+'))
Idents(sub) <- 'VectorType'
table(Idents(sub))
mark <- FindAllMarkers(sub,only.pos = T,min.pct = 0.5)
              