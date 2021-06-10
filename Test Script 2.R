require(SCC)
require(Seurat)
require(ggplot2)

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
test <- RunSCC(object = FFPE2,LR.database = 'omnipath',species = 'mouse',position.x = 'x',position.y = 'y',
               CellToCell = T,CellToSystem = T,SystemToCell = T,
               CellToCellSpatial = T,CellToNeighborhood = T,NeighborhoodToCell = T,meta.data.to.map = c('nCount_SCT','orig.ident')) #works
# Cluster
Idents(test[[1]]) <- 'VectorType'
Idents(test[[2]]) <- 'SendingType'
Idents(test[[3]]) <- 'ReceivingType'

test[[1]] <- FindVariableFeatures(test[[1]])
test[[1]] <- ScaleData(test[[1]])
test[[1]] <- RunPCA(test[[1]])
PCHeatmap(test[[1]],dims = 1:9,balanced = T,cells = 100)

# Experiment to see if neighbors have enrichment

table(Idents(test[[1]]))
table(Idents(test[[2]]))

samp1 <- subset(test[[1]],idents = 'Neural tube and notochord-Mesenchymal')
samp1$Sample <- 'NonSpatial'
samp2 <- subset(test[[2]],idents = 'Neural tube and notochord-Mesenchymal')
samp2$Sample <- 'Spatial'

merge <- merge(samp1,samp2)
Idents(merge) <- merge[['Sample']]
merge <- NormalizeData(merge)
merge <- ScaleData(merge)
markers <- FindAllMarkers(merge,min.pct = 0,logfc.threshold = 0)

# Other tests

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

# Look at distributions!
setwd("~/Box Sync/Kluger_Lab")
orig <- data.frame(prop.table(table(Idents(FFPE2))))
Idents(test[[1]]) <- 'VectorType'
autocrine <- test[[1]][,test[[1]]$SendingType == test[[1]]$ReceivingType]
nonautocrine <- test[[1]][,!(test[[1]]$SendingType == test[[1]]$ReceivingType)]
auto <- data.frame(prop.table(table(autocrine$SendingType)))
sending.nonauto <- data.frame(prop.table(table(nonautocrine$SendingType)))
receiving.nonauto <- data.frame(prop.table(table(nonautocrine$ReceivingType)))
orig$`Fraction Of` <- 'Original Data'
auto$`Fraction Of` <- 'Autocrine Signaling'
sending.nonauto$`Fraction Of` <- 'Signaling (Outgoing)'
receiving.nonauto$`Fraction Of` <- 'Signaling (Incoming)'
data <- rbind(orig,auto,sending.nonauto,receiving.nonauto)
pdf('Distributions (CelltoCell).pdf',width = 5.5,height = 4)
ggplot(data,aes(x = Var1,y = Freq,fill = `Fraction Of`)) + 
  geom_bar(stat = 'identity',position = 'dodge')+ 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  xlab('Cell Type')+
  ylab('Fraction')+
  ggtitle('Distributions (CelltoCell)')
dev.off()

# Spatial!
orig <- data.frame(prop.table(table(Idents(FFPE2))))
Idents(test[[2]]) <- 'VectorType'
autocrine <- test[[2]][,test[[2]]$SendingType == test[[2]]$ReceivingType]
nonautocrine <- test[[2]][,!(test[[2]]$SendingType == test[[2]]$ReceivingType)]
auto <- data.frame(prop.table(table(autocrine$SendingType)))
sending.nonauto <- data.frame(prop.table(table(nonautocrine$SendingType)))
receiving.nonauto <- data.frame(prop.table(table(nonautocrine$ReceivingType)))
orig$`Fraction Of` <- 'Original Data'
auto$`Fraction Of` <- 'Autocrine Signaling'
sending.nonauto$`Fraction Of` <- 'Signaling (Outgoing)'
receiving.nonauto$`Fraction Of` <- 'Signaling (Incoming)'
data <- rbind(orig,auto,sending.nonauto,receiving.nonauto)
ggplot(data,aes(x = Var1,y = Freq,fill = `Fraction Of`)) + geom_bar(stat = 'identity',position = 'dodge')+ theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
pdf('Distributions (CelltoCellSpatial).pdf',width = 5.5,height = 4)
ggplot(data,aes(x = Var1,y = Freq,fill = `Fraction Of`)) + 
  geom_bar(stat = 'identity',position = 'dodge')+ 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  xlab('Cell Type')+
  ylab('Fraction')+
  ggtitle('Distributions (CelltoCellSpatial)')
dev.off()
