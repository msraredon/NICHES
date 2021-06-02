setwd("~/Documents/GitHub/SCC/R")

load("~/Box Sync/Science Response/Base objects for computational work (Use Human GEN2 for revision work)/objects/rat_cca 2018-08-04 .Robj")

rat_cca <- UpdateSeuratObject(rat_cca)

sending <- subset(rat_cca,idents = 'ATI',downsample = 200)
receiving <- subset(rat_cca,idents = 'EC_cap',downsample = 200)
DoHeatmap(sending,c('Vegfa'))
DoHeatmap(receiving,'Kdr')
downstream <- c('Map3k5','Raf1','Rock1','Rock2')
DoHeatmap(receiving,downstream)

m1 <- sending@assays[['RNA']]@data['Vegfa',]
m2 <- receiving@assays[['RNA']]@data['Kdr',]
m3 <- receiving@assays[['RNA']]@data[downstream,]
m1*m2*colMeans(m3)
receiving$test <- m1*m2*colMeans(m3)

