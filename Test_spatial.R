devtools::document()
devtools::check()
devtools::install()


# Scripts from Seurat vignette
library(Seurat)
packageVersion("Seurat") '4.0.3'
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
install.packages("R.utils")
remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
install.packages("viridis")
library(viridis)

library("NICHES")

InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")

# Normalization 
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
# Dim reduction with all cells
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)

# Only use the cortex region data
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# Removal of some cells
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

# Integration
allen_reference <- readRDS("/data/github_proj/cell_interaction/spatial_experiments/seurat/allen_cortex.rds")
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30)
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE, 
                                  weight.reduction = cortex[["pca"]],dims = 1:30)
cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions"
cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "markvariogram", 
                                        features = rownames(cortex), r.metric = 5, slot = "data")

#################################################################################
#### My step: assign each cell with the most likely cell type prediction ######
cell_type_pos <- sapply(c(1:ncol(cortex@assays$predictions)), function(x) which.max(cortex@assays$predictions[1:23,x]))
cell_type_names <- rownames(cortex@assays$predictions)[1:23]
names(cell_type_names) <- c(1:23)
cell_type_spatial <- cell_type_names[cell_type_pos]

cortex@meta.data$cell_type <- as.factor(cell_type_spatial)
Idents(cortex) <- cortex@meta.data$cell_type

cell_unique_colors <- c(
  "mediumpurple","yellowgreen","darksalmon","slateblue4","aquamarine4","dodgerblue1","firebrick4","ivory3","royalblue3","red","seashell")
tmp_cell_unique_color <- cell_unique_colors
for(ind in 1:length(cell_unique_colors)){
  tmp_cell_unique_color[ind] <- gplots::col2hex(cell_unique_colors[ind])
}
DefaultAssay(cortex) <- "SCT"
cortex <- RunTSNE(cortex)

########
### Run SCC on Imputed Slot
#######

cortex@meta.data$x <- cortex@images$anterior1@coordinates$row
cortex@meta.data$y <- cortex@images$anterior1@coordinates$col

#DefaultAssay(cortex) <- "SCT"
DefaultAssay(cortex) <- "Spatial"
cortex <- NormalizeData(cortex)

cortex <- RunALRA(cortex)


cell_type_colors <- c("darkorchid2","green4","mediumspringgreen","#8DA0CB",
                      "orangered3","dodgerblue4","burlywood3","darkblue",
                      "cyan2" , "maroon4")
names(cell_type_colors) <- levels(cortex@meta.data$cell_type)
SpatialDimPlot(cortex, crop = TRUE, label = F,pt.size.factor = 2.2,cols = cell_type_colors)+NoLegend()

library(cowplot)
get_legend(SpatialDimPlot(cortex, crop = TRUE, label = F,pt.size.factor = 10,cols = cell_type_colors)
                        +guides(fill = guide_legend(title="Cell types",ncol = 1,override.aes = list(size = 4))))


cortex@assays$alra@counts <- cortex@assays$Spatial@counts
niche_obj <- RunNICHES(object = cortex,LR.database = "fantom5",species = "mouse",assay = "alra",
                       position.x = 'x', 
                       position.y = 'y',
                       rad.set = 2, #I think this is correct for the geometry
                       min.cells.per.ident = 1, # 1 will return error
                       min.cells.per.gene = 10,
                       meta.data.to.map = c('orig.ident','cell_type'),
                       CellToCell = F,CellToSystem = F,SystemToCell = F,
                       CellToCellSpatial = T,CellToNeighborhood = F,NeighborhoodToCell = T)


pair_scc <- niche_obj[[1]]
pair_scc <- ScaleData(pair_scc)
pair_scc <- FindVariableFeatures(pair_scc,selection.method = "disp")
pair_scc <- RunPCA(pair_scc)
pair_scc <- RunTSNE(
  pair_scc, tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:30, perplexity = 100,
  fast_tsne_path="/data/software/FIt-SNE/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)
Idents(pair_scc) <- "cell_type.Joint"

#Limit analysis to populations > 100 measurements
POI <- names(table(Idents(pair_scc))[table(Idents(pair_scc))>100])
pair_scc <- subset(pair_scc,idents = POI)
set.seed(31)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
sampled_colors_paircells <- sample(color, length(POI))

DimPlot(pair_scc,reduction = 'tsne',pt.size=0.5,shuffle = T,cols=sampled_colors_paircells) +ggtitle('Cell-Cell Signaling')+NoLegend()

legend_pair_tsne <- get_legend(DimPlot(pair_scc,reduction = 'tsne',pt.size=0.5,shuffle = T,cols=sampled_colors_paircells) +ggtitle('Cell-Cell Signaling')
                               + guides(col = guide_legend(title = "Type",nrow = 6,override.aes = list(size = 4)))     
                               +theme(legend.background = element_rect(color = NA),legend.position = "bottom")
                               + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))

plot_grid(legend_pair_tsne)


niche <- niche_obj[[2]]
Idents(niche) <- niche[['ReceivingType']]
#Limit analysis to populations > 50 measurements
POI_niche <- names(table(Idents(niche))[table(Idents(niche))>50])
niche <- subset(niche,idents = POI_niche)
# Scale and visualize
niche <- ScaleData(niche)
niche <- FindVariableFeatures(niche,selection.method = "disp")
niche <- RunPCA(niche)
niche <- RunTSNE(
  niche, tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:30, perplexity = 100,
  fast_tsne_path="/data/software/FIt-SNE/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)

DimPlot(niche,reduction = 'tsne',pt.size = 0.5,shuffle = T,cols = cell_type_colors) +ggtitle('Niche Signaling')+NoLegend()

legend_niche_tsne <- get_legend(DimPlot(niche,reduction = 'tsne',pt.size = 0.5,shuffle = T,cols = cell_type_colors) +ggtitle('Niche Signaling')
                                + guides(col = guide_legend(nrow = 3,override.aes = list(size = 4)))     
                                +theme(legend.background = element_rect(color = NA),legend.position = "bottom")
                                + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))

plot_grid(legend_niche_tsne)



# Find markers
mark <- FindAllMarkers(niche,min.pct = 0.25,only.pos = T,test.use = "roc")
#mark$ratio <- mark$pct.1/mark$pct.2
GOI_niche <- mark %>% group_by(cluster) %>% top_n(10,myAUC)

niche@meta.data$POI_receiving_type <- factor(niche@meta.data$ReceivingType,levels = POI_niche)
Idents(niche) <- niche@meta.data$POI_receiving_type
DoHeatmap(niche,features = unique(GOI_niche$gene),group.by = "POI_receiving_type") + scale_fill_gradientn(colors = c("grey","white", "blue"))

DefaultAssay(cortex) <- 'alra'
SpatialFeaturePlot(cortex, crop = TRUE, features = c("Fgf1","Fgfr2","Gnai2","Unc5b"),slot = "data",min.cutoff =  'q1',
                   max.cutoff = 'q99')+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
