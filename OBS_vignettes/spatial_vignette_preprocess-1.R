# Wrapper function to reproduce the processed 10x visium mouse anterior brain data
# based on the result in https://satijalab.org/seurat/articles/spatial_vignette.html ()

spatial_vignette_preprocessing <- function(){
  InstallData("stxBrain",force.reinstall = T)
  brain <- LoadData("stxBrain", type = "anterior1")
  
  brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
  
  brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
  brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
  brain <- FindClusters(brain, verbose = FALSE)
  
  brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
 
  cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
  cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
  cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
  cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)
  
  allen_reference <- readRDS("/Users/mbr29/Box Sync/Kluger_Lab/sConnectome Paper/Figure Drafting/visium")
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
  
  # Assign the most probable cell type as the label
  cell_type_pos <- sapply(c(1:ncol(cortex@assays$predictions)), function(x) which.max(cortex@assays$predictions[1:23,x]))
  cell_type_names <- rownames(cortex@assays$predictions)[1:23]
  names(cell_type_names) <- c(1:23)
  cell_type_spatial <- cell_type_names[cell_type_pos]
  
  cortex@meta.data$cell_type <- as.factor(cell_type_spatial)
  Idents(cortex) <- cortex@meta.data$cell_type
  
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  cell_unique_colors <- c(
    "mediumpurple","yellowgreen","darksalmon","slateblue4","aquamarine4","dodgerblue1","firebrick4","ivory3","royalblue3","red","seashell")
  tmp_cell_unique_color <- cell_unique_colors
  for(ind in 1:length(cell_unique_colors)){
    tmp_cell_unique_color[ind] <- gplots::col2hex(cell_unique_colors[ind])
  }
  

  return(cortex)
}
