# devtools::document()
# devtools::check()
# devtools::install()

library(NICHES)

# Use Seurat data as examples

library(Seurat)
data1 <-  Read10X(data.dir = "/data/jyc/test_data/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = data1, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



pbmc@meta.data$cell_types <- Idents(pbmc)
niche_obj <- RunNICHES(object = pbmc,assay = "RNA",LR.database = "fantom5",species = "human",
                       min.cells.per.ident = 20, 
                       min.cells.per.gene = 100,
                       meta.data.to.map = colnames(pbmc@meta.data),
                       CellToCell = F,CellToSystem = T,SystemToCell = T,
                       CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F)

niche_obj1_tmp <- RunNICHES(object = pbmc,assay = "RNA",cell_types = "cell_types",LR.database = "fantom5",species = "human",
                       min.cells.per.ident = 20, 
                       min.cells.per.gene = 100,
                       meta.data.to.map = colnames(pbmc@meta.data),
                       CellToCell = F,CellToSystem = T,SystemToCell = T,
                       CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F)


niche_obj1_tmp <- RunNICHES(object = pbmc,assay = "asasas",LR.database = "fantom5",species = "human",
                       min.cells.per.ident = 20, 
                       min.cells.per.gene = 100,
                       meta.data.to.map = colnames(pbmc@meta.data),
                       CellToCell = T,CellToSystem = T,SystemToCell = T,
                       CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F)

niche_obj2 <- RunNICHES(object = as.matrix(pbmc@assays$RNA@data),LR.database = "fantom5",species = "human",
                        meta.data.df = pbmc@meta.data,
                        cell_types = "cell_types",
                        min.cells.per.ident = 20, 
                        min.cells.per.gene = 100,
                        CellToCell = T,CellToSystem = T,SystemToCell = T,
                        CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F)

niche_obj2_tmp <- RunNICHES(object = pbmc@assays$RNA@data,LR.database = "fantom5",species = "human",
                            meta.data.df = pbmc@meta.data,
                        
                        min.cells.per.ident = 20, 
                        min.cells.per.gene = 100,
                        CellToCell = F,CellToSystem = T,SystemToCell = T,
                        CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F)


niche_obj3 <- RunNICHES(object = as.matrix(pbmc@assays$RNA@data),LR.database = "fantom5",species = "human",
                        meta.data.df = pbmc@meta.data,
                        cell_types = "cell_types",
                        min.cells.per.ident = 20, 
                        min.cells.per.gene = 100,
                        CellToCell = F,CellToSystem = T,SystemToCell = T,
                        CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F)



niche_obj4 <- RunNICHES(object = as.matrix(pbmc@assays$RNA@data),LR.database = "fantom5",species = "human",
                        meta.data.to.map = pbmc@meta.data,
                        cell_types = "cell_types",
                        min.cells.per.ident = 20, 
                        min.cells.per.gene = 100,
                        CellToCell = T,CellToSystem = T,SystemToCell = T,
                        CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F,output_format = "raw")


niche_obj5 <- RunNICHES(object = pbmc,assay = "RNA",LR.database = "fantom5",species = "human",
                       min.cells.per.ident = 20, 
                       min.cells.per.gene = 100,
                       meta.data.to.map = colnames(pbmc@meta.data),
                       CellToCell = T,CellToSystem = T,SystemToCell = T,
                       CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F,output_format = "raw")


custom_db <- data.frame("L" = ncomms8866_human$Ligand.ApprovedSymbol[1:100],"R" = ncomms8866_human$Receptor.ApprovedSymbol[1:100])

niche_obj2 <- RunNICHES(object = pbmc,assay = "RNA",LR.database = "custom",species = "human",custom_LR_database=custom_db,
                       min.cells.per.ident = 20, 
                       min.cells.per.gene = 100,
                       meta.data.to.map = NULL,
                       CellToCell = T,CellToSystem = T,SystemToCell = T,
                       CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F)

custom_db3 <- data.frame("L" = ncomms8866_human$Pair.Name,"R" = ncomms8866_human$Receptor.ApprovedSymbol)

niche_obj3 <- RunNICHES(object = pbmc,assay = "RNA",LR.database = "custom",species = "human",custom_LR_database=custom_db3,
                        min.cells.per.ident = 20, 
                        min.cells.per.gene = 100,
                        meta.data.to.map = NULL,
                        CellToCell = T,CellToSystem = T,SystemToCell = T,
                        CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F)


# Visualize it

pairs <- niche_obj[["CellToCell"]]
pairs <- ScaleData(pairs)
pairs <- FindVariableFeatures(pairs,selection.method = "disp")
pairs <- RunPCA(pairs)
pairs <- RunTSNE(
  pairs, tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:30, perplexity = 100,
  fast_tsne_path="/data/software/FIt-SNE/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)
DimPlot(pairs,reduction = "tsne",group.by = "VectorType") + NoLegend()


niches <- niche_obj[['SystemToCell']]
niches <- ScaleData(niches)
niches <- RunPCA(niches,features = rownames(niches))
niches <- RunTSNE(
  niches, tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:30, perplexity = 100,
  fast_tsne_path="/data/software/FIt-SNE/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)
DimPlot(niches,reduction = "tsne",group.by = "ReceivingType")



