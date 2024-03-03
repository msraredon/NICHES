# JC: Test the scaling of different knn methods
# current method: brute-force
# target method: kd-tree search

library(Seurat)
library(Matrix)
devtools::load_all()
## first simulate a giant ST dataset
n_cells <- 400
n_genes <- 20000
# nb parameters
mu <- 1
size <- 20
# assuming a gaussian 
set.seed(42)

count_data <- as(matrix(rnbinom(n = (n_cells*n_genes),mu=mu,size=size),nrow = n_genes,ncol = n_cells),"sparseMatrix")
genes <- read.csv("/data/test_data/test_data/filtered_gene_bc_matrices/hg19/genes.tsv",sep = '\t',header = F)
rownames(count_data) <- unique(genes$V2)[1:n_genes]

st_data <- CreateSeuratObject(counts = count_data)

st_data@meta.data$x <- c(1:200)
st_data@meta.data$y <- c(1:200)
st_data@meta.data$cell_types <- c(rep("cp1",n_cells/2),rep("cp2",n_cells/2))
st_data <- NormalizeData(st_data)

NICHES_output <- RunNICHES(object = st_data,
                           LR.database = "fantom5",
                           species = "human",
                           assay = "RNA",
                           position.x = 'x',
                           position.y = 'y',
                           k = 4, 
                           cell_types = "cell_types",
                           min.cells.per.ident = 0,
                           min.cells.per.gene = NULL,
                           meta.data.to.map = NULL,
                           CellToCell = T,CellToSystem = F,SystemToCell = F,
                           CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = F)


