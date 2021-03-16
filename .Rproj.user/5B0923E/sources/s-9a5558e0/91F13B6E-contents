#### TODO: add parameter check functions, e.g. positions have to be valid

#### Currently in order to run pair_spatial
#### I need user to specify "position" in the metadata and I'll also check the validity
#### add the neighborhood_radius parameter specifically for pair_spatial and niche_spatial

#' An S4 class to represent SCC hierarchy.
#'
#' @slot organizations A list containing different organizations
#' @export
setClass("SCC", representation(organizations="list"))


#' RunSCC
#'
#' Single-Cell Connectome (SCC) is a framework
#' that builds cell-cell ligand-receptor interactoin profiles
#' based on (spatial) single-cell transcriptomic data.
#'
#'
#' @param seu_obj input Seurat Object
#' @param organizations character vector, default: c("pair","pair_spatial","niche_spatial")
#' The type of SCC organizations to build
#' @param assay character, default: NULL
#' The assay to use from the input Seurat Object, NULL leads to the default assay
#' @param metadata character vector, default: NULL
#' The metadata to aggregate: e.g. cell_types, positions
#' @param species character, default: "human"
#' the species of the input scRNA-seq data
#' @param downsampling logical, default: TRUE
#' whether to do dowsampling on the cells
#' @param sample_size_n integer, default: 1000
#' the number of sampled cells if do downsampling
#' @param autocrine logical, default:TRUE
#' whether to consider cell interactions within one cell type
#' (including interact with themselves),always TRUE when no
#' cell_types provided
#' @param neighborhood_radius float, default:1
#' The euclidean distance threshold to consider as neighbors
#' @param n_threads integer,optional, default:8
#' the number of threads to use when constructing the matrix
#'
#' @return an SCC object which includes the specified organizations
#' @export
#'
#' @examples
#'
#' scc_object <- runSCC(seu_obj = spatial_data,organizations=c("pair","pair_spatial","niche_spatial"),
#'  assay='alra',metadata = c("cell_types","position","cluster"),neighborhood_radius = 3,
#'  species = "mouse",downsampling = FALSE,sample_size_n = 10,autocrine =TRUE,n_threads = 8)
#'
#'
runSCC <- function(
  seu_obj,
  organizations=c("pair","pair_spatial","niche_spatial"),
  assay=NULL,
  metadata=NULL,
  species="human",
  downsampling = TRUE,
  sample_size_n=1000,
  autocrine=TRUE,

  neighborhood_radius = 1,

  n_threads=8){
  #### Organizations can't be empty or contain strings other than "pair","pair_spatial","niche_spatial"
  if(is.null(organizations) | ((!"pair" %in% organizations)&(!"pair_spatial" %in% organizations)&(!"niche_spatial" %in% organizations))){
    stop("Error: Organization format isn't correct, which shouldn't be empty or contains organizations\
    other than 'pair','pair_spatial','niche_spatial'")
  }
  #### Return an Seurat object containing the SCC base matrix
  scc_base_obj <- scc_base(seu_obj,assay,metadata,species,downsampling,sample_size_n,autocrine,n_threads)

  #### Based on the input organizations, switch among cases to build SCC_obj
  #### double check if this the right way (if there is a better way) to set up class in the package

  org_list <- vector(mode="list",length = 3)
  #org_names_token <- c(FALSE,FALSE,FALSE)
  #names(org_names_token) <- c("pair","pair_spatial","niche_spatial")

  if("pair" %in% organizations){
    org1 <- scc_base_obj
    org_list[[1]] <-org1
    names(org_list)[1] <- "pair"
    #org_names_token['pair'] <- TRUE
  }
  if("pair_spatial" %in% organizations){
    #### only further requirement is the position input as meta data for each cell
    if(!("position" %in% metadata)) stop("Unable to run pair_spatial SCC because position isn't in metadata")
    position_aggregate <- as.vector(scc_base_obj@meta.data$position_aggregate)
    direct_neighbor_token <- sapply(position_aggregate, function(x){
      pos1 <- strsplit(x,"\\|")[[1]][1]
      pos2 <- strsplit(x,"\\|")[[1]][2]
      pos1 <- c(as.numeric(strsplit(pos1,"x")[[1]][1]),as.numeric(strsplit(pos1,"x")[[1]][2]))
      pos2 <- c(as.numeric(strsplit(pos2,"x")[[1]][1]),as.numeric(strsplit(pos2,"x")[[1]][2]))
      if(sqrt(sum((pos1 - pos2) ^ 2)) <= neighborhood_radius){
        return(TRUE)
      }
      else{return(FALSE)}
    })
    names(direct_neighbor_token) <- colnames(position_aggregate)
    org2 <- subset(scc_base_obj,cells = colnames(scc_base_obj)[direct_neighbor_token])
    #org_names_token['pair_spatial'] <- TRUE
    org_list[[2]] <-org2
    names(org_list)[2] <- "pair_spatial"

  }
  if("niche_spatial" %in% organizations){
    #### still need to run the spatial organization
    if(!("position" %in% metadata)) stop("Unable to run pair_spatial SCC because position isn't in metadata")
    position_aggregate <- as.vector(scc_base_obj@meta.data$position_aggregate)
    direct_neighbor_token <- sapply(position_aggregate, function(x){
      pos1 <- strsplit(x,"\\|")[[1]][1]
      pos2 <- strsplit(x,"\\|")[[1]][2]
      pos1 <- c(as.numeric(strsplit(pos1,"x")[[1]][1]),as.numeric(strsplit(pos1,"x")[[1]][2]))
      pos2 <- c(as.numeric(strsplit(pos2,"x")[[1]][1]),as.numeric(strsplit(pos2,"x")[[1]][2]))
      if(sqrt(sum((pos1 - pos2) ^ 2)) <= neighborhood_radius){
        return(TRUE)
      }
      else{return(FALSE)}
    })
    names(direct_neighbor_token) <- colnames(position_aggregate)

    spatial_org_tmp <- subset(scc_base_obj,cells = colnames(scc_base_obj)[direct_neighbor_token])
    niche_mat <- rowsum(t(as.matrix(spatial_org_tmp@assays$RNA@data)), group=spatial_org_tmp@meta.data$receiver_cell)

    ### what other metadata should I aggregate
    input_meta = NULL
    if("cell_types" %in% metadata){
        cell_type_mapper <- as.character(seu_obj@meta.data$cell_types)
        names(cell_type_mapper) <- colnames(seu_obj)

        input_meta <- as.data.frame(cell_type_mapper[rownames(niche_mat)])
        colnames(input_meta) <- c("cell_types")
    }
    org3 <- Seurat::CreateSeuratObject(counts=t(niche_mat),meta.data=input_meta)
    #org_names_token['niche_spatial'] <- TRUE
    org_list[[3]] <-org3
    names(org_list)[3] <- "niche_spatial"

  }

  SCC_obj <- new("SCC",organizations =org_list)

  return(SCC_obj)
}
