#' Load ligand-receptor annotations
#'
#' Load ligand-receptor annotations for different species
#'
#' @param species character, default: "human"
#' accepted species are "human","mouse","rat","pig"
#' @return a dataframe containing ligand-receptor from fantom5
#'
#' @export
#'
####load the ligands receptors of the specified species
load_lr_anno <- function(species="human"){
  switch(species,
         human = {return(get(data(ncomms8866_human)))},
         mouse = {return(get(data(ncomms8866_mouse)))},
         rat = {return(get(data(ncomms8866_rat)))},
         pig = {return(get(data(ncomms8866_pig)))},
         stop("input species not recognized"))
}
#' Geosketch downsampling
#'
#' Downsample the cells using a published method
#' called Geometric Sketching (details see Hie, Brian, et al.
#'  "Geometric sketching compactly summarizes the
#'  single-cell transcriptomic landscape." Cell systems 8.6 (2019): 483-493.)
#' @param exp_mat, scRNA-seq expression matrix
#' @param sample_size_n integer, optional, default: 1000
#' the number of sampled cells if do downsampling
#' @return
#' an index of the sampled cells
#' @export
#'
geo_downsample <- function(exp_mat,sample_size_n=1000){
  geosketch <- reticulate::import('geosketch')
  s <- rsvd::rsvd(t(exp_mat), k=10)
  X.pcs <- s$u %*% diag(s$d)
  sketch.size <- as.integer(sample_size_n)
  return(unlist(geosketch$gs(X.pcs, sketch.size,one_indexed=TRUE)))
}
#' Calculating column positions
#'
#' Calculating the column positions of the cell-pair in the clif matrix
#' based on their indices
#'
#' @param x integer
#' the index of the cell expressing the ligand gene
#' @param y integer
#' the index of the cell expressing the receptor gene
#' @param n integer
#' the total number of (downsampled) cells
#'
#' @return
#' the column position of the cell-pair
#'
#' @export
#'
###this func distinguishes cell1,cell2 and cell2,cell1
cal_pos_coldirection <- function(x,y,n){
  return((x-1)*n+y)
}
#' Calculating column positions in a vectorized way
#'
#' Calculating the column positions of the cell-pair in the clif matrix
#' in a vectorized way based on their indices
#'
#' @param x integer vector
#' the indices of the cells expressing the ligand gene
#' @param y integer vector
#' the indices of the cells expressing the receptor gene
#' @param n integer
#' the total number of (downsampled) cells
#'
#' @return
#' a vector of the column positions of the cell-pairs
#'
#' @export
#'
cal_pos_coldirection_vec <- Vectorize(cal_pos_coldirection,vectorize.args=c("x","y"))

#' Building the clif matrix
#'
#' building the clif matrix where rows are
#' ligand-recetor pairs, columns are cell-pairs,
#' entries are the expression level
#'
#' @param X
#' an initialized clif matrix
#' @param ind integer
#' index of a certain ligand-receptor pair
#' @param real_pairs character vector
#' a vector of ligand-receptor pairs
#' @param exp_data
#' the scRNA-seq expression matrix
#' @param ligand_cell_list integer list
#' a list of cell indices that express the ligands
#' @param ligand_index integer vector
#' vector of the indices of the ligands
#' @param recepts_index integer vector
#' vector of the indices of the receptors
#' @param receptor_cell_list integer list
#' a list of cell indices that express the receptors
#' @param gene_name_index integer vector
#' vector of the indices of the genes
#' @param cal_pos_coldirection_vec function
#' calculating column positions in a vectorized way
#' @param cell_names character vector
#' vector of the (downsampled cell names)
#'
#' @export
#'
builder <- function(X,ind,real_pairs,exp_data,ligand_cell_list,ligand_index,recepts_index,receptor_cell_list,gene_name_index,cal_pos_coldirection_vec,cell_names)
{

  LR <- unlist(strsplit(real_pairs[ind],split = "_"))

  ligan <- LR[1]
  recp <- LR[2]

  L_cells <- ligand_cell_list[[ligand_index[ligan]]]
  R_cells <- receptor_cell_list[[recepts_index[recp]]]
  ####only calculate the cell pairs that express both ligand and the receptor

  if(!is.null(L_cells)&!is.null(R_cells)){
    L_expressions <- as.vector(exp_data[gene_name_index[ligan],L_cells])
    R_expressions <- as.vector(exp_data[gene_name_index[recp],R_cells])
    L_length <- length(L_expressions)
    R_length <- length(R_expressions)
    elong_R_expressions <- c(rep(NA,L_length),R_expressions,rep(NA,L_length))
    L_index <- as.vector(L_cells)
    R_index <- as.vector(R_cells)
    elong_R_index <-  c(rep(NA,L_length),R_index,rep(NA,L_length))

    for(j in 1:(L_length+R_length-1)){

      ###the other operation is to take multiplication: NA is still preserved
      mul_vec <- L_expressions*elong_R_expressions[(j+1):(j+L_length)]
      pos_vec <- cal_pos_coldirection_vec(L_index,elong_R_index[(j+1):(j+L_length)],n = length(cell_names))

      X[ind,pos_vec[!is.na(pos_vec)]] <- X[ind,pos_vec[!is.na(pos_vec)]] + mul_vec[!is.na(mul_vec)]
    }

  }

  NULL
}
#' Vectorized builder the clif matrix
#'
#' building the clif matrix in a vectorized way
#' (vectorize the ligand-receptor pairs) to run in
#' parallel
#'
#' @param builder function
#' base function to build the clif matrix
#' @param X
#' an initialized clif matrix
#' @param ind integer vector
#' indices of ligand-receptor pairs
#' @param real_pairs character vector
#' a vector of ligand-receptor pairs
#' @param exp_data
#' the scRNA-seq expression matrix
#' @param ligand_cell_list integer list
#' a list of cell indices that express the ligands
#' @param ligand_index integer vector
#' vector of the indices of the ligands
#' @param recepts_index integer vector
#' vector of the indices of the receptors
#' @param receptor_cell_list integer list
#' a list of cell indices that express the receptors
#' @param gene_name_index integer vector
#' vector of the indices of the genes
#' @param cal_pos_coldirection_vec function
#' calculating column positions in a vectorized way
#' @param cell_names character vector
#' vector of the (downsampled cell names)
#'
#' @export
#'
coldirec_builder_vec <- Vectorize(builder,vectorize.args = c("ind"))
