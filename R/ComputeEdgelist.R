
#' Compute a node-to-node edgelist based on the spatial coordinates
#'
#' @param node.object A filtered Seurat object. The active identity will be used to define populations for connectomic sampling and crossings.
#' @param position.x string. Optional. Default: NULL. The name of the meta.data column specifying location on the spatial x-axis. Only required for spatial omics data.
#' @param position.y string. Optional. Default: NULL. The name of the meta.data column specifying location on the spatial y-axis. Only required for spatial omics data.
#' @param k integer. Optional. Default: 4. Number of neighbors in a knn graph. Used to compute a mutual nearest neighbor graph based on the spatial coordinates of the spatial transcriptomic datasets.  
#' @param rad.set numeric. Optional. Default: NULL. The radius threshold to define neighbors based on the spatial coordinates of the spatial transcriptomic datasets. Ignored when 'k' is provided.
#' @param nn.methods string. Optional. Default: NULL. Method to define nearest neighbors. If NULL, defaults to legacy technique, which is inefficient for very large datasets due to full edgelist construction followed by downsampling.
#' @export
#'
ComputeEdgelist <- function(node.object,
                             position.x,
                             position.y,
                             k=4,
                             rad.set=NULL,
                             nn.method=NULL){
  
  ### CREATE MAPPING ###
  if(is.null(nn.method)){ #### MSBR 2024-06-16
    
    
    # Create adjacency matrix
    # Adapted from :: https://stackoverflow.com/questions/16075232/how-to-create-adjacency-matrix-from-grid-coordinates-in-r
    # Setup numbering and labeling
    # jc: possible bug, change object to node.object
    df <- data.frame(x = node.object[[position.x]], y = node.object[[position.y]])
    df$barcode <- rownames(df)
    df$x <- as.character(df$x)
    df$y <- as.character(df$y)
    df$x <- as.numeric(df$x)
    df$y <- as.numeric(df$y)
    df <- df[,c('x','y')]
    
    # Compute the euclidean distance matrix
    distance_mat <- apply(df, 1, function(pt)
      (sqrt(abs(pt["x"] - df$x)^2 + abs(pt["y"] - df$y)^2))
    )
    # keep the cell names
    rownames(distance_mat) <- colnames(distance_mat)
    
    if(!is.null(k) & is.null(nn.method)){
      message("Compute edgelist based on mutual nearest neighbors.")
      if(!is.null(rad.set))warning("'k' is not NULL. Parameter 'rad.set' will be ignored.")
      # neighbor_mat: the indices of k nearest neighbors  
      neighbor_mat <- apply(distance_mat,1,function(dis_vec){
        order(dis_vec)[1:(k+1)]
      })
      
      # k nearest neighbor adjacency matrix
      adj_mat <- matrix(data=0,nrow=nrow(distance_mat),ncol = ncol(distance_mat))
      rownames(adj_mat) <- rownames(distance_mat)
      colnames(adj_mat) <- colnames(distance_mat)
      for(cell in colnames(neighbor_mat)) adj_mat[cell,neighbor_mat[,cell]] <- 1
      # mutual nearest neighbor adjacency matrix
      adj_mat_final <- 1*(adj_mat & t(adj_mat))
      
    }
    else if(!is.null(rad.set)){
      message("\n Compute edgelist based on spatial radius threshold.")
      # Alternatively the adjacency matrix is defined by the absolute radius
      adj_mat_final <- 1*(distance_mat <= rad.set)
    }
    else{
      stop("Both k and rad.set are NULL.")
    }
    
    # Convert adj matrix to edgelist
    edgelist <- igraph::graph.adjacency(adj_mat_final)
    edgelist <- igraph::get.data.frame(edgelist)
  }
  
  #### MSBR 2024-06-16
  #### Revised MSBR 2025-11-22 to fix k inheritance bug
  
  if(!(is.null(k)) & nn.method == 'aoz'){
    # df <- data.frame(x = node.object[[position.x]], y = node.object[[position.y]])
    # df$barcode <- rownames(df)
    # df$x <- as.character(df$x)
    # df$y <- as.character(df$y)
    # df$x <- as.numeric(df$x)
    # df$y <- as.numeric(df$y)
    # df <- df[,c('x','y')] 
    # 
    # coords <- as.matrix(df)#cbind(data.list[[i]]$x,data.list[[i]]$y)
    coords <- as.matrix(cbind(node.object$x,node.object$y))
    ord <- order(coords[,1])
    
    n.neighbors <- k # inherits input 'k'
    n.omp.t <- 1
    
    n <- nrow(coords)
    
    
    u.search.type <- 2 ##2 is very fast, 1 is slightly faster than 0, and 0 is the orginal slow one (2 and 0 should match, 1 is also corrected just different opposite sorting among the children)
    
    indx <- spNNGP:::mkNNIndxCB(coords, n.neighbors, n.omp.t)
    
    nn.indx <- indx$nnIndx
    nn.indx.lu <- indx$nnIndxLU
    nn.indx.run.time <- indx$run.time
    
    storage.mode(nn.indx) <- "integer"
    storage.mode(nn.indx.lu) <- "integer"
    
    n.indx <- spNNGP:::mk.n.indx.list(nn.indx, n, n.neighbors)
    names(n.indx) <- 1:nrow(coords)
    
    # build sparse Adjacency Matrix
    n.indx.names <- names(n.indx)
    # this unusual code let's me reference a list element's 
    # name INSIDE the apply function, neat
    # THIS NEXT STEP AOZ WROTE - IT'S SLOW AND COULD BE SPED UP
    adj.ij <- lapply(setNames(n.indx.names, n.indx.names), function(nameindex) {
      data.frame(i=rep(as.integer(nameindex), length(n.indx[[nameindex]])),
                 j=n.indx[[nameindex]])})
    adj.ij <- data.table::rbindlist(adj.ij)
    # clean and make it mutual
    adj.ij <- na.omit(adj.ij)
    adj.ij <- rbind(adj.ij, 
                    data.frame(i=adj.ij$j,
                               j=adj.ij$i))
    # map back to original coord indexing
    adj.ij$i <- ord[adj.ij$i]
    adj.ij$j <- ord[adj.ij$j]
    adj.ij <- as.data.frame(adj.ij)
    colnames(adj.ij) <- c("from", "to")
    edgelist <- adj.ij
    
    
  }
  return(edgelist)
  
}