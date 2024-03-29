#' Load FANTOM5 as ground truth
#'
#' @param species species reference in the fantom 5 database
#' @export

LoadFantom5 <- function(species){
  # jc: to get rid of the global variable check
  #utils::globalVariables(c("ncomms8866_human", "ncomms8866_mouse", "ncomms8866_rat","ncomms8866_pig"))
  if(is.null(species)){
    stop("\nPlease select species for FANTOM5 mapping. Allows 'human','mouse','rat', or 'pig' ")}
  # Load ground-truth database (FANTOM5, species-converted as appropriate, per methodlogy in Raredon et al 2019, DOI: 10.1126/sciadv.aaw3851)
  else{
    switch(species,
           human = {fantom <- get(utils::data(ncomms8866_human))},
           mouse = {fantom <- get(utils::data(ncomms8866_mouse))},
           rat = {fantom <- get(utils::data(ncomms8866_rat))},
           pig = {fantom <- get(utils::data(ncomms8866_pig))},
           stop("input species not recognized"))
  }
  
  fantom$mechanism <- paste(fantom$Ligand.ApprovedSymbol,fantom$Receptor.ApprovedSymbol,sep = '-')
  source.subunits <- as.matrix(data.frame(source_1 = fantom$Ligand.ApprovedSymbol)) #allows duplicate rownames
  rownames(source.subunits) <- fantom$Ligand.ApprovedSymbol
  target.subunits <- as.matrix(data.frame(target_1 = fantom$Receptor.ApprovedSymbol)) #allows duplicate rownames
  rownames(target.subunits) <- fantom$Receptor.ApprovedSymbol
  
  ground.truth <- list('source.subunits' = source.subunits,
                       'target.subunits' = target.subunits)
  return(ground.truth)
}
