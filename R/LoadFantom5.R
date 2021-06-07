#' Load FANTOM5 as ground truth
#'
#' @export

LoadFantom5 <- function(species){
  
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
  source.subunits <- data.frame(source_1 = fantom$Ligand.ApprovedSymbol)
  target.subunits <- data.frame(target_1 = fantom$Receptor.ApprovedSymbol)
  ground.truth <- list('source.subunits' = source.subunits,
                       'target.subunits' = target.subunits)
}