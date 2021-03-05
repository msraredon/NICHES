

test <- RunSCC(FFPE2,species = 'mouse',position.x = 'x',position.y = 'y',
               CellCell = T,CellSystem = F,SystemCell = F,
               CellCellSpatial = F,CellNeighborhood = F,NeighborhoodCell = F) #works
test <- RunSCC(FFPE2,species = 'mouse',position.x = 'x',position.y = 'y',
               CellCell = F,CellSystem = T,SystemCell = F,
               CellCellSpatial = F,CellNeighborhood = F,NeighborhoodCell = F) #works
test <- RunSCC(FFPE2,species = 'mouse',position.x = 'x',position.y = 'y',
               CellCell = F,CellSystem = F,SystemCell = T,
               CellCellSpatial = F,CellNeighborhood = F,NeighborhoodCell = F) #works
test <- RunSCC(FFPE2,species = 'mouse',position.x = 'x',position.y = 'y',
               CellCell = F,CellSystem = F,SystemCell = F,
               CellCellSpatial = T,CellNeighborhood = F,NeighborhoodCell = F) #works
test <- RunSCC(FFPE2,species = 'mouse',position.x = 'x',position.y = 'y',
               CellCell = F,CellSystem = F,SystemCell = F,
               CellCellSpatial = F,CellNeighborhood = T,NeighborhoodCell = F)
test <- RunSCC(FFPE2,species = 'mouse',position.x = 'x',position.y = 'y',
               CellCell = F,CellSystem = F,SystemCell = F,
               CellCellSpatial = F,CellNeighborhood = F,NeighborhoodCell = T)
