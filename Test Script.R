# Test script

# Load data
table(Idents(panc8))

# Set up custom list
custom.list <- Connectome::ncomms8866_human
custom.list <- custom.list[,c(2,4)]

# Test SCC
test <- RunSCC(panc8,species = 'human',min.cells.per.ident = 1000)

test.imp <- RunSCC(panc8,species = 'human',min.cells.per.ident = 1000,assay = 'alra')

test.imp <- RunSCC(panc8,LR.database = custom.list,species = 'human',min.cells.per.ident = 1000,assay = 'alra')


# Test Milieu
test <- RunMilieu(panc8,species = 'human',min.cells.per.ident = 1000)

test.imp <- RunMilieu(panc8,species = 'human',min.cells.per.ident = 1000,assay = 'alra')

test.imp <- RunMilieu(panc8,LR.database = custom.list,species = 'human',min.cells.per.ident = 1000,assay = 'alra')

