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

test <- RunMilieu(panc8,species = 'human',min.cells.per.ident = 1000,blend = 'sum')
test <- RunMilieu(panc8,species = 'human',min.cells.per.ident = 1000,blend = 'mean')


# Test Contribution
test <- RunContribution(panc8,species = 'human')

test.imp <- RunContribution(panc8,species = 'human',assay = 'alra')

test.imp <- RunContribution(panc8,LR.database = custom.list,species = 'human',min.cells.per.ident = 1000,assay = 'alra')

test <- RunContribution(panc8,species = 'human',min.cells.per.ident = 1000,blend = 'sum')
test <- RunContribution(panc8,species = 'human',min.cells.per.ident = 1000,blend = 'mean')
