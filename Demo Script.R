# Set WD
setwd("~/Box Sync/Kluger_Lab/SCC Functions/SCC Milieu Contribution Demonstration")

# Demonstration of SCC
require(Seurat)
require(Connectome)
require(SCC)

# Load data
load("~/Box Sync/Science Response/Base objects for computational work (Use Human GEN2 for revision work)/objects/rat_cca 2018-08-04 .Robj")
# Update
rat <- UpdateSeuratObject(rat_cca)
# Define small system for demonstration
table(Idents(rat))
rat_sub <- subset(rat,idents = c('ATI','ATII','Mac_alv'))
rat_sub <- subset(rat_sub,downsample = 1000)
table(Idents(rat_sub))
# Impute and also run SCT, for options later
rat_sub <- RunALRA(rat_sub)
rat_sub <- SCTransform(rat_sub)

#RunSCC
scc.rna <- RunSCC(rat_sub,species = 'rat',assay = 'RNA')
scc.alra <- RunSCC(rat_sub,species = 'rat',assay = 'alra')
scc.sct <- RunSCC(rat_sub,species = 'rat',assay = 'SCT')
Idents(scc.rna) <- scc.rna[['VectorType']]
Idents(scc.alra) <- scc.alra[['VectorType']]
Idents(scc.sct) <- scc.sct[['VectorType']]

#RunMilieu
mil.rna.sum <- RunMilieu(rat_sub,species = 'rat',assay = 'RNA',blend = 'sum')
mil.alra.sum <- RunMilieu(rat_sub,species = 'rat',assay = 'alra',blend = 'sum')
mil.sct.sum <- RunMilieu(rat_sub,species = 'rat',assay = 'SCT',blend = 'sum')
mil.rna.mean <- RunMilieu(rat_sub,species = 'rat',assay = 'RNA',blend = 'mean')
mil.alra.mean <- RunMilieu(rat_sub,species = 'rat',assay = 'alra',blend = 'mean')
mil.sct.mean <- RunMilieu(rat_sub,species = 'rat',assay = 'SCT',blend = 'mean')
Idents(mil.rna.sum) <- mil.rna.sum[['ReceivingType']]
Idents(mil.alra.sum) <- mil.alra.sum[['ReceivingType']]
Idents(mil.sct.sum) <- mil.sct.sum[['ReceivingType']]
Idents(mil.rna.mean) <- mil.rna.mean[['ReceivingType']]
Idents(mil.alra.mean) <- mil.alra.mean[['ReceivingType']]
Idents(mil.sct.mean) <- mil.sct.mean[['ReceivingType']]

#RunContribution
cont.rna.sum <- RunContribution(rat_sub,species = 'rat',assay = 'RNA',blend = 'sum')
cont.alra.sum <- RunContribution(rat_sub,species = 'rat',assay = 'alra',blend = 'sum')
cont.sct.sum <- RunContribution(rat_sub,species = 'rat',assay = 'SCT',blend = 'sum')
cont.rna.mean <- RunContribution(rat_sub,species = 'rat',assay = 'RNA',blend = 'mean')
cont.alra.mean <- RunContribution(rat_sub,species = 'rat',assay = 'alra',blend = 'mean')
cont.sct.mean <- RunContribution(rat_sub,species = 'rat',assay = 'SCT',blend = 'mean')
Idents(cont.rna.sum) <- cont.rna.sum[['SendingType']]
Idents(cont.alra.sum) <- cont.alra.sum[['SendingType']]
Idents(cont.sct.sum) <- cont.sct.sum[['SendingType']]
Idents(cont.rna.mean) <- cont.rna.mean[['SendingType']]
Idents(cont.alra.mean) <- cont.alra.mean[['SendingType']]
Idents(cont.sct.mean) <- cont.sct.mean[['SendingType']]

# Source object embedding
rat_sub <- ScaleData(rat_sub)
rat_sub <- FindVariableFeatures(rat_sub)
rat_sub <- RunPCA(rat_sub)
ElbowPlot(rat_sub)
rat_sub <- RunUMAP(rat_sub,dims = 1:5)
pdf(file = 'Source Object.pdf',width = 5,height = 5)
UMAPPlot(rat_sub)+ggtitle('Source Object')
dev.off()

# SCC object embedding (RNA)
scc.rna <- ScaleData(scc.rna)
scc.rna <- FindVariableFeatures(scc.rna)
scc.rna <- RunPCA(scc.rna)
ElbowPlot(scc.rna)
scc.rna <- RunUMAP(scc.rna,dims = 1:10)
pdf(file = 'SCC.RNA Object.pdf',width = 5,height = 5)
UMAPPlot(scc.rna)+ggtitle('SCC on RNA')
dev.off()

# SCC object embedding (ALRA)
scc.alra <- ScaleData(scc.alra)
scc.alra <- FindVariableFeatures(scc.alra)
scc.alra <- RunPCA(scc.alra)
ElbowPlot(scc.alra)
scc.alra <- RunUMAP(scc.alra,dims = 1:15)
pdf(file = 'SCC.ALRA Object.pdf',width = 5,height = 5)
UMAPPlot(scc.alra)+ggtitle('SCC on ALRA')
dev.off()

# SCC object embedding (SCT)
scc.sct <- ScaleData(scc.sct)
scc.sct <- FindVariableFeatures(scc.sct)
scc.sct <- RunPCA(scc.sct)
ElbowPlot(scc.sct)
scc.sct <- RunUMAP(scc.sct,dims = 1:15)
pdf(file = 'SCC.SCT Object.pdf',width = 5,height = 5)
UMAPPlot(scc.sct)+ggtitle('SCC on SCT')
dev.off()

# Milieu object embedding (RNA,SUM)
mil.rna.sum <- ScaleData(mil.rna.sum)
mil.rna.sum <- FindVariableFeatures(mil.rna.sum)
mil.rna.sum <- RunPCA(mil.rna.sum)
ElbowPlot(mil.rna.sum)
mil.rna.sum <- RunUMAP(mil.rna.sum,dims = 1:6)
pdf(file = 'MIL.RNA.SUM Object.pdf',width = 5,height = 5)
UMAPPlot(mil.rna.sum)+ggtitle('Milieu on RNA')
dev.off()

# Milieu object embedding (ALRA,SUM)
mil.alra.sum <- ScaleData(mil.alra.sum)
mil.alra.sum <- FindVariableFeatures(mil.alra.sum)
mil.alra.sum <- RunPCA(mil.alra.sum)
ElbowPlot(mil.alra.sum)
mil.alra.sum <- RunUMAP(mil.alra.sum,dims = 1:6)
pdf(file = 'MIL.ALRA.SUM Object.pdf',width = 5,height = 5)
UMAPPlot(mil.alra.sum)+ggtitle('Milieu on ALRA')
dev.off()

# Milieu object embedding (SCT,SUM)
mil.sct.sum <- ScaleData(mil.sct.sum)
mil.sct.sum <- FindVariableFeatures(mil.sct.sum)
mil.sct.sum <- RunPCA(mil.sct.sum)
ElbowPlot(mil.sct.sum)
mil.sct.sum <- RunUMAP(mil.sct.sum,dims = 1:6)
pdf(file = 'MIL.SCT.SUM Object.pdf',width = 5,height = 5)
UMAPPlot(mil.sct.sum)+ggtitle('Milieu on SCT')
dev.off()

# Milieu object embedding (RNA,MEAN)
mil.rna.mean <- ScaleData(mil.rna.mean)
mil.rna.mean <- FindVariableFeatures(mil.rna.mean)
mil.rna.mean <- RunPCA(mil.rna.mean)
ElbowPlot(mil.rna.mean)
mil.rna.mean <- RunUMAP(mil.rna.mean,dims = 1:6)
pdf(file = 'MIL.RNA.MEAN Object.pdf',width = 5,height = 5)
UMAPPlot(mil.rna.mean)+ggtitle('Milieu on RNA w/ Mean')
dev.off()

# Milieu object embedding (ALRA,MEAN)
mil.alra.mean <- ScaleData(mil.alra.mean)
mil.alra.mean <- FindVariableFeatures(mil.alra.mean)
mil.alra.mean <- RunPCA(mil.alra.mean)
ElbowPlot(mil.alra.mean)
mil.alra.mean <- RunUMAP(mil.alra.mean,dims = 1:6)
pdf(file = 'MIL.ALRA.MEAN Object.pdf',width = 5,height = 5)
UMAPPlot(mil.alra.mean)+ggtitle('Milieu on ALRA w/ Mean')
dev.off()

# Milieu object embedding (SCT,MEAN)
mil.sct.mean <- ScaleData(mil.sct.mean)
mil.sct.mean <- FindVariableFeatures(mil.sct.mean)
mil.sct.mean <- RunPCA(mil.sct.mean)
ElbowPlot(mil.sct.mean)
mil.sct.mean <- RunUMAP(mil.sct.mean,dims = 1:6)
pdf(file = 'MIL.SCT.MEAN Object.pdf',width = 5,height = 5)
UMAPPlot(mil.sct.mean)+ggtitle('Milieu on SCT w/ Mean')
dev.off()


# Cont object embedding (RNA,SUM)
cont.rna.sum <- ScaleData(cont.rna.sum)
cont.rna.sum <- FindVariableFeatures(cont.rna.sum)
cont.rna.sum <- RunPCA(cont.rna.sum)
ElbowPlot(cont.rna.sum)
cont.rna.sum <- RunUMAP(cont.rna.sum,dims = 1:6)
pdf(file = 'CONT.RNA.SUM Object.pdf',width = 5,height = 5)
UMAPPlot(cont.rna.sum)+ggtitle('Contribution on RNA')
dev.off()

# Cont object embedding (ALRA,SUM)
cont.alra.sum <- ScaleData(cont.alra.sum)
cont.alra.sum <- FindVariableFeatures(cont.alra.sum)
cont.alra.sum <- RunPCA(cont.alra.sum)
ElbowPlot(cont.alra.sum)
cont.alra.sum <- RunUMAP(cont.alra.sum,dims = 1:6)
pdf(file = 'CONT.ALRA.SUM Object.pdf',width = 5,height = 5)
UMAPPlot(cont.alra.sum)+ggtitle('Contribution on ALRA')
dev.off()

# Cont object embedding (SCT,SUM)
cont.sct.sum <- ScaleData(cont.sct.sum)
cont.sct.sum <- FindVariableFeatures(cont.sct.sum)
cont.sct.sum <- RunPCA(cont.sct.sum)
ElbowPlot(cont.sct.sum)
cont.sct.sum <- RunUMAP(cont.sct.sum,dims = 1:6)
pdf(file = 'CONT.SCT.SUM Object.pdf',width = 5,height = 5)
UMAPPlot(cont.sct.sum)+ggtitle('Contribution on SCT')
dev.off()

# Cont object embedding (RNA,MEAN)
cont.rna.mean <- ScaleData(cont.rna.mean)
cont.rna.mean <- FindVariableFeatures(cont.rna.mean)
cont.rna.mean <- RunPCA(cont.rna.mean)
ElbowPlot(cont.rna.mean)
cont.rna.mean <- RunUMAP(cont.rna.mean,dims = 1:6)
pdf(file = 'CONT.RNA.MEAN Object.pdf',width = 5,height = 5)
UMAPPlot(cont.rna.mean)+ggtitle('Contribution on RNA w/ Mean')
dev.off()

# Cont object embedding (ALRA,MEAN)
cont.alra.mean <- ScaleData(cont.alra.mean)
cont.alra.mean <- FindVariableFeatures(cont.alra.mean)
cont.alra.mean <- RunPCA(cont.alra.mean)
ElbowPlot(cont.alra.mean)
cont.alra.mean <- RunUMAP(cont.alra.mean,dims = 1:6)
pdf(file = 'CONT.ALRA.MEAN Object.pdf',width = 5,height = 5)
UMAPPlot(cont.alra.mean)+ggtitle('Contribution on ALRA w/ Mean')
dev.off()

# Cont object embedding (SCT,MEAN)
cont.sct.mean <- ScaleData(cont.sct.mean)
cont.sct.mean <- FindVariableFeatures(cont.sct.mean)
cont.sct.mean <- RunPCA(cont.sct.mean)
ElbowPlot(cont.sct.mean)
cont.sct.mean <- RunUMAP(cont.sct.mean,dims = 1:6)
pdf(file = 'CONT.SCT.MEAN Object.pdf',width = 5,height = 5)
UMAPPlot(cont.sct.mean)+ggtitle('Contribution on SCT w/ Mean')
dev.off()
