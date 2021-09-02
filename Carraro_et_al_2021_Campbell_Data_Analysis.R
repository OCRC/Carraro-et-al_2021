# Analysis of Arcuate Nucleus and ME single-cell RNA sequencing data
# Original data by Campbell et al, 2017, Nature Neuroscience (https://doi.org/10.1038/nn.4495)
# Data analysis by Davi Sidarta-Oliveira as in Carraro et al., 2021

# Note 1: the final Seurat object is also available as a compressed .Rds file.
# Useful if you don't want to go through all the processing hassle and just want to 
# explore the data and check our results.

# Note 2: if you downloaded data from GEO, you should use the 'barcodes.txt'
# file provided with the code repository to guarantee you are using the exact
# same cells as we did when downloading from Single-Cell Portal

#Load libraries and setwd

reticulate::use_python('/usr/bin/python3')
library(loomR)
library(reticulate)
library(Seurat)
library(SeuratWrappers)
library(plotly)

setwd("~/Documents/Bioinfo/Carraro/")


######################################################################
# Load data, create Seurat object, QC plots
######################################################################

counts <- as.matrix(read.table('expression.txt.gz', sep = '\t', header = T, row.names = 1))
meta <- read.table('meta.txt.gz', sep = '\t', header = T, row.names = 1)
dat <- CreateSeuratObject(counts = counts, meta.data = meta)

# If downloaded from GEO:
barcodes <- read.table('barcodes.txt')
dat <- subset(dat, cells = barcodes)

mito.genes <- grep(pattern = "^mt-", x = rownames(dat@assays$RNA@counts), value = TRUE)
percent.mito <- Matrix::colSums(dat@assays$RNA@counts[mito.genes, ])/Matrix::colSums(dat@assays$RNA@counts)
dat <- AddMetaData(object = dat, metadata = percent.mito, col.name = "percent.mito")

VlnPlot(object = dat, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"))


######################################################################
# Default SCTransform workflow; PCA and UMAP for control
######################################################################

dat <- SCTransform(dat, return.only.var.genes = F, variable.features.n = 5000)

dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat, ndims = 100)

dat <- RunUMAP(dat, dims = 1:50)
UMAPPlot(dat, group.by = 'All.Cell.Subclusters')


###############################################################################
# Run the diffusion basis for dbMAP
###############################################################################

# Import python libraries
library(reticulate)
np <- reticulate::import("numpy")
pd <- reticulate::import("pandas")
sp <- reticulate::import("scipy")
dm <- reticulate::import('dbmap')

# Deal with the matrix
data <- t(as.matrix(dat@assays$SCT@data[VariableFeatures(dat),]))
a <- r_to_py(data)
b <- sp$sparse$csr_matrix(a)

# Run the diffusion algorithm
diff = dm$diffusion$Diffusor(n_components = as.integer(200),
                             n_neighbors = as.integer(30),
                             ann_dist = as.character('cosine'),
                             n_jobs = as.integer(10),
                             kernel_use = 'simple',
                             transitions = 'True',
                             norm = 'False')
diff = diff$fit(b)
mms = diff$transform(b)
res = diff$return_dict()

# The diffusion components are the eigencomponents of the structure-weighted diffusion procedure
sc = py_to_r(res$EigenVectors)
ev = py_to_r(res$EigenValues)

################
# This is to deal with the new names
rownames(sc) <- colnames(dat)
new_names <- list()
for(i in 1:length(sc)){
  new_names[i] <- paste('DC_' , as.integer(colnames(sc[i])) + 1, sep = '')
}
colnames(sc) <- as.vector(new_names)
names(ev) <- as.vector(new_names)
################

# Plot the dimensionality of the dataset (similar to PCA's elbow plot)
plot(ev)

# Going back to the Seurat object
dat[["db"]] <- CreateDimReducObject(embeddings = as.matrix(sc), key = "DC_", assay = DefaultAssay(dat))

# Run UMAP for layout
um <- reticulate::import('umap')
umapper <- um$UMAP(min_dist = as.numeric(0.4), spread = as.numeric(1.5), n_epochs = as.integer(1200))
layout <- umapper$fit_transform(sc)
rownames(layout) <- colnames(dat)
plot(layout)

# Going back to the Seurat object
dat[["dbmap"]] <- CreateDimReducObject(embeddings = as.matrix(layout), key = "dbMAP_", assay = DefaultAssay(dat))
DimPlot(dat, reduction = 'dbmap', group.by = 'All.Cell.Clusters', pt.size = 0.5)

DotPlot(dat, features = c('Vim', 'Gfap', 'Cx3cr1', 'Olig2', 'Rbfox3', 'Meg3', 'Nhlh2'), group.by = 'All.Cell.Clusters')

###############################################################################
# Subsetting neurons and producing the plots
###############################################################################

Idents(dat) <- 'All.Cell.Clusters'
neu <- subset(dat, idents = c('a13.Neurons1', 'a14.Neurons2', 'a15.Neurons3', 'a16.Neurons4', 'a17.Neurons5', 'a18.Neurons6'))
DimPlot(neu, reduction = 'dbmap')

umapper1 <- um$UMAP(min_dist = as.numeric(0.8), spread = as.numeric(1.3), n_epochs = as.integer(1000))
layout1 <- umapper1$fit_transform(neu@reductions$db@cell.embeddings)
rownames(layout1) <- colnames(neu)
plot(layout1)
# Going back to the Seurat object
neu[["dbmap"]] <- CreateDimReducObject(embeddings = as.matrix(layout1), key = "dbMAP_", assay = DefaultAssay(neu))
DimPlot(neu, reduction = 'dbmap', group.by = 'Neurons.Only.Clusters', pt.size = 0.5)

FeaturePlot(neu, reduction = 'dbmap', features = c('Pomc', 'Nhlh2'),min.cutoff = c(0,0), max.cutoff = c(3,1), blend = T, pt.size = 1, order = T)
FeaturePlot(neu, reduction = 'dbmap', features = c('Kiss1', 'Nhlh2'),min.cutoff = c(0,0), max.cutoff = c(1,1), blend = T, pt.size = 1, order = T)

FeaturePlot(dat, reduction = 'dbmap', features = c('Cartpt', 'Tac2', 'Agrp', 'Ghrh', 'Sst', 'Trh', 'Th', 'Rgs16'), pt.size = 1, order = T)

DotPlot(neu, features = c('Nhlh2', 'Insr', 'Lepr', 'Glp1r', 'Pcsk1'), group.by = 'Neurons.Only.Clusters')

# Write neurons expression values for SCENIC analysis of transcriptional networks

write.table(as.matrix(neu@assays$RNA@counts), file = 'CampNeurons.tsv', row.names = T, col.names = T)


# Also save as Loom for downstream SCENIC analysis
neu_loom <- create('CampNeurons.loom', data=neu@assays$RNA@counts)



