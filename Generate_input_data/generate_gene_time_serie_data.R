setwd("/Users/atiteyk2/Documents/CROISSANCE")
getwd()

#####-------------------------- time 1 --------------------------#####

#######
####### order cell types in terms of their id
#######

tsne.data.normal <- tsne.lung.normal@reductions$tsne@cell.embeddings #ini.lung.normal@reductions$tsne@cell.embeddings
cell.types.normal <- tsne.lung.normal@active.ident #ini.lung.normal@active.ident
save(tsne.data.normal, file = "tsne.data.normal.rdata")
save(cell.types.normal, file = "cell.types.normal.rdata")

normal.gene.cell <- tsne.lung.normal@assays$RNA@scale.data
normal.gene.cell <- data.frame(normal.gene.cell)
gene.conc <- VlnPlot(tsne.lung.normal, features="IL7")
cells.id.normal <- gene.conc$data$ident
save(cells.id.normal, file="cells.id.normal.rdata")

cell.id.char.normal <- as.character(cells.id.normal)
cell.id.char.normal <- as.list(cell.id.char.normal)
colnames(normal.gene.cell) <- c(cell.id.char.normal)
save(normal.gene.cell, file = "normal.gene.cell.rdata")

#######
####### compute gene expression level in a cell at time 1
#######

load("lung_cancer/gene_cell_data/normal.gene.cell.rdata")
normal.gene.cell <- data.frame(normal.gene.cell)
epithelial.data.time1 <- normal.gene.cell %>% select(starts_with("Epithelial_cells"))
#
row.names(epithelial.data.time1) <- 1:nrow(epithelial.data.time1)
epithelial.data.time1 <- data.frame(epithelial.data.time1)

N <- nrow(epithelial.data.time1) # Number of genes
M <- ncol(epithelial.data.time1) # Number of cells

epithelial.mat.time1 <- sapply(1:N, function(i) {matrix(epithelial.data.time1[i,1:M])}) # function to assigne data to each gene
epithelial.mat.time1 <- data.frame(epithelial.mat.time1) # convert to data frame
epithelial.mat.time1 <- as.matrix(epithelial.mat.time1) # convert data to a matrix

epithelial.time1 <- matrix(unlist(epithelial.mat.time1), ncol = dim(data.frame(epithelial.mat.time1))[2],
                              byrow = TRUE)

gene.id <- rownames(normal.gene.cell) # identity of genes
epithelial.time1 <- exp(epithelial.time1)
colnames(epithelial.time1) <- c(gene.id) # assign gene identity to the columns names

#######
####### compute the IQR of gene expression level in a cell at time 1
#######

expres.level.time1 <- lapply(1:ncol(epithelial.time1), function(j) {IQR(epithelial.time1[,j])})    # compute the interqutile range 
expres.level.time1 <- data.frame(expres.level.time1)
colnames(expres.level.time1) <- c(gene.id)
expres.level.time1 <- t(expres.level.time1)

id.time1 <- rownames(expres.level.time1)
row.names(expres.level.time1) <- 1:nrow(data.frame(expres.level.time1))
epithelial.level.time1 <- cbind(id.time1, expres.level.time1)

epithelial.level.time1 <- data.frame(epithelial.level.time1)
colnames(epithelial.level.time1)[colnames(epithelial.level.time1) == "V2"] <- "level_1"         # change column name
colnames(epithelial.level.time1)[colnames(epithelial.level.time1) == "id.time1"] <- "gene"      # change column name

#####-------------------------- time 2 --------------------------#####

#######
####### order cell types in terms of their id
#######
tsne.data.tumor <- ini.lung.tumor@reductions$tsne@cell.embeddings
cell.types.tumor <- ini.lung.tumor@active.ident
save(tsne.data.tumor, file = "tsne.data.tumor.rdata")
save(cell.types.tumor, file = "cell.types.tumor.rdata")

tumor.gene.cell <- ini.lung.tumor@assays$RNA@scale.data
tumor.gene.cell <- data.frame(tumor.gene.cell)
gene.conc <- VlnPlot(ini.lung.tumor, features="IL7")
cells.id.tumor <- gene.conc$data$ident
save(cells.id.tumor, file="cells.id.tumor.rdata")

cell.id.char.tumor <- as.character(cells.id.tumor)
cell.id.char.tumor <- as.list(cell.id.char.tumor)
colnames(tumor.gene.cell) <- c(cell.id.char.tumor)
save(tumor.gene.cell, file = "tumor.gene.cell.rdata")

#######
####### compute gene expression level in a cell at time 2
#######

load("lung_cancer/gene_cell_data/tumor.gene.cell.rdata")
tumor.gene.cell <- data.frame(tumor.gene.cell)
epithelial.data.time2 <- tumor.gene.cell %>% select(starts_with("Epithelial_cells"))
#
row.names(epithelial.data.time2) <- 1:nrow(epithelial.data.time2)
epithelial.data.time2 <- data.frame(epithelial.data.time2)

N <- nrow(epithelial.data.time2) # Number of genes
M <- ncol(epithelial.data.time2) # Number of cells

epithelial.mat.time2 <- sapply(1:N, function(i) {matrix(epithelial.data.time2[i,1:M])}) # function to assigne data to each gene
epithelial.mat.time2 <- data.frame(epithelial.mat.time2) # convert to data frame
epithelial.mat.time2 <- as.matrix(epithelial.mat.time2) # convert data to a matrix

epithelial.time2 <- matrix(unlist(epithelial.mat.time2), ncol = dim(data.frame(epithelial.mat.time2))[2],
                              byrow = TRUE)

gene.id <- rownames(tumor.gene.cell) # identity of genes
epithelial.time2 <- exp(epithelial.time2)
colnames(epithelial.time2) <- c(gene.id) # assign gene identity to the columns names

#######
####### compute the IQR of gene expression level in a cell at time 2
#######

expres.level.time2 <- lapply(1:ncol(epithelial.time2), function(j) {IQR(epithelial.time2[,j])})      # compute the interqutile range
expres.level.time2 <- data.frame(expres.level.time2)
colnames(expres.level.time2) <- c(gene.id)
expres.level.time2 <- t(expres.level.time2)

id.time2 <- rownames(expres.level.time2)
row.names(expres.level.time2) <- 1:nrow(data.frame(expres.level.time2))
epithelial.level.time2 <- cbind(id.time2, expres.level.time2)

epithelial.level.time2 <- data.frame(epithelial.level.time2)
colnames(epithelial.level.time2)[colnames(epithelial.level.time2) == "V2"] <- "level_2"          # change column name
colnames(epithelial.level.time2)[colnames(epithelial.level.time2) == "id.time2"] <- "gene"       # change column name

#####-------------------------- time 3 --------------------------#####

#######
####### order cell types in terms of their id
#######

tsne.data.meta <- ini.lung.meta@reductions$tsne@cell.embeddings
cell.types.meta <- ini.lung.meta@active.ident
save(tsne.data.meta, file = "tsne.data.meta.rdata")
save(cell.types.meta, file = "cell.types.meta.rdata")

meta.gene.cell <- ini.lung.meta@assays$RNA@scale.data
meta.gene.cell <- data.frame(meta.gene.cell)
gene.conc <- VlnPlot(ini.lung.meta, features="IL7")
cells.id.meta <- gene.conc$data$ident
save(cells.id.meta, file="cells.id.meta.rdata")

cell.id.char.meta <- as.character(cells.id.meta)
cell.id.char.meta <- as.list(cell.id.char.meta)
colnames(meta.gene.cell) <- c(cell.id.char.meta)
save(meta.gene.cell, file = "meta.gene.cell.rdata")

#######
####### compute gene expression level in a cell at time 3
#######

load("lung_cancer/gene_cell_data/meta.gene.cell.rdata")
meta.gene.cell <- data.frame(meta.gene.cell)
epithelial.data.time3 <- meta.gene.cell %>% select(starts_with("Epithelial_cells"))

row.names(epithelial.data.time3) <- 1:nrow(epithelial.data.time3)
epithelial.data.time3 <- data.frame(epithelial.data.time3)

N <- nrow(epithelial.data.time3) # Number of genes
M <- ncol(epithelial.data.time3) # Number of cells

epithelial.mat.time3 <- sapply(1:N, function(i) {matrix(epithelial.data.time3[i,1:M])}) # function to assigne data to each gene
epithelial.mat.time3 <- data.frame(epithelial.mat.time3) # convert to data frame
epithelial.mat.time3 <- as.matrix(epithelial.mat.time3) # convert data to a matrix

epithelial.time3 <- matrix(unlist(epithelial.mat.time3), ncol = dim(data.frame(epithelial.mat.time3))[2],
                              byrow = TRUE)

gene.id <- rownames(meta.gene.cell) # identity of genes
epithelial.time3 <- exp(epithelial.time3)
colnames(epithelial.time3) <- c(gene.id) # assign gene identity to the columns names

#######
####### compute the IQR of gene expression level in a cell at time 3
#######

expres.level.time3 <- lapply(1:ncol(epithelial.time3), function(j) {IQR(epithelial.time3[,j])})      # compute the interqutile range
expres.level.time3 <- data.frame(expres.level.time3)
colnames(expres.level.time3) <- c(gene.id)
expres.level.time3 <- t(expres.level.time3)

id.time3 <- rownames(expres.level.time3)
row.names(expres.level.time3) <- 1:nrow(data.frame(expres.level.time3))
epithelial.level.time3 <- cbind(id.time3, expres.level.time3)

epithelial.level.time3 <- data.frame(epithelial.level.time3)
colnames(epithelial.level.time3)[colnames(epithelial.level.time3) == "V2"] <- "level_3"         # change column name
colnames(epithelial.level.time3)[colnames(epithelial.level.time3) == "id.time3"] <- "gene"      # change column name




