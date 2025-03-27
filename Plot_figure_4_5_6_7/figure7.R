setwd("/Users/atiteyk2/Documents/CROISSANCE")
getwd()

#####################################################################################################################
##################################################################################################################### Figure 5bcde
#####################################################################################################################
library(monocle3)
library(ggplot2)
library(dplyr)
library(Matrix)
library(qlcMatrix)
library(monocle)
library(Biobase)
library(roahd)
library(uwot)
library(scater)
library(umap)

#####-------------------------- run  monocle --------------------------#####
#https://rpubs.com/WXM/684313
#http://cole-trapnell-lab.github.io/monocle-release/docs/

expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_expression.rds"))
cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_colData.rds"))
gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_rowData.rds"))
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

# How do we find the genes that are differentially expressed on the different paths through the trajectory?
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
driver.genes <- ciliated_cds_pr_test_res$gene_short_name

plot_cells(cds, genes=c("hlh-4", "gcy-8", "dac-1", "oig-8"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)



######################
###################### patient 1
###################### 
load("data_GSE_131907/LUNG.normal.rdata")
load("data_GSE_131907/LUNG.tumor.rdata")
load("data_GSE_131907/LUNG.metastatic.rdata")

LUNG.normal.mat <- LUNG.normal[,2:length(LUNG.normal)]
row.names(LUNG.normal.mat) <- LUNG.normal$data.id

LUNG.tumor.mat <- LUNG.tumor[,2:length(LUNG.tumor)]
row.names(LUNG.tumor.mat) <- LUNG.tumor$data.id

LUNG.metastatic.mat <- LUNG.metastatic[,2:length(LUNG.metastatic)]
row.names(LUNG.metastatic.mat) <- LUNG.metastatic$data.id

LUNG_patient_1 <- cbind(LUNG.normal.mat, LUNG.tumor.mat, LUNG.metastatic.mat)
LUNG_patient_1 <- data.frame(LUNG_patient_1)

ini.lung.patient_1 <- CreateSeuratObject(LUNG_patient_1)

cell.metadata <- ini.lung.patient_1@meta.data
expression.matrix <- ini.lung.patient_1@assays$RNA@layers$counts
annotation <- ini.lung.patient_1@assays$RNA@features
annotation <- data.frame(annotation)

gene_short_name <- LUNG.normal$data.id
gene_short_name <- data.frame(gene_short_name)
gene.annotation <- cbind(gene_short_name, annotation)
rownames(gene.annotation[1,]) <- gene_short_name

cds.patient_1 <- new_cell_data_set(expression.matrix,
                                   cell_metadata = cell.metadata ,
                                   gene_metadata = gene.annotation)
cds.patient_1 <- preprocess_cds(cds.patient_1, num_dim = 50)
B <- reduce_dimension(cds.patient_1, reduction_method = "UMAP", preprocess_method = "PCA")
plot_cells(B)
colnames(colData(B))
plot_cells(B, color_cells_by="Size_Factor")
C  <- cluster_cells(B)
D <- learn_graph(C)
plot_cells(D)
plot_cells(D,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
F <- order_cells(D)
E  <- graph_test(D, neighbor_graph="principal_graph", cores=4)


gene.id <- rownames(E)
driver.data.coef <- E$morans_I
driver.data.coef <- abs(driver.data.coef)
nan_indices <- is.nan(driver.data.coef)
driver.data.coef[nan_indices] <- 0
driver.data.coef <- cbind(gene.id, driver.data.coef)
driver.data.coef <- data.frame(driver.data.coef)
colnames(driver.data.coef) <- c("gene.id", "morans")
driver.data.coef$morans <- as.numeric(driver.data.coef$morans)
moran.gene.driver <- driver.data.coef[order(-driver.data.coef$morans),]
save(moran.gene.driver, file = "moran.gene.driver.rdata")
head(moran.gene.driver, 10)


######################
###################### patient 2
###################### 
load("data_GSE_19/LUNG.normal19.data.rdata")
load("LUNG.tumor19.data.rdata")
load("data_GSE_19/LUNG.metas19.data.rdata")

LUNG.normal.mat <- LUNG.normal19.data[,2:length(LUNG.normal19.data)]
row.names(LUNG.normal.mat) <- LUNG.normal19.data$data.id

LUNG.tumor.mat <- LUNG.tumor19.data[,2:length(LUNG.tumor19.data)]
row.names(LUNG.tumor.mat) <- LUNG.tumor19.data$data.id

LUNG.metastatic.mat <- LUNG.metas19.data[,2:length(LUNG.metas19.data)]
row.names(LUNG.metastatic.mat) <- LUNG.metas19.data$data.id

LUNG_patient_2 <- cbind(LUNG.normal.mat, LUNG.tumor.mat, LUNG.metastatic.mat)
LUNG_patient_2 <- data.frame(LUNG_patient_2)

ini.lung.patient_2 <- CreateSeuratObject(LUNG_patient_2)

cell.metadata <- ini.lung.patient_2@meta.data
expression.matrix <- ini.lung.patient_2@assays$RNA@layers$counts
annotation <- ini.lung.patient_2@assays$RNA@features
annotation <- data.frame(annotation)

gene_short_name <- LUNG.tumor19.data$data.id
gene_short_name <- data.frame(gene_short_name)
gene.annotation <- cbind(gene_short_name, annotation)
rownames(gene.annotation) <- gene_short_name

cds.patient_2 <- new_cell_data_set(expression.matrix,
                                   cell_metadata = cell.metadata ,
                                   gene_metadata = gene.annotation)
cds.patient_2 <- preprocess_cds(cds.patient_2, num_dim = 50)

B <- reduce_dimension(cds.patient_2, reduction_method = "UMAP", preprocess_method = "PCA")
plot_cells(B)
colnames(colData(B))
plot_cells(B, color_cells_by="Size_Factor")
C  <- cluster_cells(B)
D <- learn_graph(C)
plot_cells(D)
plot_cells(D,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
F <- order_cells(D)
E  <- graph_test(D, neighbor_graph="principal_graph", cores=4)

gene.id <- rownames(E)
driver.data.coef <- E$morans_I
driver.data.coef <- abs(driver.data.coef)
nan_indices <- is.nan(driver.data.coef)
driver.data.coef[nan_indices] <- 0
driver.data.coef <- cbind(gene.id, driver.data.coef)
driver.data.coef <- data.frame(driver.data.coef)
colnames(driver.data.coef) <- c("gene.id", "morans")
driver.data.coef$morans <- as.numeric(driver.data.coef$morans)
moran.gene.driver <- driver.data.coef[order(-driver.data.coef$morans),]
moran.gene.order <- moran.gene.driver

######################
###################### patient 3
###################### 
load("data_GSE_131907_8/LUNG.normal.8.rdata")
load("data_GSE_131907_8/LUNG.tumor.8.rdata")
load("data_GSE_131907_8/LUNG.metas.8.rdata")

LUNG.normal.mat <- LUNG.normal.8[,2:length(LUNG.normal.8)]
row.names(LUNG.normal.mat) <- LUNG.normal.8$data.id

LUNG.tumor.mat <- LUNG.tumor.8[,2:length(LUNG.tumor.8)]
row.names(LUNG.tumor.mat) <- LUNG.tumor.8$data.id

LUNG.metastatic.mat <- LUNG.metas.8[,2:length(LUNG.metas.8)]
row.names(LUNG.metastatic.mat) <- LUNG.metas.8$data.id

LUNG_patient_3 <- cbind(LUNG.normal.mat, LUNG.tumor.mat, LUNG.metastatic.mat)
LUNG_patient_3 <- data.frame(LUNG_patient_3)

ini.lung.patient_3 <- CreateSeuratObject(LUNG_patient_3)

cell.metadata <- ini.lung.patient_3@meta.data
expression.matrix <- ini.lung.patient_3@assays$RNA@layers$counts
annotation <- ini.lung.patient_3@assays$RNA@features
annotation <- data.frame(annotation)

gene_short_name <- LUNG.normal.8$data.id
gene_short_name <- data.frame(gene_short_name)
gene.annotation <- cbind(gene_short_name, annotation)
rownames(gene.annotation) <- gene_short_name

cds.patient_3 <- new_cell_data_set(expression.matrix,
                                   cell_metadata = cell.metadata ,
                                   gene_metadata = gene.annotation)
cds.patient_3 <- preprocess_cds(cds.patient_3, num_dim = 50)

B <- reduce_dimension(cds.patient_3, reduction_method = "UMAP", preprocess_method = "PCA")
plot_cells(B)
colnames(colData(B))
plot_cells(B, color_cells_by="Size_Factor")
C  <- cluster_cells(B)
D <- learn_graph(C)
plot_cells(D)
plot_cells(D,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
F <- order_cells(D)
E  <- graph_test(D, neighbor_graph="principal_graph", cores=4)

gene.id <- rownames(E)
driver.data.coef <- E$morans_I
driver.data.coef <- abs(driver.data.coef)
nan_indices <- is.nan(driver.data.coef)
driver.data.coef[nan_indices] <- 0
driver.data.coef <- cbind(gene.id, driver.data.coef)
driver.data.coef <- data.frame(driver.data.coef)
colnames(driver.data.coef) <- c("gene.id", "morans")
driver.data.coef$morans <- as.numeric(driver.data.coef$morans)
moran.gene.driver <- driver.data.coef[order(-driver.data.coef$morans),]
moran.gene.order <- moran.gene.driver

####################
LUNG.meta.mat <- LUNG.metas.8[,2:length(LUNG.metas.8)]
row.names(LUNG.meta.mat) <- LUNG.metas.8$data.id
ini.lung.meta <- CreateSeuratObject(LUNG.meta.mat)

#ge_annotation <- ini.lung.meta@assays$RNA@features # featureData
cell.metadata <- ini.lung.meta@meta.data
expression.matrix <- ini.lung.meta@assays$RNA@layers$counts

annotation <- ini.lung.meta@assays$RNA@features
annotation <- data.frame(annotation)
gene_short_name <- LUNG.metas.8$data.id
gene_short_name <- data.frame(gene_short_name)
gene.annotation <- cbind(gene_short_name, annotation)
rownames(gene.annotation) <- gene_short_name

cds.meta <- new_cell_data_set(expression.matrix,
                              cell_metadata = cell.metadata ,
                              gene_metadata = gene.annotation)
cds.meta <- preprocess_cds(cds.meta, num_dim = 50)
cds.meta <- reduce_dimension(cds.meta, reduction_method="UMAP")
cds.meta <- cluster_cells(cds.meta)
cds.meta <- learn_graph(cds.meta)
#cds.meta <- order_cells(cds.meta)
ciliated_cds.meta_pr_test_res <- graph_test(cds.meta, neighbor_graph="principal_graph", cores=4)

###
gene.id <- rownames(ciliated_cds.meta_pr_test_res)
driver.data.coef <- ciliated_cds.meta_pr_test_res$morans_I
driver.data.coef <- abs(driver.data.coef)
nan_indices <- is.nan(driver.data.coef)
driver.data.coef[nan_indices] <- 0
driver.data.coef <- cbind(gene.id, driver.data.coef)
driver.data.coef <- data.frame(driver.data.coef)
colnames(driver.data.coef) <- c("gene.id", "morans")
driver.data.coef$morans <- as.numeric(driver.data.coef$morans)
moran.gene.driver <- driver.data.coef[order(-driver.data.coef$morans),]
moran.gene.order <- moran.gene.driver
##############
load("moran.gene.order.data.rdata")
load("picdgi.gene.order.data.rdata")

rownames(picdgi.gene.order.data) <- NULL
rownames(moran.gene.order.data) <- NULL

moran.gene.order <- data.frame(moran.gene.order.data)
picdgi.gene.order <- data.frame(picdgi.gene.order.data)

NK <- cbind(moran.gene.order$NK, picdgi.gene.order$NK)
Tcell <- cbind(moran.gene.order$T.lymphocytes, picdgi.gene.order$T.lymphocytes)
Bcell <- cbind(moran.gene.order$B.lymphocytes, picdgi.gene.order$B.lymphocytes)
DC <- cbind(moran.gene.order$DC, picdgi.gene.order$DC)
Mast <- cbind(moran.gene.order$Mast, picdgi.gene.order$Mast)
NK <- NK[1:10,]

#####-------------------------- Figure 7B,C,D,E,F --------------------------#####
library(tidyr)
library(ggplot2)
##### NK
X <- NK[,1]
Y <- NK[,2]
position <- 1:10
df.NK <- data.frame(X, Y, position)
colour_palette <- c("darkred", "forestgreen")
df.NK <- gather(df.NK, event, total, X:Y)
plot <- ggplot(df.NK, aes(position, total, fill=event)) 
plot <- plot + geom_bar(stat = "identity", position = 'dodge') + theme_classic() + scale_fill_manual(values = colour_palette)
plot 

##### Tcell
X <- Tcell[,1]
Y <- Tcell[,2]
position <- 1:10
df.NK <- data.frame(X, Y, position)
colour_palette <- c("darkred", "forestgreen")
df.NK <- gather(df.NK, event, total, X:Y)
plot <- ggplot(df.NK, aes(position, total, fill=event)) 
plot <- plot + geom_bar(stat = "identity", position = 'dodge') + theme_classic() + scale_fill_manual(values = colour_palette)
plot 

##### Bcell
Y <- Bcell[,1]
X <- Bcell[,2]
position <- 1:10
df.NK <- data.frame(X, Y, position)
colour_palette <- c("darkred", "forestgreen")
df.NK <- gather(df.NK, event, total, X:Y)
plot <- ggplot(df.NK, aes(position, total, fill=event)) 
plot <- plot + geom_bar(stat = "identity", position = 'dodge') + theme_classic() + scale_fill_manual(values = colour_palette)
plot 

##### DC
X <- DC[,1]
Y <- DC[,2]
position <- 1:10
df.NK <- data.frame(X, Y, position)
colour_palette <- c("darkred", "forestgreen")
df.NK <- gather(df.NK, event, total, X:Y)
plot <- ggplot(df.NK, aes(position, total, fill=event)) 
plot <- plot + geom_bar(stat = "identity", position = 'dodge') + theme_classic() + scale_fill_manual(values = colour_palette)
plot 

##### Mast
Y <- Mast[,1]
X <- Mast[,2]
position <- 1:10
df.NK <- data.frame(X, Y, position)
colour_palette <- c("darkred", "forestgreen")
df.NK <- gather(df.NK, event, total, X:Y)
plot <- ggplot(df.NK, aes(position, total, fill=event)) 
plot <- plot + geom_bar(stat = "identity", position = 'dodge') + theme_classic() + scale_fill_manual(values = colour_palette)
plot 

