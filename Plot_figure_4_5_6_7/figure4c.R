setwd("/Users/atiteyk2/Documents/CROISSANCE")
getwd()

library(dplyr)
library(Seurat)
library(patchwork)
#rm(list=ls())
library(knitr)
library(ggplot2)
library(tidyverse)
library(tidyr)

########## CROISSANCE
load("tsne_reduced_data/tsne.lung.normal.rdata")
load("tsne_reduced_data/tsne.lung.tumor.rdata")
load("tsne_reduced_data/tsne.lung.meta.rdata")

########## EVOLUTION
load("tsne_data/tsne.lung.normal.rdata")
load("tsne_data/tsne.lung.tumor.rdata")
load("tsne_data/tsne.lung.meta.rdata")

########## EVOLUTION2
load("tsne_data/tsne.lung.normal.rdata")
load("tsne_data/tsne.lung.tumor.rdata")
load("tsne_data/tsne.lung.meta.rdata")

##########################################################################
########################################################################## clusters ids 
##########################################################################
#####-------------------------- patient 1
# normal
new.cluster.ids <- c("T_lymphocytes", "Dendritic_cells", "NK_cells", "Epithelial_cells", "Fibroblasts", 
                     "Mast_cells", "Endothelial_cells", "B_lymphocytes", "Ependymal_cells", "Endothelial_cells")
# tumor
new.cluster.ids <- c("NK_cells", "T_lymphocytes", "Dendritic_cells", "Endothelial_cells", "Mast_cells",  "Epithelial_cells", "Fibroblasts")
# metastasis
new.cluster.ids <- c("NK_cells", "T_lymphocytes", "Dendritic_cells", "Epithelial_cells", "Fibroblasts",
                     "Mast_cells", "B_lymphocytes", "Ependymal_cells") 

#####-------------------------- patient 2
# normal
new.cluster.ids <- c("T_lymphocytes", "Dendritic_cells", "B_lymphocytes", "Epithelial_cells", "Mast_cells",
                     "NK_cells", "NK_cells", "B_lymphocytes", "Epithelial_cells", "Fibroblasts",
                     "Endothelial_cells", "Epithelial_cells", "Ependymal_cells")
# tumor
new.cluster.ids <- c("T_lymphocytes", "Dendritic_cells", "B_lymphocytes", "Mast_cells",
                     "NK_cells", "Epithelial_cells", "Endothelial_cells", "Epithelial_cells")
# metastasis
new.cluster.ids <- c("T_lymphocytes", "B_lymphocytes", "NK_cells", "Dendritic_cells", "Epithelial_cells", "Epithelial_cells",
                     "Fibroblasts", "Mast_cells" )

#####-------------------------- patient 3
#normal
new.cluster.ids <- c("T_lymphocytes", "Dendritic_cells" , "Epithelial_cells", "Fibroblasts",
                     "Epithelial_cells", "Epithelial_cells", "B_lymphocytes", "B_lymphocytes",
                     "Oligodendrocytes", "Mast_cells", "Epithelial_cells", "Fibroblasts")
# tumor
new.cluster.ids <- c( "Dendritic_cells" , "Epithelial_cells", "T_lymphocytes", "B_lymphocytes", "Epithelial_cells",
                      "Epithelial_cells")
# metastasis
new.cluster.ids <- c( "Epithelial_cells", "Epithelial_cells", "Dendritic_cells" , "T_lymphocytes", "Epithelial_cells",
                      "Epithelial_cells","Fibroblasts", "B_lymphocytes", "Oligodendrocytes")

##########################################################################
########################################################################## select clusters data
##########################################################################
#####-------------------------- normal
names(new.cluster.ids) <- levels(tsne.lung.normal)
normal.gene.cell <- RenameIdents(tsne.lung.normal, new.cluster.ids)
normal.gene.cell <- tsne.lung.normal@assays$RNA$counts
normal.gene.cell <- data.frame(normal.gene.cell)
cell.types <- tsne.lung.normal@active.ident
levels(cell.types) <- new.cluster.ids
C.types <- data.frame(cell.types)
colnames(normal.gene.cell) <- C.types$cell.types

normal.gene.cell <- data.frame(normal.gene.cell)
T_lymphocytes.time1 <- normal.gene.cell %>% select(starts_with("T_lymphocytes"))
DC.time1 <- normal.gene.cell %>% select(starts_with("Dendritic_cells"))
NK.time1 <- normal.gene.cell %>% select(starts_with("NK_cells"))
Epithelial.time1 <- normal.gene.cell %>% select(starts_with("Epithelial_cells"))
Fibroblasts.time1 <- normal.gene.cell %>% select(starts_with("Fibroblasts"))
Mast.time1 <- normal.gene.cell %>% select(starts_with("Mast_cells"))
Endothelial.time1 <- normal.gene.cell %>% select(starts_with("Endothelial_cells"))
B_lymphocytes.time1 <- normal.gene.cell %>% select(starts_with("B_lymphocyte"))
Ependymal.time1 <- normal.gene.cell %>% select(starts_with("Ependymal_cells"))

#####-------------------------- tumor
tumor.gene.cell <- tsne.lung.tumor@assays$RNA$counts
tumor.gene.cell <- data.frame(tumor.gene.cell)
cell.types <- tsne.lung.tumor@active.ident
C.types <- data.frame(cell.types)
colnames(tumor.gene.cell) <- C.types$cell.types
tumor.gene.cell <- data.frame(tumor.gene.cell)

T_lymphocytes.time2 <- tumor.gene.cell %>% select(starts_with("T_lymphocytes"))
DC.time2 <- tumor.gene.cell %>% select(starts_with("Dendritic_cells"))
B_lymphocytes.time2 <- tumor.gene.cell %>% select(starts_with("B_lymphocyte"))
Epithelial.time2 <- tumor.gene.cell %>% select(starts_with("Epithelial_cells"))
Mast.time2 <- tumor.gene.cell %>% select(starts_with("Mast_cells"))
NK.time2 <- tumor.gene.cell %>% select(starts_with("NK_cells"))
Fibroblasts.time2 <- tumor.gene.cell %>% select(starts_with("Fibroblasts"))
Endothelial.time2 <- tumor.gene.cell %>% select(starts_with("Endothelial_cells"))
Ependymal.time2 <- tumor.gene.cell %>% select(starts_with("Ependymal_cells"))
Oligodendrocytes.time2 <- tumor.gene.cell %>% select(starts_with("Oligodendrocytes"))

#####-------------------------- metastasis
meta.gene.cell <- tsne.lung.meta@assays$RNA$counts
meta.gene.cell <- data.frame(meta.gene.cell)
cell.types <- tsne.lung.meta@active.ident
C.types <- data.frame(cell.types)
colnames(meta.gene.cell) <- C.types$cell.types
meta.gene.cell <- data.frame(meta.gene.cell)

T_lymphocytes.time3 <- meta.gene.cell %>% select(starts_with("T_lymphocytes"))
DC.time3 <- meta.gene.cell %>% select(starts_with("Dendritic_cells"))
Epithelial.time3 <- meta.gene.cell %>% select(starts_with("Epithelial_cells"))
Fibroblasts.time3 <- meta.gene.cell %>% select(starts_with("Fibroblasts"))
B_lymphocytes.time3 <- meta.gene.cell %>% select(starts_with("B_lymphocyte"))
Mast.time3 <- meta.gene.cell %>% select(starts_with("Mast_cells"))

##########################################################################
########################################################################## compute cell normal and tumor fractions
##########################################################################

###### T lymphocyte cell
non_maligant.T <- subset(T_lymphocytes.time1, row.names(T_lymphocytes.time1) %in% "CD3D") # select normal cell
non_maligant.T <- as.numeric(non_maligant.T)
non_maligant.T <- length(non_maligant.T) -  sum(non_maligant.T == 0)
non_maligant.T <- non_maligant.T / length(T_lymphocytes.time1)

tumor.T <- subset(T_lymphocytes.time1, row.names(T_lymphocytes.time1) %in% "EPCAM") # select cancer cell
tumor.T <- as.numeric(tumor.T)
tumor.T <- length(tumor.T) -  sum(tumor.T == 0)
tumor.T <- tumor.T / length(T_lymphocytes.time1)

###### DC cell
non_maligant.DC <- subset(DC.time1, row.names(DC.time1) %in% "MS4A7") # select normal cell
non_maligant.DC <- as.numeric(non_maligant.DC)
non_maligant.DC <- length(non_maligant.DC) -  sum(non_maligant.DC == 0)
non_maligant.DC <- non_maligant.DC / length(DC.time1)

tumor.DC <- subset(DC.time1, row.names(DC.time1) %in% "EPCAM") # select cancer cell
tumor.DC <- as.numeric(tumor.DC)
tumor.DC <- length(tumor.DC) -  sum(tumor.DC == 0)
tumor.DC <- tumor.DC / length(DC.time1)

###### NK cell
non_maligant.NK <- subset(NK.time1, row.names(NK.time1) %in% "KLRF1") # select normal cell
non_maligant.NK <- as.numeric(non_maligant.NK)
non_maligant.NK <- length(non_maligant.NK) -  sum(non_maligant.NK == 0)
non_maligant.NK <- non_maligant.NK / length(NK.time1)

tumor.NK <- subset(NK.time1, row.names(NK.time1) %in% "EPCAM") # select cancer cell
tumor.NK <- as.numeric(tumor.NK)
tumor.NK <- length(tumor.NK) -  sum(tumor.NK == 0)
tumor.NK <- tumor.NK / length(NK.time1)

###### Epithelial cell
non_maligant.epi <- subset(Epithelial.time1, row.names(Epithelial.time1) %in% "SFTA2") # select normal cell
non_maligant.epi <- as.numeric(non_maligant.epi)
non_maligant.epi <- length(non_maligant.epi) -  sum(non_maligant.epi == 0)
non_maligant.epi <- non_maligant.epi / length(Epithelial.time1)

tumor.epi <- subset(Epithelial.time1, row.names(Epithelial.time1) %in% "EPCAM") # select cancer cell
tumor.epi <- as.numeric(tumor.epi)
tumor.epi <- length(tumor.epi) -  sum(tumor.epi == 0)
tumor.epi <- tumor.epi / length(Epithelial.time1)

###### Fibroblasts cell
non_maligant.fi <- subset(Fibroblasts.time1, row.names(Fibroblasts.time1) %in% "COL1A2") # select normal cell
non_maligant.fi <- as.numeric(non_maligant.fi)
non_maligant.fi <- length(non_maligant.fi) -  sum(non_maligant.fi == 0)
non_maligant.fi <- non_maligant.fi / length(Fibroblasts.time1)

tumor.fi <- subset(Fibroblasts.time1, row.names(Fibroblasts.time1) %in% "EPCAM") # select cancer cell
tumor.fi <- as.numeric(tumor.fi)
tumor.fi <- length(tumor.fi) -  sum(tumor.fi == 0)
tumor.fi <- tumor.fi / length(Fibroblasts.time1)

###### Mast cell
non_maligant.mast <- subset(Mast.time1, row.names(Mast.time1) %in% "MS4A2") # select normal cell
non_maligant.mast <- as.numeric(non_maligant.mast)
non_maligant.mast <- length(non_maligant.mast) -  sum(non_maligant.mast == 0)
non_maligant.mast <- non_maligant.mast / length(Mast.time1)

tumor.mast <- subset(Mast.time1, row.names(Mast.time1) %in% "EPCAM") # select cancer cell
tumor.mast <- as.numeric(tumor.mast)
tumor.mast <- length(tumor.mast) -  sum(tumor.mast == 0)
tumor.mast <- tumor.mast / length(Mast.time1)

###### Endothelial cell
non_maligant.endo <- subset(Endothelial.time1, row.names(Endothelial.time1) %in% "CLDN5") # select normal cell
non_maligant.endo <- as.numeric(non_maligant.endo)
non_maligant.endo <- length(non_maligant.endo) -  sum(non_maligant.endo == 0)
non_maligant.endo <- non_maligant.endo / length(Endothelial.time1)

tumor.endo <- subset(Endothelial.time1, row.names(Endothelial.time1) %in% "EPCAM") # select cancer cell
tumor.endo <- as.numeric(tumor.endo)
tumor.endo <- length(tumor.endo) -  sum(tumor.endo == 0)
tumor.endo <- tumor.endo / length(Endothelial.time1)

###### B_lymphocytes cell
non_maligant.B <- subset(B_lymphocytes.time1, row.names(B_lymphocytes.time1) %in% "CD79A") # select normal cell
non_maligant.B <- as.numeric(non_maligant.B)
non_maligant.B <- length(non_maligant.B) -  sum(non_maligant.B == 0)
non_maligant.B <- non_maligant.B / length(B_lymphocytes.time1)

tumor.B <- subset(B_lymphocytes.time1, row.names(B_lymphocytes.time1) %in% "EPCAM") # select cancer cell
tumor.B <- as.numeric(tumor.B)
tumor.B <- length(tumor.B) -  sum(tumor.B == 0)
tumor.B <- tumor.B / length(B_lymphocytes.time1)

###### Ependymal cell
non_maligant.epen <- subset(Ependymal.time1, row.names(Ependymal.time1) %in% "FAM183A") # select normal cell
non_maligant.epen <- as.numeric(non_maligant.epen)
non_maligant.epen <- length(non_maligant.epen) -  sum(non_maligant.epen == 0)
non_maligant.epen <- non_maligant.epen / length(Ependymal.time1)

tumor.epen <- subset(Ependymal.time1, row.names(Ependymal.time1) %in% "EPCAM") # select cancer cell
tumor.epen <- as.numeric(tumor.epen)
tumor.epen <- length(tumor.epen) -  sum(tumor.epen == 0)
tumor.epen <- tumor.epen / length(Ependymal.time1)

###### Oligodentrocytes
non_maligant.oli <- subset(Oligodendrocytes.time1, row.names(Oligodendrocytes.time1) %in% "FAM183A") # select normal cell
non_maligant.oli <- as.numeric(non_maligant.oli)
non_maligant.oli <- length(non_maligant.oli) -  sum(non_maligant.oli== 0)
non_maligant.oli <- non_maligant.oli / length(Oligodendrocytes.time1)

tumor.oli <- subset(Oligodendrocytes.time1, row.names(Oligodendrocytes.time1) %in% "EPCAM") # select cancer cell
tumor.oli <- as.numeric(tumor.oli)
tumor.oli <- length(tumor.oli) -  sum(tumor.oli == 0)
tumor.oli <- tumor.oli / length(Oligodendrocytes.time1)

##########################################################################
########################################################################## combine the data
##########################################################################
#####-------------------------- normal
non_tumor_1 <- c(non_maligant.T, non_maligant.DC, non_maligant.NK, non_maligant.epi, non_maligant.fi, non_maligant.mast, non_maligant.endo, non_maligant.B, non_maligant.epen)
tumor_1 <- c(tumor.T, tumor.DC, tumor.NK, tumor.epi, tumor.fi, tumor.mast, tumor.endo, tumor.B, tumor.epen)
cells.fraction.time1 <- cbind(non_tumor_1, tumor_1)
cells.fraction.time1 <- data.frame(cells.fraction.time1)
rownames(cells.fraction.time1) <- c("T_lymphocytes", "Dendritic_cells", "NK_cells", "Epithelial_cells", "Fibroblasts", 
                                    "Mast_cells", "Endothelial_cells", "B_lymphocytes", "Ependymal_cells")
save(cells.fraction.time1, file = "cells.fraction.time1.rdata")

#####-------------------------- tumor
non_tumor_2 <- c(non_maligant.T, non_maligant.DC, non_maligant.NK, non_maligant.epi, non_maligant.fi, non_maligant.mast,
               non_maligant.endo, non_maligant.B, non_maligant.epen, non_maligant_oli)
tumor_2 <- c(tumor.T, tumor.DC, tumor.NK, tumor.epi, tumor.fi, tumor.mast, tumor.endo, tumor.B, tumor.epen, tumor_oli)
cells.fraction.time2 <- cbind(non_tumor_2, tumor_2)
cells.fraction.time2 <- data.frame(cells.fraction.time2)
rownames(cells.fraction.time2) <- c("T_lymphocytes", "Dendritic_cells", "NK_cells", "Epithelial_cells", "Fibroblasts", 
                                    "Mast_cells", "Endothelial_cells", "B_lymphocytes", "Ependymal_cells", "Oligodendrocytes")
save(cells.fraction.time2, file = "cells.fraction.time2.rdata")

#####-------------------------- metastasis
non_tumor_3 <- c(non_maligant.T, non_maligant.DC, 0, non_maligant.epi, non_maligant.fi, non_maligant.mast,
               0, non_maligant.B, 0, non_maligant.oli)
tumor_3 <- c(tumor.T, tumor.DC, 0, tumor.epi, tumor.fi, tumor.mast, 0, tumor.B, 0, tumor.oli)
cells.fraction.time3 <- cbind(non_tumor_3, tumor_3)
cells.fraction.time3 <- data.frame(cells.fraction.time3)
rownames(cells.fraction.time3) <- c("T_lymphocytes", "Dendritic_cells", "NK_cells", "Epithelial_cells", "Fibroblasts", 
                                    "Mast_cells", "Endothelial_cells", "B_lymphocytes", "Ependymal_cells", "Oligodendrocytes")
save(cells.fraction.time3, file = "cells.fraction.time3.rdata")

##########################################################################
########################################################################## plot
##########################################################################
load("cells.fraction.time1.rdata")
load("cells.fraction.time2.rdata")
load("cells.fraction.time3.rdata")

#####-------------------------- patient 1

rownames(cells.fraction.time1) <- c("T_lymphocytes.1", "Dendritic_cells.1", "NK_cells.1", "Epithelial_cells.1", "Fibroblasts.1", 
                                    "Mast_cells.1", "Endothelial_cells.1", "B_lymphocytes.1", "Ependymal_cells.1", "Oligodendrocytes.1")
rownames(cells.fraction.time2) <- c("T_lymphocytes.2", "Dendritic_cells.2", "NK_cells.2", "Epithelial_cells.2", "Fibroblasts.2", 
                                    "Mast_cells.2", "Endothelial_cells.2", "B_lymphocytes.2", "Ependymal_cells.2", "Oligodendrocytes.2")
rownames(cells.fraction.time3) <- c("T_lymphocytes.3", "Dendritic_cells.3", "NK_cells.3", "Epithelial_cells.3", "Fibroblasts.3", 
                                    "Mast_cells.3", "Endothelial_cells.3", "B_lymphocytes.3", "Ependymal_cells.3", "Oligodendrocytes.3")
cell.fraction.patient1 <- rbind(cells.fraction.time1[1,], cells.fraction.time2[1,], cells.fraction.time3[1,],
                       cells.fraction.time1[2,], cells.fraction.time2[2,], cells.fraction.time3[2,],
                       cells.fraction.time1[3,], cells.fraction.time2[3,], cells.fraction.time3[3,],
                       cells.fraction.time1[4,], cells.fraction.time2[4,], cells.fraction.time3[4,],
                       cells.fraction.time1[5,], cells.fraction.time2[5,], cells.fraction.time3[5,],
                       cells.fraction.time1[6,], cells.fraction.time2[6,], cells.fraction.time3[6,],
                       cells.fraction.time1[7,], cells.fraction.time2[7,], cells.fraction.time3[7,],
                       cells.fraction.time1[8,], cells.fraction.time2[8,], cells.fraction.time3[8,],
                       cells.fraction.time1[9,], cells.fraction.time2[9,], cells.fraction.time3[9,],
                       cells.fraction.time1[10,], cells.fraction.time2[10,], cells.fraction.time3[10,])
save(cell.fraction.patient1, file = "cell.fraction.patient1.rdata")

A <- c("T_lymphocytes.1", "T_lymphocytes.2", "T_lymphocytes.3", 
       "Dendritic_cells.1", "Dendritic_cells.2","Dendritic_cells.3",
       "NK_cells.1", "NK_cells.2", "NK_cells.3",
       "Epithelial_cells.1","Epithelial_cells.2","Epithelial_cells.3",
       "Fibroblasts.1", "Fibroblasts.2", "Fibroblasts.3",
       "Mast_cells.1", "Mast_cells.2", "Mast_cells.3", 
       "Endothelial_cells.1", "Endothelial_cells.2","Endothelial_cells.3",
       "B_lymphocytes.1", "B_lymphocytes.2", "B_lymphocytes.3", 
       "Ependymal_cells.1", "Ependymal_cells.2", "Ependymal_cells.3",
       "Oligodentrocytes.1", "Oligodentrocytes.2", "Oligodentrocytes.3",
       "T_lymphocytes.1", "T_lymphocytes.2", "T_lymphocytes.3", 
       "Dendritic_cells.1", "Dendritic_cells.2","Dendritic_cells.3",
       "NK_cells.1", "NK_cells.2", "NK_cells.3",
       "Epithelial_cells.1","Epithelial_cells.2","Epithelial_cells.3",
       "Fibroblasts.1", "Fibroblasts.2", "Fibroblasts.3",
       "Mast_cells.1", "Mast_cells.2", "Mast_cells.3", 
       "Endothelial_cells.1", "Endothelial_cells.2","Endothelial_cells.3",
       "B_lymphocytes.1", "B_lymphocytes.2", "B_lymphocytes.3", 
       "Ependymal_cells.1", "Ependymal_cells.2", "Ependymal_cells.3",
       "Oligodentrocytes.1", "Oligodentrocytes.2", "Oligodentrocytes.3")

B <- c(cell.fraction.patient1$non_tumor, cell.fraction$tumor)

D <- c("non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant",
       "non_maligant","non_maligant", "non_maligant","non_maligant", "non_maligant","non_maligant","non_maligant","non_maligant",
       "non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant",
       "non_maligant","non_maligant","non_maligant", "non_maligant","non_maligant","non_maligant",
       "tumor", "tumor", "tumor","tumor", "tumor", "tumor","tumor", "tumor", "tumor","tumor", "tumor", "tumor","tumor", "tumor", "tumor",
       "tumor", "tumor", "tumor","tumor", "tumor", "tumor","tumor", "tumor", "tumor","tumor", "tumor", "tumor", "tumor", "tumor", "tumor")

cell.fraction.data <- data.frame(A, B, D)
names(cell.fraction.data) <- c("Cell", "Percentage", "Condition")

p <- ggplot(cell.fraction.data, aes(x=Cell, y=Percentage, fill=Condition)) +
  geom_bar (position = "fill", stat="identity") +
  scale_fill_manual(values = c("darkolivegreen","firebrick3")) +
  coord_flip()
p + theme_classic()
save(cell.fraction.data, file = "cell.fraction.data.rdata")

#####-------------------------- patient 2

rownames(cells.fraction.time1) <- c("T_lymphocytes.1", "Dendritic_cells.1", "NK_cells.1", "Epithelial_cells.1", "Fibroblasts.1", 
                                    "Mast_cells.1", "Endothelial_cells.1", "B_lymphocytes.1")
rownames(cells.fraction.time2) <- c("T_lymphocytes.2", "Dendritic_cells.2", "NK_cells.2", "Epithelial_cells.2", "Fibroblasts.2", 
                                    "Mast_cells.2", "Endothelial_cells.2", "B_lymphocytes.2")
rownames(cells.fraction.time3) <- c("T_lymphocytes.3", "Dendritic_cells.3", "NK_cells.3", "Epithelial_cells.3", "Fibroblasts.3", 
                                    "Mast_cells.3", "Endothelial_cells.3", "B_lymphocytes.3")
cell.fraction.patient2 <- rbind(cells.fraction.time1[1,], cells.fraction.time2[1,], cells.fraction.time3[1,],
                       cells.fraction.time1[2,], cells.fraction.time2[2,], cells.fraction.time3[2,],
                       cells.fraction.time1[3,], cells.fraction.time2[3,], cells.fraction.time3[3,],
                       cells.fraction.time1[4,], cells.fraction.time2[4,], cells.fraction.time3[4,],
                       cells.fraction.time1[5,], cells.fraction.time2[5,], cells.fraction.time3[5,],
                       cells.fraction.time1[6,], cells.fraction.time2[6,], cells.fraction.time3[6,],
                       cells.fraction.time1[7,], cells.fraction.time2[7,], cells.fraction.time3[7,],
                       cells.fraction.time1[8,], cells.fraction.time2[8,], cells.fraction.time3[8,])
save(cell.fraction.patient2, file = "cell.fraction.patient2.rdata")


A <- c("T_lymphocytes.1", "T_lymphocytes.2", "T_lymphocytes.3", 
       "Dendritic_cells.1", "Dendritic_cells.2","Dendritic_cells.3",
       "NK_cells.1", "NK_cells.2", "NK_cells.3",
       "Epithelial_cells.1","Epithelial_cells.2","Epithelial_cells.3",
       "Fibroblasts.1", "Fibroblasts.2", "Fibroblasts.3",
       "Mast_cells.1", "Mast_cells.2", "Mast_cells.3", 
       "Endothelial_cells.1", "Endothelial_cells.2","Endothelial_cells.3",
       "B_lymphocytes.1", "B_lymphocytes.2", "B_lymphocytes.3", 
       "T_lymphocytes.1", "T_lymphocytes.2", "T_lymphocytes.3", 
       "Dendritic_cells.1", "Dendritic_cells.2","Dendritic_cells.3",
       "NK_cells.1", "NK_cells.2", "NK_cells.3",
       "Epithelial_cells.1","Epithelial_cells.2","Epithelial_cells.3",
       "Fibroblasts.1", "Fibroblasts.2", "Fibroblasts.3",
       "Mast_cells.1", "Mast_cells.2", "Mast_cells.3", 
       "Endothelial_cells.1", "Endothelial_cells.2","Endothelial_cells.3",
       "B_lymphocytes.1", "B_lymphocytes.2", "B_lymphocytes.3")

B <- c(cell.fraction.patient2$non_tumor, cell.fraction$tumor)

D <- c("non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant",
       "non_maligant","non_maligant", "non_maligant","non_maligant", "non_maligant","non_maligant","non_maligant","non_maligant",
       "non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant",
       "tumor", "tumor", "tumor","tumor", "tumor", "tumor","tumor", "tumor", "tumor","tumor", "tumor", "tumor","tumor", "tumor", "tumor",
       "tumor", "tumor", "tumor","tumor", "tumor", "tumor","tumor", "tumor", "tumor")

cell.fraction.data <- data.frame(A, B, D)
names(cell.fraction.data) <- c("Cell", "Percentage", "Condition")
p <- ggplot(cell.fraction.data, aes(x=Cell, y=Percentage, fill=Condition)) +
  geom_bar (position = "fill", stat="identity") +
  scale_fill_manual(values = c("darkolivegreen","firebrick3")) +
  coord_flip()
p + theme_classic()
save(cell.fraction.data, file = "cell.fraction.data.rdata")



#####-------------------------- patient 3

rownames(cells.fraction.time1) <- c("T_lymphocytes.1", "Dendritic_cells.1", "NK_cells.1", "Epithelial_cells.1", "Fibroblasts.1", 
                                    "Mast_cells.1", "B_lymphocytes.1", "Ependymal_cells.1", "Oligodentrocytes.1")
rownames(cells.fraction.time2) <- c("T_lymphocytes.2", "Dendritic_cells.2", "NK_cells.2", "Epithelial_cells.2", "Fibroblasts.2", 
                                    "Mast_cells.2", "B_lymphocytes.2", "Ependymal_cells.2", "Oligodentrocytes.2")
rownames(cells.fraction.time3) <- c("T_lymphocytes.3", "Dendritic_cells.3", "NK_cells.3", "Epithelial_cells.3", "Fibroblasts.3", 
                                    "Mast_cells.3", "B_lymphocytes.3", "Ependymal_cells.3", "Oligodentrocytes.3")
cell.fraction.patient3 <- rbind(cells.fraction.time1[1,], cells.fraction.time2[1,], cells.fraction.time3[1,],
                       cells.fraction.time1[2,], cells.fraction.time2[2,], cells.fraction.time3[2,],
                       cells.fraction.time1[3,], cells.fraction.time2[3,], cells.fraction.time3[3,],
                       cells.fraction.time1[4,], cells.fraction.time2[4,], cells.fraction.time3[4,],
                       cells.fraction.time1[5,], cells.fraction.time2[5,], cells.fraction.time3[5,],
                       cells.fraction.time1[6,], cells.fraction.time2[6,], cells.fraction.time3[6,],
                       cells.fraction.time1[7,], cells.fraction.time2[7,], cells.fraction.time3[7,],
                       cells.fraction.time1[8,], cells.fraction.time2[8,], cells.fraction.time3[8,],
                       cells.fraction.time1[9,], cells.fraction.time2[9,], cells.fraction.time3[9,])
save(cell.fraction.patient3, file = "cell.fraction.patient3.rdata")

cell.fraction

A <- c("T_lymphocytes.1", "T_lymphocytes.2", "T_lymphocytes.3", 
       "Dendritic_cells.1", "Dendritic_cells.2","Dendritic_cells.3",
       "NK_cells.1", "NK_cells.2", "NK_cells.3",
       "Epithelial_cells.1","Epithelial_cells.2","Epithelial_cells.3",
       "Fibroblasts.1", "Fibroblasts.2", "Fibroblasts.3",
       "Mast_cells.1", "Mast_cells.2", "Mast_cells.3", 
       "B_lymphocytes.1", "B_lymphocytes.2", "B_lymphocytes.3", 
       "Ependymal_cells.1", "Ependymal_cells.2", "Ependymal_cells.3",
       "Oligodentrocytes.1", "Oligodentrocytes.2", "Oligodentrocytes.3",
       
       "T_lymphocytes.1", "T_lymphocytes.2", "T_lymphocytes.3", 
       "Dendritic_cells.1", "Dendritic_cells.2","Dendritic_cells.3",
       "NK_cells.1", "NK_cells.2", "NK_cells.3",
       "Epithelial_cells.1","Epithelial_cells.2","Epithelial_cells.3",
       "Fibroblasts.1", "Fibroblasts.2", "Fibroblasts.3",
       "Mast_cells.1", "Mast_cells.2", "Mast_cells.3", 
       "B_lymphocytes.1", "B_lymphocytes.2", "B_lymphocytes.3", 
       "Ependymal_cells.1", "Ependymal_cells.2", "Ependymal_cells.3",
       "Oligodentrocytes.1", "Oligodentrocytes.2", "Oligodentrocytes.3")

B <- c(cell.fraction$non_tumor, cell.fraction$tumor)

D <- c("non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant",
       "non_maligant","non_maligant", "non_maligant","non_maligant", "non_maligant","non_maligant","non_maligant","non_maligant",
       "non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant","non_maligant",
       "non_maligant","non_maligant","non_maligant",
       "tumor", "tumor", "tumor","tumor", "tumor", "tumor","tumor", "tumor", "tumor","tumor", "tumor", "tumor","tumor", "tumor", "tumor",
       "tumor", "tumor", "tumor","tumor", "tumor", "tumor","tumor", "tumor", "tumor","tumor", "tumor", "tumor")

cell.fraction.data <- data.frame(A, B, D)
names(cell.fraction.data) <- c("Cell", "Percentage", "Condition")
p <- ggplot(cell.fraction.data, aes(x=Cell, y=Percentage, fill=Condition)) +
  geom_bar (position = "fill", stat="identity") +
  scale_fill_manual(values = c("darkolivegreen","firebrick3")) +
  coord_flip()
p + theme_classic()
save(cell.fraction.data, file = "cell.fraction.data.rdata")








