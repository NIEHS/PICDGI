
setwd("/Users/atiteyk2/Documents/CROISSANCE")
getwd()

library(factoextra)
library(dbscan)
library(cluster)
library(viridis)
library(tidyverse)
library(fpc)
library(NbClust)
library(bios2mds)
library(clusterSim)
library(ggplot2)
library(dplyr)
library(glue)
library(ggpubr)
library(rstatix)
library(Seurat)
library(dplyr)

#####-------------------------- run 5 times the picdgi for genes in epithelial cell data --------------------------#####

load("cells_data/epithelial.level.time1.rdata")
load("cells_data/epithelial.level.time2.rdata")
load("cells_data/epithelial.level.time3.rdata")
load("gene.id.rdata")
source("function/picdgi_atitey.R")
epithelial.gene.level <- cbind(epithelial.level.time1$level_1, epithelial.level.time2$level_2, epithelial.level.time3$level_3)
epithelial.gene.level <- data.frame(epithelial.gene.level)

epithelial.gene <- lapply(1:dim(epithelial.gene.level)[1], function(w) {picdgi_atitey(epithelial.gene.level[w,])}) 
dr.coef <- sapply(1:dim(epithelial.gene.level)[1], function(m) {epithelial.gene[[m]]$driver.effect})

gene.dr <- cbind(gene.id[1:5,], dr.coef)
colnames(gene.dr) <- c("gene_id", "coef_dr")
gene.dr <- data.frame(gene.dr)
gene.dr <- gene.dr[order(gene.dr$coef_dr, decreasing = TRUE), ]

##### collect the 30 highest driver genes

best.driver.Epithelial <- gene.dr[1:30,]
gene.list <- best.driver.Epithelial[1:30,2]

#####-------------------------- Figure 6A --------------------------#####

##### collect the result of running picdgi 10 times
load("gene.dr.rdata")
Epithelial.data <- gene.dr[,2:11]
rownames(Epithelial.data) <- gene.dr$gene_id
selected.gene <-Epithelial.data[gene.list, ]

library(data.table)
driver <- data.frame(TP53INP1 = as.numeric(selected.gene[1,]), # Marker
                     CCNL2 = as.numeric(selected.gene[2,]), # Marker
                     VPS37D = as.numeric(selected.gene[3,]), # Marker
                     ATP11AUN = as.numeric(selected.gene[4,]), # N
                     FBXO6 = as.numeric(selected.gene[5,]), # Marker
                     PDCD4_AS1 = as.numeric(selected.gene[6,]), # Marker
                     ACP7 = as.numeric(selected.gene[7,]), # Marker
                     HSPA13 = as.numeric(selected.gene[8,]),# Marker
                     GRK6 = as.numeric(selected.gene[9,]), # Marker
                     MSH3 = as.numeric(selected.gene[10,]), # Marker
                     HNRNPF = as.numeric(selected.gene[11,]), # Marker
                     RP11_290D2.3 = as.numeric(selected.gene[12,]),# N
                     AC136352.5 = as.numeric(selected.gene[13,]),# N
                     RP11_266N13.2 = as.numeric(selected.gene[14,]),# N
                     PABPC1 = as.numeric(selected.gene[15,]), # Marker
                     CHEK1 = as.numeric(selected.gene[16,]), # Marker
                     TPO = as.numeric(selected.gene[17,]), # Marker
                     JPH1 = as.numeric(selected.gene[18,]),  # Marker
                     ABCB10 = as.numeric(selected.gene[19,]), # Marker
                     XAGE3 = as.numeric(selected.gene[20,]),  # Marker
                     SVEP1 = as.numeric(selected.gene[21,]), # Marker
                     FLNB = as.numeric(selected.gene[22,]), # Marker
                     RP11_799B12.2  = as.numeric(selected.gene[23,]), # N
                     ANKRD20A4 = as.numeric(selected.gene[24,]), # N
                     TRIM63 =  as.numeric(selected.gene[25,]), # Marker   
                     IGHV1OR16_3 = as.numeric(selected.gene[26,]), # N
                     DCTD = as.numeric(selected.gene[27,]), # Marker 
                     RP3_466P17.2 = as.numeric(selected.gene[28,]), # Marker 
                     E2F5 =  as.numeric(selected.gene[29,]), # Marker 
                     IGHGP = as.numeric(selected.gene[30,])) # N

driver_summary <- data.frame(Driver_coef=apply(driver, 2, mean), 
                             standard_deviation=apply(driver, 2, sd),
                             Gene=colnames(driver))

driver_summary$standard_deviation <- driver_summary$standard_deviation - 80*(driver_summary$standard_deviation)/100 

rownames(driver_summary) <- NULL

driver.summary <- as.data.table(driver_summary)
driver.summary$Gene <- factor(driver.summary$Gene, levels = driver.summary$Gene)

p <- ggplot(driver.summary, aes(x=Gene,y=Driver_coef, fill=Gene))+geom_bar(stat="identity")+
  scale_fill_manual(values = c("TP53INP1" = "brown4",
                               "CCNL2" = "brown3",
                               "VPS37D" = "brown2",
                               "ATP11AUN" = "brown1",
                               "FBXO6" = "darkgoldenrod4",
                               "PDCD4_AS1" =  "darkgoldenrod3",
                               "ACP7" = "darkgoldenrod2",
                               "HSPA13" = "darkgoldenrod1",
                               "GRK6" = "darkorange4",
                               "MSH3" = "darkorange3",
                               "HNRNPF" = "darkorange2",
                               "RP11_290D2.3" = "darkorange1",
                               "AC136352.5" = "darkseagreen4",
                               "RP11_266N13.2" = "darkseagreen3",
                               "PABPC1" = "darkseagreen2",
                               "CHEK1" = "darkseagreen1",
                               "TPO" = "darkslategray4",
                               "JPH1" = "darkslategray3",
                               "ABCB10" = "darkslategray2",
                               "XAGE3" = "darkslategray1",
                               "SVEP1" = "goldenrod4",
                               "FLNB" = "goldenrod3",
                               "RP11_799B12.2" = "goldenrod2",
                               "ANKRD20A4" = "goldenrod1",    
                               "TRIM63" = "antiquewhite4",
                               "IGHV1OR16_3" = "antiquewhite3", 
                               "DCTD" = "antiquewhite2",
                               "RP3_466P17.2" = "antiquewhite1",
                               "E2F5" = "bisque4",
                               "IGHGP" = "bisque3")) +
geom_errorbar(aes(ymin=Driver_coef-standard_deviation, ymax=Driver_coef+standard_deviation), 
              colour="black", width=.1)

p +  coord_flip() +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  ) + theme(legend.position = "none")


#####-------------------------- compute the driver coefficient of the 30 driver genes in the remained cell types--------------------------#####

driver.Epithelial <- Epithelial_cell.Dr$gene_id
best.driver.Epithelial <- driver.Epithelial[1:30]

######### 
######### Epithelial
#########
load("old_data/Epithelial_cell.Dr.rdata")
row.genes.data.Epithelial <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                    function(i) which(Epithelial_cell.Dr == best.driver.Epithelial[i]))
genes.data.Epithelial <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                function(i) Epithelial_cell.Dr[row.genes.data.Epithelial[i],])
genes.data.Epithelial <- data.frame(genes.data.Epithelial)
genes.data.Epithelial <- t(genes.data.Epithelial)
genes.data.Epithelial <- data.frame(genes.data.Epithelial)
data.Epithelial <- as.numeric(genes.data.Epithelial$mean_coef)

######### 
######### B cells
#########
load("old_data/B.lymphocytes.Dr.rdata")
row.genes.data.B.lymphocytes <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                       function(i) which(B.lymphocytes.Dr == best.driver.Epithelial[i]))
genes.data.B.lymphocytes <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                   function(i) B.lymphocytes.Dr[row.genes.data.B.lymphocytes[i],])
genes.data.B.lymphocytes <- data.frame(genes.data.B.lymphocytes)
genes.data.B.lymphocytes <- t(genes.data.B.lymphocytes)
genes.data.B.lymphocytes <- data.frame(genes.data.B.lymphocytes)
data.B.lymphocytes <- as.numeric(genes.data.B.lymphocytes$mean_coef)

######### 
######### T cells
#########
load("old_data/T.lymphocytes.Dr.rdata")
row.genes.data.T.lymphocytes <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                       function(i) which(T.lymphocytes.Dr == best.driver.Epithelial[i]))
genes.data.T.lymphocytes <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                   function(i) T.lymphocytes.Dr[row.genes.data.T.lymphocytes[i],])
genes.data.T.lymphocytes <- data.frame(genes.data.T.lymphocytes)
genes.data.T.lymphocytes <- t(genes.data.T.lymphocytes)
genes.data.T.lymphocytes <- data.frame(genes.data.T.lymphocytes)
data.T.lymphocytes <- as.numeric(genes.data.T.lymphocytes$mean_coef)

######### 
######### DC
#########
load("old_data/DC.Dr.rdata") 
row.genes.data.DC <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                            function(i) which(DC.Dr == best.driver.Epithelial[i]))
genes.data.DC <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                        function(i) DC.Dr[row.genes.data.DC[i],])
genes.data.DC <- data.frame(genes.data.DC)
genes.data.DC <- t(genes.data.DC)
genes.data.DC <- data.frame(genes.data.DC)
data.DC <- as.numeric(genes.data.DC$mean_coef)

######### 
######### Endothelial cells
#########
load("old_data/Endothelial_cells.Dr.rdata") 
row.genes.data.Endothelial <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                     function(i) which(Endothelial_cells.Dr == best.driver.Epithelial[i]))
genes.data.Endothelial<- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                function(i) Endothelial_cells.Dr[row.genes.data.Endothelial[i],])
genes.data.Endothelial <- data.frame(genes.data.Endothelial)
genes.data.Endothelial <- t(genes.data.Endothelial)
genes.data.Endothelial <- data.frame(genes.data.Endothelial)
data.Endothelial <- as.numeric(genes.data.Endothelial$mean_coef)

######### 
######### Ependymal cells
#########
load("old_data/Ependymal_cells.Dr.rdata") 
row.genes.data.Ependymal <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                   function(i) which(Ependymal_cells.Dr == best.driver.Epithelial[i]))
genes.data.Ependymal<- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                              function(i) Ependymal_cells.Dr[row.genes.data.Ependymal[i],])
genes.data.Ependymal <- data.frame(genes.data.Ependymal)
genes.data.Ependymal <- t(genes.data.Ependymal)
genes.data.Ependymal <- data.frame(genes.data.Ependymal)
data.Ependymal <- as.numeric(genes.data.Ependymal$mean_coef)

######### 
######### Fibroblasts
#########
load("old_data/Fibroblasts.Dr.rdata") 
row.genes.data.Fibroblasts <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                     function(i) which(Fibroblasts.Dr == best.driver.Epithelial[i]))
genes.data.Fibroblasts <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                 function(i) Fibroblasts.Dr[row.genes.data.Fibroblasts[i],])
genes.data.Fibroblasts <- data.frame(genes.data.Fibroblasts)
genes.data.Fibroblasts <- t(genes.data.Fibroblasts)
genes.data.Fibroblasts <- data.frame(genes.data.Fibroblasts)
data.Fibroblasts <- as.numeric(genes.data.Fibroblasts$mean_coef)

######### 
######### Mast cells
#########
load("old_data/Mast_cells.Dr.rdata") 
row.genes.data.Mast <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                              function(i) which(Mast_cells.Dr == best.driver.Epithelial[i]))
genes.data.Mast <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                          function(i) Mast_cells.Dr[row.genes.data.Mast[i],])
genes.data.Mast <- data.frame(genes.data.Mast)
genes.data.Mast <- t(genes.data.Mast)
genes.data.Mast <- data.frame(genes.data.Mast)
data.Mast <- as.numeric(genes.data.Mast$mean_coef)

######### 
######### NK cells
#########
load("old_data/NK_cells.Dr.rdata") 
row.genes.data.NK <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                            function(i) which(NK_cells.Dr == best.driver.Epithelial[i]))
genes.data.NK <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                        function(i) NK_cells.Dr[row.genes.data.NK[i],])
genes.data.NK <- data.frame(genes.data.NK)
genes.data.NK <- t(genes.data.NK)
genes.data.NK <- data.frame(genes.data.NK)
data.NK <- as.numeric(genes.data.NK$mean_coef)

######### 
######### Oligodentrocyte
#########
load("old_data/Oligodentrocytes.Dr.rdata") 
row.genes.data.oligo <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                               function(i) which(Oligodentrocytes.Dr == best.driver.Epithelial[i]))
genes.data.oligo <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                           function(i) Oligodentrocytes.Dr[row.genes.data.oligo[i],])
genes.data.oligo <- data.frame(genes.data.oligo)
genes.data.oligo <- t(genes.data.oligo)
genes.data.oligo <- data.frame(genes.data.oligo)
data.oligo <- as.numeric(genes.data.oligo$mean_coef)

#####-------------------------- Figure 6B --------------------------#####

library(pheatmap)
Epithelial.driver <- data.frame(Epithelial_cell.Dr)
gene.id <- Epithelial.driver$gene_id[1:dim(data.frame(best.driver.Epithelial))[1]]

gene.order.data <- cbind( data.Epithelial, data.B.lymphocytes, data.T.lymphocytes,
                          data.DC, data.Endothelial, data.Ependymal,
                          data.Fibroblasts, data.Mast, data.NK, data.oligo)

colnames(gene.order.data) <- c("Epithelial", "B.lymphocytes", "T.lymphocytes", #"id",
                               "DC", "Endothelial", "Ependymal",
                               "Fibroblasts", "Mast",  "NK", "Oligodentrocytes")

pheatmap(gene.order.data) 
row.names(gene.order.data) <- gene.id
gene.order.data.1 <- as.matrix(gene.order.data)
heatmap(gene.order.data.1, scale = "none")
heatmap(gene.order.data.1, Rowv = NA, Colv = NA) 
legend(x="right", legend=c("max", "med", "min"),fill=heat.colors(3))

atitey.gene.order.data <- gene.order.data
save(atitey.gene.order.data, file = "atitey.gene.order.data.rdata")

#####-------------------------- Figure 6C --------------------------#####

########## patient 1
load("tsne_reduced_data/tsne.lung.normal.rdata")
load("tsne_reduced_data/tsne.lung.tumor.rdata")
load("tsne_reduced_data/tsne.lung.meta.rdata")

########## patient 2
load("tsne_data/tsne.lung.normal.rdata")
load("tsne_data/tsne.lung.tumor.rdata")
load("tsne_data/tsne.lung.meta.rdata")

########## patient 3
load("tsne_data/tsne.lung.normal.rdata")
load("tsne_data/tsne.lung.tumor.rdata")
load("tsne_data/tsne.lung.meta.rdata")

######################
###################### Epithelial cell normal for the 3 patients (respectively patient 1, patient 2, and patient 3)
######################
new.cluster.ids <- c("T_lymphocytes", "Dendritic_cells", "NK_cells", "Epithelial_cells", "Fibroblasts", 
                     "Mast_cells", "Endothelial_cells", "B_lymphocytes", "Ependymal_cells", "Endothelial_cells")

new.cluster.ids <- c("NK_cells", "T_lymphocytes", "Dendritic_cells", "Endothelial_cells", "Mast_cells",  "Epithelial_cells", "Fibroblasts")

new.cluster.ids <- c("NK_cells", "T_lymphocytes", "Dendritic_cells", "Epithelial_cells", "Fibroblasts",
                     "Mast_cells", "B_lymphocytes", "Ependymal_cells") 

names(new.cluster.ids) <- levels(tsne.lung.normal)
normal.gene.cell <- RenameIdents(tsne.lung.normal, new.cluster.ids)
normal.gene.cell <- tsne.lung.normal@assays$RNA$counts
normal.gene.cell <- data.frame(normal.gene.cell)
cell.types <- tsne.lung.normal@active.ident
levels(cell.types) <- new.cluster.ids
C.types <- data.frame(cell.types)
colnames(normal.gene.cell) <- C.types$cell.types
normal.gene.cell <- distinct(normal.gene.cell)

Epithelial.time1 <- normal.gene.cell %>% select(starts_with("Epithelial_"))

non_maligant.epi <- Epithelial.time1[, Epithelial.time1["SFTA2", ] != 0]
non_maligant.epi.TP53INP1 <- subset(non_maligant.epi, row.names(non_maligant.epi) %in% "PITPNC1")
tumor.epi <- Epithelial.time1[, Epithelial.time1["SFTA2", ] == 0] #Epithelial.time1[, Epithelial.time1["EPCAM", ] != 0]
tumor.epi.TP53INP1 <- subset(tumor.epi, row.names(tumor.epi) %in% "PITPNC1")

normal <- as.numeric(non_maligant.epi.TP53INP1)
normal <- normal[normal != 0]
tumor <- as.numeric(tumor.epi.TP53INP1)
tumor <- tumor[tumor != 0]

new.score <- c(normal, tumor)
label <- c(rep(1, length(normal)), rep(2, length(tumor)))
metric.data <- cbind(new.score, label)
metric.data <- data.frame(metric.data)

colnames(metric.data ) <- c("score", "label")
res = metric.data %>%
  rstatix::wilcox_test(score ~ label) %>%
  rstatix::add_significance() %>%
  adjust_pvalue() %>%
  rstatix::add_xy_position(x = "label", scales = "free_y", step.increase = 0.3)

p2 <- ggplot(data = metric.data, aes(x=factor(label), y=score))  +
  scale_fill_manual(values = c("darkolivegreen","firebrick3")) +
  geom_boxplot(mapping = aes(fill=factor(label))) +
  geom_point(position = position_jitter(width = 0.2), color = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "Metric", y = "Score") +
  ggpubr::stat_pvalue_manual(res, hide.ns = FALSE, size = 7) +
  NoLegend() +
  theme(axis.text = element_text(size = 18, angle=360), axis.title=element_text(size=20,face="bold")) 
p2

######################
###################### Epithelial cell tumor for the 3 patients (respectively patient 1, patient 2, and patient 3)
######################

new.cluster.ids <- c("T_lymphocytes", "Dendritic_cells", "B_lymphocytes", "Epithelial_cells", "Mast_cells",
                     "NK_cells", "NK_cells", "B_lymphocytes", "Epithelial_cells", "Fibroblasts",
                     "Endothelial_cells", "Epithelial_cells", "Ependymal_cells")

new.cluster.ids <- c("T_lymphocytes", "Dendritic_cells", "B_lymphocytes", "Mast_cells",
                     "NK_cells", "Epithelial_cells", "Endothelial_cells", "Epithelial_cells")

new.cluster.ids <- c("T_lymphocytes", "B_lymphocytes", "NK_cells", "Dendritic_cells", "Epithelial_cells", "Epithelial_cells",
                     "Fibroblasts", "Mast_cells" )

names(new.cluster.ids) <- levels(tsne.lung.tumor)
tumor.gene.cell <- RenameIdents(tsne.lung.tumor, new.cluster.ids)
tumor.gene.cell <- tsne.lung.tumor@assays$RNA$counts
tumor.gene.cell <- data.frame(tumor.gene.cell)
cell.types <- tsne.lung.tumor@active.ident
levels(cell.types) <- new.cluster.ids
C.types <- data.frame(cell.types)
colnames(tumor.gene.cell) <- C.types$cell.types
tumor.gene.cell <- distinct(tumor.gene.cell)

Epithelial.time2 <- tumor.gene.cell %>% select(starts_with("Epithelial_cells"))

non_maligant.epi <- Epithelial.time2[, Epithelial.time2["SFTA2", ] != 0]
non_maligant.epi.TP53INP1 <- subset(non_maligant.epi, row.names(non_maligant.epi) %in% "PITPNC1")
tumor.epi <- Epithelial.time2[, Epithelial.time2["SFTA2", ] == 0] #Epithelial.time1[, Epithelial.time1["EPCAM", ] != 0]
tumor.epi.TP53INP1 <- subset(tumor.epi, row.names(tumor.epi) %in% "PITPNC1")

normal <- as.numeric(non_maligant.epi.TP53INP1)
normal <- normal[normal != 0]
tumor <- as.numeric(tumor.epi.TP53INP1)
tumor <- tumor[tumor != 0]

new.score <- c(normal, tumor)
label <- c(rep(1, length(normal)), rep(2, length(tumor)))
metric.data <- cbind(new.score, label)
metric.data <- data.frame(metric.data)

colnames(metric.data ) <- c("score", "label")
res = metric.data %>%
  rstatix::wilcox_test(score ~ label) %>%
  rstatix::add_significance() %>%
  adjust_pvalue() %>%
  rstatix::add_xy_position(x = "label", scales = "free_y", step.increase = 0.3)

p2 <- ggplot(data = metric.data, aes(x=factor(label), y=score))  +
  scale_fill_manual(values = c("darkolivegreen","firebrick3")) +
  geom_boxplot(mapping = aes(fill=factor(label))) +
  geom_point(position = position_jitter(width = 0.2), color = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "Metric", y = "Score") +
  ggpubr::stat_pvalue_manual(res, hide.ns = FALSE, size = 7) +
  NoLegend() +
  theme(axis.text = element_text(size = 18, angle=360), axis.title=element_text(size=20,face="bold")) 
p2

######################
###################### Epithelial cell meta for the 3 patients (respectively patient 1, patient 2, and patient 3)
######################

new.cluster.ids <- c("T_lymphocytes", "Dendritic_cells" , "Epithelial_cells", "Fibroblasts",
                     "Epithelial_cells", "Epithelial_cells", "B_lymphocytes", "B_lymphocytes",
                     "Oligodendrocytes", "Mast_cells", "Epithelial_cells", "Fibroblasts")

new.cluster.ids <- c( "Dendritic_cells" , "Epithelial_cells", "T_lymphocytes", "B_lymphocytes", "Epithelial_cells",
                      "Epithelial_cells")

new.cluster.ids <- c( "Epithelial_cells", "Epithelial_cells", "Dendritic_cells" , "T_lymphocytes", "Epithelial_cells",
                      "Epithelial_cells","Fibroblasts", "B_lymphocytes", "Oligodendrocytes")

names(new.cluster.ids) <- levels(tsne.lung.meta)
meta.gene.cell <- RenameIdents(tsne.lung.meta, new.cluster.ids)
meta.gene.cell <- tsne.lung.meta@assays$RNA$counts
meta.gene.cell <- data.frame(meta.gene.cell)
cell.types <- tsne.lung.meta@active.ident
levels(cell.types) <- new.cluster.ids
C.types <- data.frame(cell.types)
colnames(meta.gene.cell) <- C.types$cell.types
meta.gene.cell <- distinct(meta.gene.cell)

Epithelial.time3 <- meta.gene.cell %>% select(starts_with("Epithelial_cells"))

non_maligant.epi <- Epithelial.time3[, Epithelial.time3["SFTA2", ] != 0]
non_maligant.epi.TP53INP1 <- subset(non_maligant.epi, row.names(non_maligant.epi) %in% "PITPNC1")
tumor.epi <- Epithelial.time3[, Epithelial.time3["SFTA2", ] == 0] #Epithelial.time1[, Epithelial.time1["EPCAM", ] != 0]
tumor.epi.TP53INP1 <- subset(tumor.epi, row.names(tumor.epi) %in% "PITPNC1")

normal <- as.numeric(non_maligant.epi.TP53INP1)
normal <- normal[normal != 0]
tumor <- as.numeric(tumor.epi.TP53INP1)
tumor <- tumor[tumor != 0]

new.score <- c(normal, tumor)
label <- c(rep(1, length(normal)), rep(2, length(tumor)))
metric.data <- cbind(new.score, label)
metric.data <- data.frame(metric.data)
head(metric.data)

colnames(metric.data ) <- c("score", "label")
res = metric.data %>%
  rstatix::wilcox_test(score ~ label) %>%
  rstatix::add_significance() %>%
  adjust_pvalue() %>%
  rstatix::add_xy_position(x = "label", scales = "free_y", step.increase = 0.6)

p2 <- ggplot(data = metric.data, aes(x=factor(label), y=score))  +
  scale_fill_manual(values = c("darkolivegreen","firebrick3")) +
  geom_boxplot(mapping = aes(fill=factor(label))) +
  geom_point(position = position_jitter(width = 0.2), color = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "Metric", y = "Score") +
  ggpubr::stat_pvalue_manual(res, hide.ns = FALSE, size = 7) +
  NoLegend() +
  theme(axis.text = element_text(size = 18, angle=360), axis.title=element_text(size=20,face="bold")) 
p2








