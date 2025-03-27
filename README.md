# PICDGI: Predict Immunosuppressive Cancer Driver Gene using gene-gene Interaction feature

$~~$

**Traditional frequency-based methods identify driver mutations by evaluating their recurrence across tumors, operating under the assumption that mutations occurring more frequently are functionally significant. However, these approaches often overlook rare driver mutations, particularly those with low mutation frequencies.**

**To overcome this limitation, we present PICDGI, a computational framework designed to predict cancer driver genes (CDGs) by leveraging single-cell RNA sequencing (scRNA-seq) data to model gene-gene interactions, thereby capturing cellular heterogeneity within tumors.**

**Unlike conventional methods, PICDGI utilizes variational Bayesian inference to model gene expression dynamics and quantify the influence of gene-gene interactions on tumor progression. It identifies CDGs by detecting their immunosuppressive effects on immune cells through interaction-driven mechanisms during cancer development.**

**Applied to LUAD scRNA-seq data, PICDGI identified top-ranked CDGs, with 62% overlapping known OGs and TSGs, while 38% are novel candidates. Compared to Moran’s I test in Monocle3, PICDGI-identified genes showed higher expression, reinforcing their role in tumor growth. By capturing gene expression dynamics and interactions, PICDGI offers insights into potential therapeutic targets for personalized cancer treatment.**

$~~$

# Steps in the PICDGI Framework for Identifying Cancer Driver Genes

**Computational objective:** Model the impact of gene interactions on other genes in driving cancer progression in patients.

![](Figure/github_2.png)

**Step 1:** PICDGI calculates the average expression level of each gene across all cells within each cluster at different time points corresponding to the scRNA-seq data collection. This generates a pseudo-time-series gene expression dataset with dimensions N × T × C, where N represents the number of genes, T the number of time points, and C the number of clusters. Then PICDGI determines the cancer cell fraction within each cell to identify cancer progenitor cells using appropriate markers.

![](Figure/github_1.png)

**Step 2:** PICDGI models gene mutation dynamics using a time-varying Autoregressive Moving Average (ARMA) process. 

**Step 3:** PICDGI applies Bayes' theorem along with variational Bayesian inference to estimate the probability density function that captures gene mutation dynamics during cancer progression, accounting for gene interactions.

**Step 4:** PICDGI calculates the gene driver coefficient as the highest density interval (HDI) of the posterior probability density function, which represents the range of the most probable true gene effects on other genes.

![](Figure/github_3.png)

$~~$

# Run PICDGI

**R Packages**<br>
libs <- c("MASS", "tidyr", "dplyr", "caret", "viridis", "magrittr", "plyr")<br>
lapply(libs, library, character.only=TRUE)<br>

########<br>
######## load time series data of the cancer progenitor cells (epithelial cell)<br>
########<br>
load("dataset/epithelial.level.time1.rdata")<br>
load("dataset/epithelial.level.time2.rdata")<br>
load("dataset/epithelial.level.time3.rdata")<br>

########<br>
########load genes id<br>
########<br>
load("gene.id.rdata")<br>

########<br>
######## call picdgi function<br>
########<br>
source("function/picdgi_atitey.R")<br>
epithelial.gene.level <- cbind(epithelial.level.time1$level_1, epithelial.level.time2$level_2, epithelial.level.time3$level_3)<br>
epithelial.gene.level <- data.frame(epithelial.gene.level)<br>

epithelial.gene <- lapply(1:dim(epithelial.gene.level)[1], function(w) {picdgi_atitey(epithelial.gene.level[w,])}) <br>
dr.coef <- sapply(1:dim(epithelial.gene.level)[1], function(m) {epithelial.gene[[m]]$driver.effect})<br>

gene.dr <- cbind(gene.id[1:5,], dr.coef)<br>
colnames(gene.dr) <- c("gene_id", "coef_dr")<br>
gene.dr <- data.frame(gene.dr)<br>

########<br>
######## order genes in decreasing order of their coefficient drivers<br>
########<br>
gene.dr <- gene.dr[order(gene.dr$coef_dr, decreasing = TRUE), ]<br>



