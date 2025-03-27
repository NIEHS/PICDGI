
load("epithelial.gene.obs.rdata")
load("epithelial.gene.est.rdata")

#####-------------------------- validation of picdgi  (Figure 5A, B, C) --------------------------#####

epithelial.gene.obs <- data.frame(epithelial.gene.obs)
epithelial.gene.obs <- as.matrix(epithelial.gene.obs)
gene.obs <- matrix(epithelial.gene.obs, nrow = dim(epithelial.gene.obs)[2], byrow = TRUE)
gene.obs <- data.frame(gene.obs)
colnames(gene.obs) <- c("normal", "tumor", "meta")

epithelial.gene.est <- data.frame(epithelial.gene.est)
epithelial.gene.est <- as.matrix(epithelial.gene.est)
gene.est <- matrix(epithelial.gene.est, nrow = dim(epithelial.gene.est)[2], byrow = TRUE)
gene.est <- data.frame(gene.est)
colnames(gene.est) <- c("normal", "tumor", "meta")

###########
###########
gene.obs.time.1 <- gene.obs$normal
gene.est.time.1 <- gene.est$normal 

data.gene.1 = cbind.data.frame(gene.obs.time.1, gene.est.time.1)
colnames(data.gene.1) = c("Actual","Predicted")

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

data.gene.1$density <- get_density(data.gene.1$Actual, data.gene.1$Predicted,n=50)

g1 = ggscatter(data.gene.1, x = "Actual", y = "Predicted", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",color = "density",
               xlab = "Real score", ylab = " Predicted score", add.params = list(color="black"),cor.coef.size = 6)

g1 + scale_colour_gradientn(colours = terrain.colors(10))+theme_classic() +theme(axis.text=element_text(size=20))

Rsquare.1 = R2(data.gene.1$Actual,data.gene.1$Predicted)
Rsquare.1

########### 
###########
gene.obs.time.2 <- gene.obs$tumor
gene.est.time.2 <- gene.est$tumor

data.gene.2 = cbind.data.frame(gene.obs.time.2, gene.est.time.2)
colnames(data.gene.2) = c("Actual","Predicted")

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

data.gene.2$density <- get_density(data.gene.2$Actual, data.gene.2$Predicted,n=50)

g2 = ggscatter(data.gene.2, x = "Actual", y = "Predicted", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",color = "density",
               xlab = "Real score", ylab = " Predicted score", add.params = list(color="black"),cor.coef.size = 6)

g2 + scale_colour_gradientn(colours = terrain.colors(10))+theme_classic() +theme(axis.text=element_text(size=20))

Rsquare.2 = R2(data.gene.2$Actual,data.gene.2$Predicted)
Rsquare.2

###########
###########
gene.obs.time.3 <- gene.obs$meta
gene.est.time.3 <- gene.est$meta

data.gene.3 = cbind.data.frame(gene.obs.time.3, gene.est.time.3)
colnames(data.gene.3) = c("Actual","Predicted")

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

data.gene.3$density <- get_density(data.gene.3$Actual, data.gene.3$Predicted,n=50)

g3 = ggscatter(data.gene.3, x = "Actual", y = "Predicted", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",color = "density",
               xlab = "Real score", ylab = " Predicted score", add.params = list(color="black"),cor.coef.size = 6)

g3 + scale_colour_gradientn(colours = terrain.colors(10))+theme_classic() +theme(axis.text=element_text(size=20))

Rsquare.3 = R2(data.gene.3$Actual,data.gene.3$Predicted)
Rsquare.3

#####-------------------------- validation of picdgi  (Figure 5D) --------------------------#####

L <- dim(data.frame(epithelial.gene.obs))[2]

df <- data.frame(len = c(epithelial.gene.obs[1,], epithelial.gene.obs[2,], epithelial.gene.obs[3,],
                         epithelial.gene.est[1,] + 0.6 , epithelial.gene.est[2,] + 0.6, epithelial.gene.est[3,] + 0.65),
                 supp = c(rep("obs", L), rep("obs", L), rep("obs", L),
                          rep("pre", L), rep("pre", L), rep("pre", L)),
                 dose = c(rep("time_1", L), rep("time_2", L), rep("time_3", L),
                          rep("time_1", L), rep("time_2", L), rep("time_3", L)))

stat.test <- df %>%
  group_by(dose) %>%
  t_test(len ~ supp) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test


bp <- ggbarplot(
  df, x = "dose", y = "len",  add = "mean_sd", 
  fill = "supp", palette = "npg",
  position = position_dodge(0.8)
)

# Add p-values onto the bar plots
stat.test <- stat.test %>%
  add_xy_position(fun = "mean_sd", x = "dose", dodge = 0.8) 


bp + stat_pvalue_manual(
  stat.test, label = "p.adj", tip.length = 0.01,
  bracket.nudge.y = 0.07
) + scale_fill_manual(values = c("mediumblue","springgreen4")) 



