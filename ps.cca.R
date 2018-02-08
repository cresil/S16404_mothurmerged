#### libraries and path ####
library("phyloseq")
library("ggplot2")
library("vegan")
library(dada2)
library(gplots)
library(vegan)

rm(list=ls())
setwd("C:/Users/mdavi/Documents/Trouw/S16404_mothurmerged")
ps <- readRDS("ps.S16404.rds")
ps <- rarefy_even_depth(ps, sample.size = 10000)
#library(MASS)
#library(vegan3d)

#bla <- lapply(data.frame(ps@sam_data),class)
#sam_data.nf <- data.frame(ps@sam_data[,!unlist(lapply(data.frame(ps@sam_data),class))=="factor"])
#sam_data.f <- data.frame(ps@sam_data[,unlist(bla)=="factor"])

vare.cca <- cca(t(ps@otu_table) ~ ., as.data.frame(unclass(sample_data(ps))))
#vare.cca <- cca(t(ps@otu_table) ~ pH +LAcol , sam_data.nf)

plot(vare.cca)
colvec <- rainbow(4)
colvec <- colvec[ps@sam_data[rownames(vare.cca$CCA$wa)]$Treat]
text(vare.cca$CCA$wa[,1],vare.cca$CCA$wa[,2],rownames(vare.cca$CCA$wa))
points(vare.cca$CCA$wa[,1:2],cex=2, pch=18, col=colvec, )
legend(x="topright",legend=c("Control","AXOS","CELL","A+C"), fill = rainbow(4))


bla <- lapply(data.frame(ps@sam_data),class)
sam_data.nf <- data.frame(ps@sam_data[,!unlist(lapply(data.frame(ps@sam_data),class))=="factor"])
sam_data.f <- data.frame(ps@sam_data[,unlist(bla)=="factor"])

vare.cca <- cca(t(ps@otu_table) ~ ., cbind(sam_data.nf, sam_data.f))



## colon


ps.c <- prune_samples(ps@sam_data$Sample.Type!="Feaces",ps)

bla <- lapply(data.frame(ps.c@sam_data),class)
sam_data.nf <- data.frame(ps.c@sam_data[,!unlist(lapply(data.frame(ps.c@sam_data),class))=="factor"])
sam_data.f <- data.frame(ps.c@sam_data[,unlist(bla)=="factor"])

vare.cca <- cca(t(ps.c.c@otu_table) ~ ., as.data.frame(unclass(sample_data(ps.c.c))))
vare.cca <- cca(t(ps.c@otu_table) ~ ., cbind(sam_data.nf, sam_data.f)) 
plot(vare.cca)
colvec <- rainbow(4)
colvec <- colvec[ps.c@sam_data[rownames(vare.cca$CCA$wa)]$Treat]
#text(vare.cca$CCA$wa[,1],vare.cca$CCA$wa[,2],rownames(vare.cca$CCA$wa))
points(vare.cca$CCA$wa[,1:2],cex=2, pch=18, col=colvec, )
legend(x="right",legend=c(1:4), fill = rainbow(4))



#### rel change ####


ps.cea <- prune_samples(ps@sam_data$Sample.Type!="Feaces",ps)
ps.col <- prune_samples(ps@sam_data$Sample.Type=="Feaces",ps)
ps.rel <- ps.cea


colotu <- as.matrix(unclass(ps.col@otu_table[,match(as.character(ps.col@sam_data$animal.number), as.character(ps.cea@sam_data$animal.number))]))
ceaotu <- as.matrix(unclass(ps.cea@otu_table))

attributes(colotu) <- NULL 
colotu <- matrix(colotu, ncol = 40)

attributes(ceaotu) <- NULL 
ceaotu <- matrix(ceaotu, ncol = 40)



(ps.cea@otu_table + ps.col@otu_table[match(as.character(ps.col@sam_data$animal.number), as.character(ps.cea@sam_data$animal.number)),])

oturel <- colotu/(ceaotu + colotu)
colnames(oturel) <- sample_names(ps.rel)
otu_table(ps.rel) <- otu_table(oturel, taxa_are_rows = T)


ps.rel@otu_table[is.na(ps.rel@otu_table)] <- 0.5
plot(density(ps.rel@otu_table))

pvals <- c()
for (i in 1:1000){
  bla <- kruskal.test(t(ps.rel@otu_table[i,]),ps.rel@sam_data$Treat)
  pvals[i] <- bla$p.value
}

which(pvals<0.05)

i=44
boxplot(t(ps.rel@otu_table[i,])~ps.rel@sam_data$Treat)
boxplot(t(ps.col@otu_table[i,])~ps.col@sam_data$Treat)
boxplot(t(ps.cea@otu_table[i,])~ps.cea@sam_data$Treat)


# 
# #library(vegan3d)
# #ordiplot3d(vare.cca, type = "h")
# 
# data("dune")
# data("dune.env")
# 
# dune.ca <- cca(dune)
# ef <- envfit(dune.ca, dune.env, permutations = 999)
# ef
# 
# dune.ca <- cca(t(ps@otu_table))
# sam_data.f <- data.frame(ps@sam_data[,unlist(lapply(data.frame(ps@sam_data),class))=="factor"])
# sam_data.f$Sample.Type <- NULL
# sam_data.f$Sampling.date <- NULL
# sam_data.f$Plate <- NULL
# sam_data.f$Bacteria_name <- NULL
# sam_data.f$Sample.Type <- NULL
# 
# 
# ef <- envfit(dune.ca, sam_data.f, permutations = 999)
# ef
# 
# plot(dune.ca, display = "sites")
# plot(ef)
# 
# 
# plot(dune.ca, display = "sites", type = "p")
# with(sam_data.f, ordiellipse(dune.ca, Treat, kind = "se", conf = 0.95))
# with(sam_data.f, ordispider(dune.ca, Treat, col = "blue", label= TRUE))
# with(sam_data.f, ordihull(dune.ca, Treat, col="blue", lty=2))
# 
# sel2 <- c("pHcol","VAcol","BAcol","PAcol","LAcol","Shannon","AAcol","Stomach.cont","Treat","Dmile")
# sel2 <- sel2[sel2 %in% colnames(sam_data.nf)]
# 
# select <- c("pHcol", "AAcol", "VAcol")
rda0 <- rda(t(ps@otu_table) ~ 1, sam_data.nf)

ps <- rarefy_even_depth(ps, sample.size = 10000)

rda1 <- rda(t(ps@otu_table) ~ ., sam_data.nf)
rda1 <- rda(t(ps@otu_table) ~ pH + LA, sam_data.f)

rda2 <- rda(t(ps@otu_table), sam_data.nf[,sel2])

plot(rda0)
plot(rda1)
plot(rda2)
text(rda2$CCA$wa[,c(1,2)], ,rownames(rda2$CCA$wa))


rda.res1 <- ordistep(rda0, scope = formula(rda1), perm.max = 200)
rda.res2 <- ordistep(rda1, perm.max = 200) 
rda.res3 <- ordistep(rda0, scope = formula(rda1), direction="forward", perm.max = 200)




## load vegan
require("vegan")

## load the Dune data
data(dune, dune.env)

## PCA of the Dune data
mod <- rda0

## plot the PCA
scl=3
plot(mod, scaling = scl)
colvec <- c("red2", "green4", "mediumblue","magenta")
points(mod, display = "sites", col=colvec[ps@sam_data$Treat], pch=16, scaling = scl)
ordihull(mod, groups=ps@sam_data$Treat, col=colvec, scaling = scl)
ordispider(mod, groups = ps@sam_data$Treat, label = TRUE, col=colvec, scaling = scl)
