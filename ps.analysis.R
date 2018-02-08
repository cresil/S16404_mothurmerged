#### init ####

setwd("C:/Users/mdavi/Documents/Trouw/S16404_mothurmerged")

library(phyloseq)
library(ggplot2)
library(sas7bdat)

source("/Users/mdavi/Documents/R_functions/taxa_facet_barplot_asv.R")
source("/Users/mdavi/Documents/R_functions/calculate_adonis.R")
source("/Users/mdavi/Documents/R_functions/calculate_rarefactioncurves.R")

ps <- phyloseq::import_mothur(mothur_shared_file = "processed.agc.unique_list.shared",
                               mothur_constaxonomy_file = "processed.agc.unique_list.0.03.cons.taxonomy")

samdata_colon <- read.csv("Sample_data_libary_gatc3.txt", sep="\t", row.names = 1)
samdata_colon <- samdata_colon[gsub("_.*","",colnames(ps@otu_table))[1:40],]

samdata_ceacum <- data.frame(animal.number=paste("11",sprintf("%02d",as.numeric(gsub(".*_","",colnames(ps@otu_table))[41:80])),sep=""))
samdata_ceacum$Sample.Type <- as.factor(rep("Ceacum",40))

samdata <-data.frame(Sample.Type=c(as.character(samdata_colon$Sample.Type), as.character(samdata_ceacum$Sample.Type)),
           animal.number=c(as.character(samdata_colon$animal.number), as.character(samdata_ceacum$animal.number))
)

rownames(samdata) <- colnames(ps@otu_table)
sample_data(ps) <- sample_data(samdata)

sampleid.piglet <- read.csv("sample_piglet.csv",sep=";")
rownames(sampleid.piglet) <- sampleid.piglet[,1]

ps@sam_data$piglet <- as.factor(sampleid.piglet[as.vector(ps@sam_data$animal.number),2])

Acids <- read.csv("S16404.01_SCFA_NH4_pH_DM_Final.txt",sep="\t")
rownames(Acids) <- Acids$piglet
Acids <- Acids[as.vector(ps@sam_data$piglet),]

Necropsy <- read.csv("Necropsy_measurements.txt",sep="\t")
rownames(Necropsy) <- Necropsy$Piglet
Necropsy <- Necropsy[as.vector(ps@sam_data$piglet),]

sample_data(ps) <- sample_data(cbind(ps@sam_data,Acids,Necropsy))
ps@sam_data$Treat <- as.factor(ps@sam_data$Treat)

sasdata <- sas7bdat::read.sas7bdat(file = "datamark.sas7bdat")
sasdata2 <- sas7bdat::read.sas7bdat(file="extramark.sas7bdat")
colnames(sasdata2)[colnames(sasdata2)=="Pen"] <- "pen"
colnames(sasdata2)[colnames(sasdata2)=="Piglet"] <- "Bignr"

sasdata <- sasdata[,colnames(sasdata) %in% colnames(sasdata2)]
sasdata2 <- sasdata2[,colnames(sasdata2) %in% colnames(sasdata)]
sasdata <- rbind(sasdata,sasdata2)
rownames(sasdata)<-sasdata$Bignr

sample_data(ps) <- sample_data(cbind(cbind(ps@sam_data,sasdata[as.character(ps@sam_data$piglet),])))

ps@sam_data$Bignr <-NULL
ps@sam_data$sow <- NULL
ps@sam_data$Piglet <- NULL
ps@sam_data$Sample.id <- NULL

ps@sam_data$Sample_id <- as.factor(ps@sam_data$Sample_id)
ps@sam_data$Day <- as.factor(ps@sam_data$Day)
ps@sam_data$pen <- as.factor(ps@sam_data$pen)


ps@sam_data$ES.Treitz..mg.ml <- as.numeric(ps@sam_data$ES.Treitz..mg.ml)
ps@sam_data$ES.ileum..mg.ml <- as.numeric(ps@sam_data$ES.ileum..mg.ml)
ps@sam_data$Full.weight.LI..kg <- as.numeric(ps@sam_data$Full.weight.LI..kg)
ps@sam_data$Empty.weight.LI..kg <- as.numeric(ps@sam_data$Empty.weight.LI..kg)

write.csv(cbind(rownames(ps@sam_data),paste0(ps@sam_data$Sample.Type,ps@sam_data$Treat)),file="stamp.sample.txt")

#### modify acid profiles

ps@sam_data$AA <- c(ps@sam_data$AAcol[ps@sam_data$Sample.Type=="Feaces"],ps@sam_data$Aacae[ps@sam_data$Sample.Type!="Feaces"])
ps@sam_data$pH <- c(ps@sam_data$pHcol[ps@sam_data$Sample.Type=="Feaces"],ps@sam_data$pHcae[ps@sam_data$Sample.Type!="Feaces"])
ps@sam_data$VA <- c(ps@sam_data$VAcol[ps@sam_data$Sample.Type=="Feaces"],ps@sam_data$VAcae[ps@sam_data$Sample.Type!="Feaces"])
ps@sam_data$PA <- c(ps@sam_data$PAcol[ps@sam_data$Sample.Type=="Feaces"],ps@sam_data$PAcae[ps@sam_data$Sample.Type!="Feaces"])
ps@sam_data$BA <- c(ps@sam_data$BAcol[ps@sam_data$Sample.Type=="Feaces"],ps@sam_data$Bacae[ps@sam_data$Sample.Type!="Feaces"])
ps@sam_data$LA <- c(ps@sam_data$LAcol[ps@sam_data$Sample.Type=="Feaces"],ps@sam_data$Lacae[ps@sam_data$Sample.Type!="Feaces"])

ps@sam_data$dm <- c(ps@sam_data$Dmcol[ps@sam_data$Sample.Type=="Feaces"],ps@sam_data$DMcae[ps@sam_data$Sample.Type!="Feaces"])
ps@sam_data$IBU <- c(ps@sam_data$IBUcol[ps@sam_data$Sample.Type=="Feaces"],ps@sam_data$IBUcae[ps@sam_data$Sample.Type!="Feaces"])
ps@sam_data$IVA <- c(ps@sam_data$IVAcol[ps@sam_data$Sample.Type=="Feaces"],ps@sam_data$IVAcae[ps@sam_data$Sample.Type!="Feaces"])

ps <- prune_taxa(names(!which(rowSums(ps@otu_table!=0)>1)),ps)
ps@sam_data$Sow <- as.factor(ps@sam_data$Sow)

sample_data(ps) <- sample_data(ps@sam_data[,!duplicated(colnames(ps@sam_data))])
lapply(ps@sam_data, class)

ps@sam_data$Day <- as.numeric(ps@sam_data$Day) + 22
saveRDS(ps, file = "ps.S16404.rds")

#### plot to check data ####

ps <- readRDS("ps.S16404.rds")
p <- plot_richness(ps.rare)
p + facet_grid(variable ~ Treat, drop=TRUE,scale="free",space="free_x")


p <- plot_ordination(physeq = ps.rare, ordination = ordinate(physeq = ps.rare, method = "PCoA", distance = "bray"), color="pH", shape="Sample.Type", label="Treat")
#p <- 
p <- p + geom_line(aes(group=animal.number))
p + scale_color_gradient(low="red", high="green")

ps.rare@sam_data$Sow <- as.factor(ps.rare@sam_data$Sow)

p <- plot_ordination(physeq = ps.rare, ordination = ordinate(physeq = ps.rare, method = "PCoA", distance = "bray"), color="Sow", shape="Sample.Type", label="Treat")
#p <- 
p <- p + geom_line(aes(group=animal.number))
p


p <- plot_ordination(physeq = ps, ordination = ordinate(physeq = ps, method = "PCoA", distance = "bray"), color="Sow", shape="Sample.Type", label="Treat")
#p <- 
p <- p + geom_line(aes(group=Sow))
p

p <- plot_ordination(physeq = ps, ordination = ordinate(physeq = ps, method = "PCoA", distance = "bray"), color="Day", shape="Sample.Type", label="Treat")
#p <- 
p <- p + geom_line(aes(group=Day))
p



ord.ps.rare.pcoa.bray=ordinate(physeq = ps.rare, method = "PCoA", distance = "bray" )
var="BiW"
p.bray <- plot_ordination(ps.rare, ord.ps.rare.pcoa.bray, color = var)
#p.bray <- p.bray + geom_line(aes(group = animal.number, color=Treat))
plot(p.bray)

p.jac <- plot_ordination(ps.rare, ordinate(physeq = ps.rare, method = "PCoA", distance = "jaccard", binary=T ), color = "Sample.Type")
p.jac <- p.jac + geom_line(aes(group = animal.number, color=Treat))
plot(p.jac)

#### cca/rda ####

ps <- readRDS("ps.S16404.rds")
ps <- rarefy_even_depth(ps)
ps <- prune_taxa(names(which(rowSums(!ps@otu_table==0)>1)),ps)
library(vegan)
#library(MASS)
#library(vegan3d)

bla <- lapply(data.frame(ps@sam_data),class)
sam_data.nf <- data.frame(ps@sam_data[,!unlist(lapply(data.frame(ps@sam_data),class))=="factor"])
sam_data.f <- data.frame(ps@sam_data[,unlist(bla)=="factor"])

vare.cca <- cca(t(ps@otu_table) ~ pH + AA + BA + PA + VA + IBU + IVA , sam_data.nf)
vare.cca <- cca(t(ps@otu_table) ~ . , sam_data.nf)

plot(vare.cca)
colvec <- rainbow(4)
colvec <- colvec[ps@sam_data[rownames(vare.cca$CCA$wa)]$Treat]

points(vare.cca$CCA$wa[,1:2],cex=2, pch=18, col=colvec, )
legend(x="topright",legend=c(1:4),col=colvec)


vare.cca <- cca(t(ps@otu_table) ~ depart + pHcol + Sample.Type, cbind(sam_data.f, sam_data.nf))
plot(vare.cca)
colvec <- rainbow(4)
colvec <- colvec[ps@sam_data[rownames(vare.cca$CCA$wa)]$Treat]

points(vare.cca$CCA$wa[,1:2],cex=2, pch=18, col=colvec, )
legend(x="topright",legend=c(1:4),col=colvec)


rda1 <- rda(t(ps@otu_table) ~ AX + CELL, sam_data.f)
rda1 <- rda(t(ps@otu_table) ~ pH + piglet, data.frame(ps@sam_data))

mod <- rda1

## plot the PCA
scl=3
plot(mod, scaling = scl)
colvec <- c("red2", "green4", "mediumblue","magenta")
points(mod, display = "sites", col=colvec[ps@sam_data$Treat], pch=16, scaling = scl)
ordihull(mod, groups=ps@sam_data$Treat, col=colvec, scaling = scl)
ordispider(mod, groups = ps@sam_data$Treat, label = TRUE, col=colvec, scaling = scl)


#### rel change table ####

ps.rel <- prune_taxa(names(which(rowSums(ps@otu_table)>50)),ps)
ps.rel <- ps

ps.colon <- prune_samples(ps@sam_data$Sample.Type=="Feaces", ps.rel)
ps.ceacum <- prune_samples(ps@sam_data$Sample.Type!="Feaces", ps.rel)

otu_tab.colon <- ps.colon@otu_table
otu_tab.ceacum <- ps.ceacum@otu_table[,match(as.character(ps.ceacum@sam_data$animal.number),as.character(ps.colon@sam_data$animal.number))]

bla <- otu_tab.ceacum/(otu_tab.colon+otu_tab.ceacum)
bla[is.na(bla)] <- 0.5

source("/Users/mdavi/Documents/R_functions/taxa_facet_barplot_asv.R")
source("/Users/mdavi/Documents/R_functions/calculate_adonis.R")
source("/Users/mdavi/Documents/R_functions/calculate_rarefactioncurves.R")

### deseq2 AX ####

library("DESeq2")
packageVersion("DESeq2")

colnames(ps@tax_table) <- c("Kingdom","Phylum","Class","Order","Family","Genus")

ps.ceacum <- prune_samples(ps@sam_data$Sample.Type=="Ceacum", ps)
ps.ceacum <- prune_samples(ps@sam_data$Sample.Type=="Feaces", ps)

ps.ceacum.2 <- prune_samples(ps.ceacum@sam_data$Treat==1 | ps.ceacum@sam_data$Treat==2, ps.ceacum)
ps.ceacum.3 <- prune_samples(ps.ceacum@sam_data$Treat==1 | ps.ceacum@sam_data$Treat==3, ps.ceacum)
ps.ceacum.4 <- prune_samples(ps.ceacum@sam_data$Treat==1 | ps.ceacum@sam_data$Treat==4, ps.ceacum)

ps.ceacum.23 <- prune_samples(ps.ceacum@sam_data$Treat==2 | ps.ceacum@sam_data$Treat==3, ps.ceacum)

diagdds = phyloseq_to_deseq2(ps.ceacum.4, ~ Treat)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.1
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.ceacum.4)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(aes(size=baseMean)) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#### ordistep ####

library(vegan)

data(dune)
data(dune.env)

ps <- readRDS("ps.S16404.rds")
ps.rare <- rarefy_even_depth(ps)

otu <- as.matrix(ps@otu_table)
df <- as.data.frame(unclass(ps@sam_data))

df <- df[,sapply(df, function(x) length(levels(x)))!=1]

df

mod1 <- rda(t(otu) ~ ., df)
plot(mod1)
ordistep(mod1)
ordistep(mod2, scope = formula(mod2))

mod2 <- rda(t(otu) ~ AA, df)
plot(mod2)
mod3 <- rda(t(otu) ~ AA + pH, df)

plot(mod3)
points(mod3$CCA$wa[,1:2],cex=2, pch=18, col=colvec)
legend(x="topright",legend=c(1:4),col=colvec)

plot(mod3,type = c("text"))

mod4 <- rda(t(otu) ~ Sample.Type + pH, df)
plot(mod4,type = c("text"))

ordiellipse(ord, Management, col=1:4, kind = "ehull", lwd=3)
ordiellipse(ord, Management, col=1:4, draw="polygon")
ordispider(ord, Management, col=1:4, label = TRUE)


plot_ordination(physeq = )

colvec <- rainbow(4)
colvec <- colvec[ps@sam_data[rownames(mod3$CCA$wa)]$Treat]

points(mod3$CCA$wa[,1:2],cex=2, pch=18, col=colvec, )
legend(x="topright",legend=c(1:4),col=colvec)




plot(mod2, type = c("text"))
plot(mod3, type = c("text"))
plot(mod4, type = c("text"))

data(dune)
mod <- rda(dune, scale = TRUE)
biplot(mod2, scaling = "symmetric")

## different type for species and site scores
biplot(mod1, scaling = "symmetric", type = c("text", "points"))
biplot(mod1, scaling = "symmetric", type = c("points", "text"))

biplot(mod2, scaling = "symmetric", type = c("text", "points"))
biplot(mod2, scaling = "symmetric", type = c("points", "text"))

biplot(mod3, scaling = "symmetric", type = c("text", "points"))
biplot(mod3, scaling = "symmetric", type = c("points", "text"))

biplot(mod4, scaling = "symmetric", type = c("text", "points"))
biplot(mod4, scaling = "symmetric", type = c("points", "text"))

GP <- ps

ordcap = ordinate(GP, "CAP", "bray", ~pH)
plot_ordination(ps, ordcap, "samples", color="Treat")

p <- plot_ordination(ps, ordcap, "samples", color="AA") + ggplot2::scale_color_gradient(low="green", high="red")
p
p <- plot_ordination(ps, ordcap, "samples", color="Treat")
p

RDA = ordinate(GP, "RDA", "bray", ~pH)
plot_ordination(ps, RDA, "samples", color="Treat")


p <- plot_ordination(ps, RDA, "samples", color="AA") + ggplot2::scale_color_gradient(low="green", high="red")
p
p <- plot_ordination(ps, RDA, "samples", color="Treat")
p

#### phyloseq RDA ####

ord <- ordinate(physeq = ps, formula = ~ pH + Sample.Type, method = "RDA", distance = "bray")
# <- ordinate(physeq = ps, formula = . ~ pH + Sample.Type, method = "RDA", distance = "bray")
#ord <- ordinate(physeq = ps, formula = ~pH, method = "RDA", distance = "bray")
ord <- ordinate(physeq = ps, formula = . ~ pH, method = "RDA", distance = "bray")


#ord <- ordinate(physeq = ps, formula = . ~ pH, method = "RDA", distance = "bray")
plot_ordination(physeq = ps, ordination = ord, axes = c(3,4), color = "AA", shape="Treat") + ggplot2::scale_color_gradient(low="red", high="green")
plot_ordination(physeq = ps, ordination = ord, axes = c(3,4), color = "AA", shape="Treat", type="biplot", label="Genus") + ggplot2::scale_color_gradient(low="red", high="green")
plot_ordination(physeq = ps, ordination = ord, axes = c(1,2), color = "AA", shape="Treat", type="biplot", label="Genus") + ggplot2::scale_color_gradient(low="red", high="green")
plot_ordination(physeq = ps, ordination = ord, axes = c(3,4), color = "AA", shape="Treat", type="biplot", label="Genus") + ggplot2::scale_color_gradient(low="red", high="green")


plot_ordination(physeq = ps, ordination = ord, axes = c(3,4), color = "Treat", shape="Sample.Type", label="piglet") + geom_line(aes(group=animal.number))
plot_ordination(physeq = ps, ordination = ord, axes = c(3,4), color = "Treat", shape="Sample.Type", type="biplot", label="Genus") + geom_line(aes(group=animal.number))
plot_ordination(physeq = ps, ordination = ord, axes = c(3,4), color = "Treat", shape="Sample.Type", type="biplot", label="animal.number")


#### metadata correlations #####

library(Maaslin)
InputTSV <- system.file('extdata','maaslin_demo2.tsv', package="Maaslin")
InputConfig <-system.file('extdata','maaslin_demo2.read.config', package="Maaslin")
Maaslin(InputTSV,'maaslin_example_output',strInputConfig=InputConfig)


ps.ce <- prune_samples(ps@sam_data$Sample.Type!="Ceacum", ps)
taxa_names(ps.ce) <- make.unique(ps.ce@tax_table[,6])


library(gplots)

cormat <- sapply(ps.ce@sam_data[,unlist(lapply(ps.ce@sam_data, FUN=is.numeric))],function(x){cor(x, t(ps.ce@otu_table[taxa_sums(ps.ce)>100,]),method = "spearman")})
cormat[is.na(cormat)] <- 0
rownames(cormat) <- make.unique(ps.ce@tax_table[,6])[taxa_sums(ps.ce)>100]

heatmap.2(t(cormatred), trace="none", margins = c(10,15), col=redgreen(30))

cormatred <- cormat[rowSums(cormat<(0.3*-1))>0 & rowSums(cormat>0.3)>0,]
cormatred <- cormatred[,colSums(cormatred<(0.5*-1))>0 & colSums(cormatred>0.5)>0]
#cormatred <- cormat

cormatred

#### pls Ceacum ####

ps <- readRDS("ps.S16404.rds")
ps <- rarefy_even_depth(ps)
ps <- prune_samples(ps@sam_data$Sample.Type=="Ceacum", ps)
ps <- prune_taxa(names(which(rowSums(!ps@otu_table==0)>2)),ps)
taxa_names(ps) <- make.unique(ps@tax_table[,6])

library(mixOmics)

Y <- t(as.matrix(unclass(as.matrix(ps@otu_table))))
Y <- Y[,apply(X = Y, MARGIN = 2,FUN = sd)!=0]
X <- as.matrix(ps@sam_data[,unlist(lapply(ps@sam_data, FUN=is.numeric))])
X <- X[,apply(X = X, MARGIN = 2,FUN = sd)!=0]
X <- X[,grep("col",colnames(X), value = T, invert = T)]
X <- X[,grep("cae",colnames(X), value = T, invert = T)]
X <- X[,grep("ile",colnames(X), value = T, invert = T)]
X <- X[,grep("BiW.1",colnames(X), value = T, invert = T)]
X <- X[,grep("piglet.1",colnames(X), value = T, invert = T)]
X <- X[,grep("pH",colnames(X), value = T, invert = T)]

pca.otu <- pca(Y, ncomp = 10, center = TRUE, scale = TRUE)
pca.otu
plot(pca.otu)

pca.meta <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
pca.meta
plot(pca.meta)

plotIndiv(pca.otu, comp = c(1, 2), group = ps@sam_data$Treat,
          legend = TRUE, title = 'OTU table')

plotIndiv(pca.meta, comp = c(1, 2), group = ps@sam_data$Treat,
#          ind.names = as.character(ps.MRNA.overlap.T2@sam_data$Subject),
          legend = TRUE, title = 'Meta_data')

liver.pls <- pls(Y, X, ncomp = 10, mode = "regression")
liver.spls <- spls(Y, X, ncomp =10, keepX = c(20,20,20), keepY= c(20,20,20), mode = "regression")

tune.pls <- perf(liver.pls, validation = "Mfold", folds = 10, progressBar = FALSE, nrepeat = 10)
plot(tune.pls$Q2.total)
abline(h = 0.0975)


tune.spls <- perf(liver.spls, validation = "Mfold", folds = 10, progressBar = FALSE, nrepeat = 10)
plot(tune.spls$Q2.total)
abline(h = 0.0975)

plotVar(liver.spls, comp =1:2, 
        var.names = list(X.label = colnames(Y), 
                         Y.label = TRUE), cex = c(4, 5))


cim(liver.spls, comp = 1:3, xlab = "co-variables", ylab = "", 
    margins = c(10, 10))

#### pls Colon ####

ps <- readRDS("ps.S16404.rds")
ps <- rarefy_even_depth(ps)
ps <- prune_taxa(names(which(rowSums(!ps@otu_table==0)>1)),ps)

ps <- prune_samples(ps@sam_data$Sample.Type!="Ceacum", ps)

library(mixOmics)

Y <- t(as.matrix(unclass(as.matrix(ps@otu_table))))
Y <- Y[,apply(X = Y, MARGIN = 2,FUN = sd)!=0]
X <- as.matrix(ps@sam_data[,unlist(lapply(ps@sam_data, FUN=is.numeric))])
X <- X[,apply(X = X, MARGIN = 2,FUN = sd)!=0]

pca.otu <- pca(Y, ncomp = 10, center = TRUE, scale = TRUE)
pca.otu
plot(pca.otu)

pca.meta <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
pca.meta
plot(pca.meta)

plotIndiv(pca.otu, comp = c(1, 2), group = ps@sam_data$Treat,
          legend = TRUE, title = 'OTU table')

plotIndiv(pca.meta, comp = c(1, 2), group = ps@sam_data$Treat,
          #          ind.names = as.character(ps.MRNA.overlap.T2@sam_data$Subject),
          legend = TRUE, title = 'Meta_data')

liver.pls <- pls(Y, X, ncomp = 10, mode = "regression")
liver.spls <- spls(Y, X, ncomp =10, keepX = c(50,50,50), keepY= c(50,50,50), mode = "regression")

tune.pls <- perf(liver.pls, validation = "Mfold", folds = 10, progressBar = FALSE, nrepeat = 10)
plot(tune.pls$Q2.total)
abline(h = 0.0975)


tune.spls <- perf(liver.spls, validation = "Mfold", folds = 10, progressBar = FALSE, nrepeat = 10)
plot(tune.spls$Q2.total)
abline(h = 0.0975)

plotVar(liver.spls, comp =1:2, 
        var.names = list(X.label = colnames(Y), 
                         Y.label = TRUE), cex = c(4, 5))


cim(liver.spls, comp = 1:3, xlab = "clinic", ylab = "genes", 
    margins = c(7, 7))

#### pls Colon Genus ####

ps <- readRDS("ps.S16404.rds")
ps <- rarefy_even_depth(ps)
ps <- prune_taxa(names(which(rowSums(!ps@otu_table==0)>1)),ps)

ps <- prune_samples(ps@sam_data$Sample.Type!="Ceacum", ps)
ps <- tax_glom(physeq =  ps, taxrank = "Genus", NArm = F)

library(mixOmics)

Y <- t(as.matrix(unclass(as.matrix(ps@otu_table))))
Y <- Y[,apply(X = Y, MARGIN = 2,FUN = sd)!=0]
X <- as.matrix(ps@sam_data[,unlist(lapply(ps@sam_data, FUN=is.numeric))])
X <- X[,apply(X = X, MARGIN = 2,FUN = sd)!=0]

pca.otu <- pca(Y, ncomp = 10, center = TRUE, scale = TRUE)
pca.otu
plot(pca.otu)

pca.meta <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
pca.meta
plot(pca.meta)

plotIndiv(pca.otu, comp = c(1, 2), group = ps@sam_data$Treat,
          legend = TRUE, title = 'OTU table')

plotIndiv(pca.meta, comp = c(1, 2), group = ps@sam_data$Treat,
          #          ind.names = as.character(ps.MRNA.overlap.T2@sam_data$Subject),
          legend = TRUE, title = 'Meta_data')

liver.pls <- pls(Y, X, ncomp = 10, mode = "regression")
liver.spls <- spls(Y, X, ncomp =10, keepX = c(50,50,50), keepY= c(50,50,50), mode = "regression")

tune.pls <- perf(liver.pls, validation = "Mfold", folds = 10, progressBar = FALSE, nrepeat = 10)
plot(tune.pls$Q2.total)
abline(h = 0.0975)


tune.spls <- perf(liver.spls, validation = "Mfold", folds = 10, progressBar = FALSE, nrepeat = 10)
plot(tune.spls$Q2.total)
abline(h = 0.0975)

plotVar(liver.spls, comp =1:2, 
        var.names = list(X.label = colnames(Y), 
                         Y.label = TRUE), cex = c(4, 5))


cim(liver.spls, comp = 1:3, xlab = "clinic", ylab = "genes", 
    margins = c(7, 7))


#### pls Ceacum Genus ####

ps <- readRDS("ps.S16404.rds")
#ps@sam_data$GF <- ((ps@sam_data$BW.d23.24.kg/ps@sam_data$BiW)^(1/(as.numeric(as.character(ps@sam_data$Day))+22))-1)*100

ps <- rarefy_even_depth(ps)
ps <- prune_samples(ps@sam_data$Sample.Type=="Ceacum", ps)
ps <- tax_glom(physeq =  ps, taxrank = "Rank6", NArm = F)
ps <- prune_taxa(names(which(rowSums(!ps@otu_table==0)>2)),ps)
taxa_names(ps) <- make.unique(ps@tax_table[,6])

library(mixOmics)

Y <- t(as.matrix(unclass(as.matrix(ps@otu_table))))
Y <- Y[,apply(X = Y, MARGIN = 2,FUN = sd)!=0]
X <- as.matrix(ps@sam_data[,unlist(lapply(ps@sam_data, FUN=is.numeric))])
X <- X[,apply(X = X, MARGIN = 2,FUN = sd)!=0]
X <- X[,grep("col",colnames(X), value = T, invert = T)]
X <- X[,grep("cae",colnames(X), value = T, invert = T)]
X <- X[,grep("ile",colnames(X), value = T, invert = T)]
X <- X[,grep("BiW.1",colnames(X), value = T, invert = T)]
X <- X[,grep("piglet.1",colnames(X), value = T, invert = T)]
X <- X[,grep("pH",colnames(X), value = T, invert = T)]

pca.otu <- pca(Y, ncomp = 10, center = TRUE, scale = TRUE)
pca.otu
plot(pca.otu)

pca.meta <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
pca.meta
plot(pca.meta)

plotIndiv(pca.otu, comp = c(1, 2), group = ps@sam_data$Treat,
          legend = TRUE, title = 'OTU table')

plotIndiv(pca.meta, comp = c(1, 2), group = ps@sam_data$Treat,
          #          ind.names = as.character(ps.MRNA.overlap.T2@sam_data$Subject),
          legend = TRUE, title = 'Meta_data')

liver.pls <- pls(Y, X, ncomp = 10, mode = "regression")
liver.spls <- spls(Y, X, ncomp =10, keepX = c(19,19,19), keepY= c(19,19,19), mode = "regression")

tune.pls <- perf(liver.pls, validation = "Mfold", folds = 10, progressBar = FALSE, nrepeat = 10)
plot(tune.pls$Q2.total)
abline(h = 0.0975)

tune.spls <- perf(liver.spls, validation = "Mfold", folds = 10, progressBar = FALSE, nrepeat = 10)
plot(tune.spls$Q2.total)
abline(h = 0.0975)

plotVar(liver.spls, comp =1:2, 
        var.names = list(X.label = colnames(Y), 
                         Y.label = TRUE), cex = c(4, 5))

cim(liver.spls, comp = 1:3, xlab = "clinic", ylab = "genes", 
    margins = c(7, 7))

#### adonis model Colon ####

library(vegan)

ps <- readRDS("ps.S16404.rds")
ps <- rarefy_even_depth(ps)
ps <- prune_taxa(names(which(rowSums(!ps@otu_table==0)>1)),ps)

ps <- prune_samples(ps@sam_data$Sample.Type=="Ceacum", ps)

ps@sam_data$Sow <- as.factor(ps@sam_data$Sow)

df = as(sample_data(ps), "data.frame")
d = phyloseq::distance(ps, "bray")
adonis(d ~ AX*CELL, df)
adonis(d ~ Day+Sow, df)
adonis(d ~ AX*CELL, df, strata=ps@sam_data$Day)
adonis(d ~ AX*CELL, df, strata=ps@sam_data$Sow)
adonis(d ~ AX*CELL, df)



betadisper(d = d,group = ps@sam_data$Treat)

data(varespec)

## Bray-Curtis distances between samples
dis <- vegdist(varespec)

## First 16 sites grazed, remaining 8 sites ungrazed
groups <- factor(c(rep(1,16), rep(2,8)), labels = c("grazed","ungrazed"))

## Calculate multivariate dispersions
mod <- betadisper(d, ps@sam_data$Sow)
mod

## Perform test
anova(mod)

## Permutation test for F
permutest(mod, pairwise = TRUE)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod)

## Draw a boxplot of the distances to centroid for each group
boxplot(mod)
groups[c(2,20)] <- NA
dis[c(2, 20)] <- NA
mod2 <- betadisper(dis, groups) ## warnings
mod2
permutest(mod, control = permControl(nperm = 100))
anova(mod2)
plot(mod2)
boxplot(mod2)
plot(TukeyHSD(mod2))





function(ps){
  require('vegan')
  d = phyloseq::distance(ps, "bray")
  df = as(sample_data(ps), "data.frame")
  pvals=rep("NA",length(colnames(ps@sam_data)))
  names(pvals) <- colnames(ps@sam_data)
  for (i in 1:length(colnames(ps@sam_data))){
    tryCatch({
      test.out <- adonis(d ~ get(colnames(ps@sam_data)[i]), df)
      pvals[colnames(ps@sam_data)[i]] <- test.out$aov.tab$`Pr(>F)`[1]
    }, error=function(e){})
  }
  names(pvals) <- colnames(ps@sam_data)
  sort(pvals,decreasing = T)
}


p <- plot_ordination(physeq = ps, ordination = ordinate(physeq = ps, method = "PCoA", distance = "bray"), color="Sow", shape="Sample.Type", label="Treat")
#p <- 
p <- p + geom_line(aes(group=Sow))
p

p <- plot_ordination(physeq = ps, ordination = ordinate(physeq = ps, method = "PCoA", distance = "bray"), color="Day", shape="Sample.Type", label="Treat")
#p <- 
p <- p + geom_line(aes(group=Day))
p

p <- plot_ordination(physeq = ps, ordination = ordinate(physeq = ps, method = "PCoA", distance = "bray"), color="pen", shape="Sample.Type", label="Treat")
#p <- 
p <- p + geom_line(aes(group=pen))
p


#### rel growth ####
assume exp growth first 3 weeks?

as.numeric(as.character(ps@sam_data$Day))+22

1/ps@sam_data$BW.d23.24.kg/ps@sam_data$BiW