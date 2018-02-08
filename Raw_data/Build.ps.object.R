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