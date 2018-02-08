ps <- readRDS("ps.S16404.rds")
facet1="Sample.Type"
facet2="Treat"
#common for facet facet1
sample <- "piglet"
plotname=paste("taxaplot","otu","png",sep=".")

#library("phyloseq")
library("ggplot2")
#library(gplots)

colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00");
taxonrank <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

#maximimum number of colours
N <- length(colours)
#N=20
temp <- otu_table(ps)
abund_table <- t(temp)
x<-abund_table/rowSums(abund_table)
x<-x[,order(colSums(x),decreasing=TRUE)]
taxa_list<-colnames(x)[1:N]
N<-length(taxa_list)
if (length(colnames(x))-length(taxa_list) != 1){
  new_x<-data.frame(x[,colnames(x) %in% taxa_list],Low_abundance=rowSums(x[,!colnames(x) %in% taxa_list]))
} 
if (length(colnames(x))-length(taxa_list) == 1){
  new_x<-as.data.frame(x)
}

df<-NULL
for (i in 1:dim(new_x)[2]){
  tmp<-data.frame(row.names=NULL,
                  Sample=unname(ps@sam_data[,sample]),
                  Taxa=rep(colnames(new_x)[i],
                           dim(new_x)[1]),
                  Value=new_x[,i],
                  facet1=unname(ps@sam_data[,eval(facet1)]), 
                  facet2=unname(ps@sam_data[,eval(facet2)]))
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
zotu=N+1
df$Taxa <- factor(df$Taxa,levels=rev(levels(df$Taxa)))
for (Otu in grep("Low_abundance",levels(df$Taxa),invert = T, value = T)){
  zotu=zotu-1
  levels(df$Taxa)[levels(df$Taxa)==Otu] <- paste("OTU",zotu,unname(ps@tax_table[Otu,6]))
}

p <- ggplot(df,aes(Sample,Value,fill=Taxa))
p <- p + geom_bar(stat="identity")
p <- p + facet_grid(facet1 ~ facet2, drop=TRUE,scale="free",space="free_x")
p <- p + scale_fill_manual(values=rev(colours[0:(N+1)]))
p <- p + theme_bw()+ylab("Proportions")
p <- p + scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.spacing = unit(0.3, "lines"))
p <- p + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p <- p + labs(title = "ZOTU")
p <- p + guides(fill=guide_legend(ncol=1, reverse=F))
p <- p + theme(panel.spacing = unit(1, "lines"))
p <- p + labs(x = sample)

png(width = 1200, height = 1200, file=plotname)
p = p + theme(text = element_text(size=20), axis.text.x = element_text(angle=90, hjust=1)) 
print(p)
dev.off()
