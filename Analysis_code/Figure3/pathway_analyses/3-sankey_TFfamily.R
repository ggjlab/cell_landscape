setwd("/media/ggj/ggjlab/RData/agingatlas/MCA2/all/SCENIC/results/geneset_correlation/result/")

library(reshape2)
library(networkD3)
library(stringr)
library(RColorBrewer)
up <- read.csv("./up_TF_cor_0.3.csv",row.names = 1)
length(unique(up$TF))
TF_info <- read.table("/media/ggj/ggjlab/RData/agingatlas/cross_species/TF_information/Mus_musculus_TF.txt",
                      sep="\t",header=T)

use <- TF_info[TF_info$Symbol%in%up$TF,]
table(use$Family)
length(unique(use$Family))
up <- up[,c(3,2,1)]
# development -------------------------------------------------------------


colnames(up)[3] <- "value"
HT_OUT2 <- up
HT_OUT2 <- merge(HT_OUT2,use,by.x="TF",by.y="Symbol")
HT_OUT2 <- HT_OUT2[,-c(4,5,7,8)]
HT_OUT2 <- HT_OUT2[,c(4,2,3)]
colnames(HT_OUT2)[1] <- "TF"
HT_OUT2$Cellcluster1 <- HT_OUT2$TF
HT_OUT2$Cellcluster2 <- HT_OUT2$tissue

edges<-cbind(as.character(HT_OUT2$Cellcluster1),as.character(HT_OUT2$Cellcluster2),HT_OUT2$value)
#edges<-cbind(as.character(HT_OUT2$Celltype1),as.character(HT_OUT2$Celltype2),HT_OUT$Mean_AUROC)
edges<-as.data.frame(edges)
colnames(edges)<-c("N1","N2","Value")
edges$N1 = as.character(edges$N1)    
edges$N2 = as.character(edges$N2)  
d3links <- edges
d3nodes <- data.frame(name = unique(c(edges$N1, edges$N2)), stringsAsFactors = FALSE)
d3nodes$seq <- 0:(nrow(d3nodes) - 1)

d3links <- merge(d3links, d3nodes, by.x="N1", by.y="name")
names(d3links)[4] <- "source"
d3links <- merge(d3links, d3nodes, by.x="N2", by.y="name")
names(d3links)[5] <- "target"
names(d3links)[3] <- "value"
#d3links <- merge(d3links,use,by.x='N1',by.y='Symbol')
d3links <- subset(d3links, select=c("source", "target", "value"))
d3nodes <- subset(d3nodes, select=c("name"))

tissue <- data.frame(name=unique(HT_OUT2$tissue),Family="tissue")
rownames(use) <- use$Symbol
use <- use[use$Family%in%as.character(d3nodes$name)[1:13],]
d3nodes <- use[,c(2,4)]
colnames(d3nodes) <- c("Family","name")
d3nodes$Family <- d3nodes$name
d3nodes <- d3nodes[!duplicated(d3nodes),]
d3nodes <- d3nodes[,c(2,1)]
d3nodes <- rbind(d3nodes,tissue)

col <-  colorRampPalette(brewer.pal(9,"YlOrRd"))(length(unique(d3nodes$Family))-1)
col <- data.frame(Family=unique(d3nodes$Family),color=c(col,"#A9A9A9"))
col <- merge(col,d3nodes,by="Family")
family <- c(col$name)
col1 <- c(col$color)
#write.table(family,file="./sankey_up_family.txt",row.names = F,col.names = F)
#write.table(col1,file="./sankey_up_color.txt",row.names = F,col.names = F)
color <- 'd3.scaleOrdinal().domain(["bHLH","CSD","E2F","HMG","IRF","MYB","NF-YB","Nrf1","Others","P53","TF_bZIP","THAP", "zf-C2H2", "tissue" ])
.range(["#FFFFCC", "#FFF3AE", "#FEE692", "#FED976" ,"#FEBF5A" ,"#FDA546" ,"#FD8D3C", "#FC6330" ,"#F33C25", "#E31A1C"
,"#C90822" ,"#A80026", "#800026","#A9A9A9"])'

#color setting
net<-sankeyNetwork(Links = d3links,
                   Nodes = d3nodes, 
                   Source = "source",
                   Target = "target", 
                   Value = "value", 
                   NodeID = "name",
                   NodeGroup = "Family",
                   colourScale = JS(color),
                   #LinkGroup = 'Family',
                   units = "votes",
                   height =1000,
                   width =800,
                   nodePadding = 10,
                   sinksRight = FALSE, 
                   fontSize = 15, 
                   nodeWidth = 5
)
net
#sankeyNetworkOutput()
write.csv(d3nodes,file="./sankey_up_anno.csv",quote = F,row.names = F)

col.plot <- col[-which(col$Family=="tissue"),]


ggplot(col.plot,aes(x=Family,y=name,fill=Family))+geom_bar(stat='identity')+
  scale_fill_manual(values=c("#FFFFCC", "#FFF3AE", "#FEE692", "#FED976" ,"#FEBF5A" ,"#FDA546" ,"#FD8D3C", "#FC6330" ,"#F33C25", "#E31A1C"
                             ,"#C90822" ,"#A80026", "#800026"))
ggsave("./sankey_up_nocolor_family_color_anno.pdf",width = 8,height = 6)
################down

down <- read.csv("./down_TF_cor_0.3.csv",row.names = 1)
length(unique(down$TF))
TF_info <- read.table("/media/ggj/ggjlab/RData/agingatlas/cross_species/TF_information/Mus_musculus_TF.txt",
                      sep="\t",header=T)

use <- TF_info[TF_info$Symbol%in%down$TF,]
table(use$Family)
length(unique(use$Family))
down <- down[,c(3,2,1)]


colnames(down)[3] <- "value"
HT_OUT2 <- down
HT_OUT2 <- merge(HT_OUT2,use,by.x="TF",by.y="Symbol")
HT_OUT2 <- HT_OUT2[,-c(4,5,7,8)]
HT_OUT2 <- HT_OUT2[,c(4,2,3)]
colnames(HT_OUT2)[1] <- "TF"
HT_OUT2$Cellcluster1 <- HT_OUT2$TF
HT_OUT2$Cellcluster2 <- HT_OUT2$tissue

edges<-cbind(as.character(HT_OUT2$Cellcluster1),as.character(HT_OUT2$Cellcluster2),HT_OUT2$value)
#edges<-cbind(as.character(HT_OUT2$Celltype1),as.character(HT_OUT2$Celltype2),HT_OUT$Mean_AUROC)
edges<-as.data.frame(edges)
colnames(edges)<-c("N1","N2","Value")
edges$N1 = as.character(edges$N1)    
edges$N2 = as.character(edges$N2)  
d3links <- edges
d3nodes <- data.frame(name = unique(c(edges$N1, edges$N2)), stringsAsFactors = FALSE)
d3nodes$seq <- 0:(nrow(d3nodes) - 1)

d3links <- merge(d3links, d3nodes, by.x="N1", by.y="name")
names(d3links)[4] <- "source"
d3links <- merge(d3links, d3nodes, by.x="N2", by.y="name")
names(d3links)[5] <- "target"
names(d3links)[3] <- "value"
#d3links <- merge(d3links,use,by.x='N1',by.y='Symbol')
d3links <- subset(d3links, select=c("source", "target", "value"))
d3nodes <- subset(d3nodes, select=c("name"))

tissue <- data.frame(name=unique(HT_OUT2$tissue),Family="tissue")
which(use$Symbol=="Jund")
use <- use[-28,]
rownames(use) <- use$Symbol
use <- use[use$Family%in%as.character(d3nodes$name)[1:59],]
d3nodes <- use[,c(2,4)]
colnames(d3nodes) <- c("Family","name")
d3nodes$Family <- d3nodes$name
d3nodes <- d3nodes[!duplicated(d3nodes),]
d3nodes <- d3nodes[,c(2,1)]
d3nodes <- rbind(d3nodes,tissue)

col <-  colorRampPalette(brewer.pal(9,"YlOrRd"))(length(unique(d3nodes$Family))-1)
col <- data.frame(Family=unique(d3nodes$Family),color=c(col,"#A9A9A9"))
col <- merge(col,d3nodes,by="Family")
family <- c(col$name)
col1 <- c(col$color)
#write.table(family,file="./sankey_down_family.txt",row.names = F,col.names = F)
#write.table(col1,file="./sankey_down_color.txt",row.names = F,col.names = F)

color <- 'd3.scaleOrdinal().domain(["bHLH","CP2","CTF/NFI","ESR-like","HMG" ,"Homeobox","Nrf1","RXR-like" ,"SRF"   ,  
 "STAT","TEA","TF_bZIP","THR-like","zf-C2H2" ,"tissue" ])
.range(["#FFFFD9", "#F5FBC3" ,"#EAF6B1" ,"#D5EEB2" ,"#BBE4B5", "#95D5B8" ,"#70C7BD", "#4FBBC1" ,"#35AAC2", "#2296C0"
,"#1E7DB6", "#2161A9", "#174697" ,"#0C2C84","#A9A9A9"])'

#color setting
net<-sankeyNetwork(Links = d3links,
                   Nodes = d3nodes, 
                   Source = "source",
                   Target = "target", 
                   Value = "value", 
                   NodeID = "name",
                   NodeGroup = "Family",
                   colourScale = JS(color),
                   #LinkGroup = 'Family',
                   units = "votes",
                   height =1500,
                   width =800,
                   nodePadding = 10,
                   sinksRight = FALSE, 
                   fontSize = 15, 
                   nodeWidth = 5
)
net
write.csv(d3nodes,file="./sankey_down_anno.csv",quote = F,row.names = F)

col.plot <- col[-which(col$Family=="tissue"),]


ggplot(col.plot,aes(x=Family,y=name,fill=Family))+geom_bar(stat='identity')+
  scale_fill_manual(values=c("#FFFFD9", "#F5FBC3" ,"#EAF6B1" ,"#D5EEB2" ,"#BBE4B5", "#95D5B8" ,"#70C7BD", "#4FBBC1" ,"#35AAC2", "#2296C0"
                             ,"#1E7DB6", "#2161A9", "#174697" ,"#0C2C84"))
ggsave("./sankey_dwon_nocolor_family_color_anno.pdf",width = 8,height = 6)