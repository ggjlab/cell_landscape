library(reshape2)
library(Matrix)
library(readr)
library(tidyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggdendro)
library(grid)
library(dendextend)
library(tidyverse)
library(RColorBrewer)
library(ape)
library(phylogram)
library(circlize)
message("Tree plot")
aurocs<-read.csv("./AUROC.csv",row.names = 1,check.names = F)
hh<-length(aurocs[1,])/2
aurocs<-aurocs[1:hh,(hh+1):(hh*2)]
aurocs<-as.matrix(aurocs)
rownames(aurocs) <- gsub("Sample1[|]","",rownames(aurocs))
colnames(aurocs) <- gsub("Sample2[|]","",colnames(aurocs))
rownames(aurocs) <- gsub("\\D","",rownames(aurocs))
colnames(aurocs) <- gsub("\\D","",colnames(aurocs))
aurocs <- aurocs[,as.character(rownames(aurocs))]
total.dist<-as.dist(1-aurocs)
total.tree<-hclust(total.dist)
total.tree<-as.dendrogram(total.tree)
#---------add color-------------------------------------------------------------
Clu <- read.table("397Celltype.txt",sep="\t")
Clu <- Clu[Clu$Species == "Mouse",]
Clu$Celltype <- gsub("\\D","",Clu$Celltype)
rownames(Clu) <- Clu$Celltype
Clu <- Clu[as.character(rownames(aurocs)),]
Cellcluster <- Clu$Cellcluster
Cellcluster<-as.factor(Cellcluster)

color_regions = c("#E6AB02", "#66A61E", "#D95F02", "#1B9E77", "#E7298A", "#E31A1C", "#A6761D", "#DCDCDC", "#FFFF99", "#7570B3", "#FF7F00", "#A65628", "#B3CDE3", "#BC80BD", "#A6CEE3","#984EA3", "#CCEBC5","#E41A1C","#4DAF4A","#BEBADA","#B3DE69","#CAB2D6","#FFFFB3","#33A02C","#B15928", "#6A3D9A","#FBB4AE","blue","#FB8072","#FFFF33","#CCEBC5","#A6761D","#2c7fb8","#fa9fb5","#BEBADA")
names(color_regions) = c("Secretory" ,"Muscle" ,"Neuron" ,"Immune" , "Epithelial","Glia","Proliferating","Other","Parenchymal","Stromal","Phagocytes","Pharynx","Rectum","Coelomocytes","Intestine","Hepatocyte","Germline","Endothelial","Erythroid","Testis","Unknown","Midgut","Hemocytes" ,"Hindgut","Embryo","Fat","SalivaryGland","Gastrodermis","DigFilaments","Pigment","BasementMembrane","Endoderm","Mesenchyme","FatBody","Female")
color_regions_use = color_regions[as.character(levels(Cellcluster))]
which(is.na(color_regions_use))
Clu.order <- Clu[order.dendrogram(total.tree),]
celltype_color <- color_regions_use[as.character(Clu.order$Cellcluster)]
dend<-total.tree %>%
  color_branches(k=11) %>%
  set("labels_cex",0.6)%>%
  set("leaves_pch",19)%>%
  set("leaves_col",celltype_color) %>% 
  set("labels_colors",celltype_color)

circlize_dendrogram(dend,dend_track_height = 0.8)#,facing = "inside")







