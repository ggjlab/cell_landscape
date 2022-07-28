setwd("/media/ggj/ggjlab/RData/agingatlas/MCA2/all/SCENIC/results/")

library(ggplot2)
library(pheatmap)
library(reshape2)
library(scales)
library(RColorBrewer)
library(data.table)
library(Seurat)

anno <- read.table("./anno_600gene_324810cell_0401.csv",sep=",",header=T,row.names = 1)
auc <- fread("./AUCELL.csv") 
auc <- as.data.frame(auc)  
auc[1:5,1:5]
rownames(auc) <- auc$Cell
auc <- auc[,-1]
auc <- as.data.frame(t(auc))
rownames(anno) <- anno$cellID
anno <- anno[colnames(auc),]
pbmc <- CreateSeuratObject(auc)
setdiff(colnames(auc),rownames(anno))
table(anno$celllineage)
#anno <- anno[colnames(auc),]
anno$celllineage <- gsub("Stem|Adipose|other","Other",anno$celllineage)

Idents(pbmc) <- paste0(anno$stage,"_",anno$celllineage)
avg <- AverageExpression(pbmc)
avg <- as.data.frame(avg$RNA)


###fuzzy means
##cmeans
library(e1071)
set.seed(1234)
cl<-cmeans(t(avg),15,1000,verbose=F,method="cmeans",m=1.25)
#cl<-cmeans(t(avg),15,1000,verbose=F,method="cmeans",m=1.3)
centermatrix <- cl$centers
membership <- cl$membership
centermatrix <- t(centermatrix)
col <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(10))
anno.col <- list(
  lineage=c(
    Immune="#379B7B",
    Epithelial="#CF622B",
    Stromal="#6E70AC",
    Germ="#D7B777",
    Endothelial="#6BA143",
    Erythroid="#DAA033",
    Muscle="#9E7536",
    Neuron="#77c286",
    Other="#AFA3CB",
    Secretory="#FF8C00",
    Proliferating="#A0522D"
    #Hemocyte="#416EA7"
  ),
  
  stage=c(
    `E10.5` = col[1],
    `E12.5` = col[2],
    `Fetal` = col[3],
    `Neonatal` = col[4],
    `TenDays` = col[5],
    `ThreeWeeks` = col[6],
    `Adult` = col[7],
    `OneYear` = col[8],
    `EighteenMonths` = col[9],
    `TwoYears` = col[10]
  )
  
)


type <- colsplit(rownames(centermatrix),"_",names = c("n1","n2"))
type$type <- rownames(centermatrix)
type <- type[order(type$n1),]
type <- type[order(type$n2),]

type.use <- NULL
for(i in 1:length(unique(type$n2))){
  temp <- type[type$n2== unique(type$n2)[i],]
  # if(i==3){
  #   temp <- temp[c(2,3,5,6,7,8,1,4,9),]
  # }
  if(i%in%c(4,10)){
    temp <- temp[c(3,4,6,7,1,5,2,8),]
  }
  if(i==9){
    temp <- temp[c(2,4,5,7,8,1,6,3,9),]
  }
  if(i%in%c(1,2,3,5,6,7,8,11)){
    temp <- temp[c(2,3,5,6,8,9,1,7,4,10),]
  }
  type.use <- rbind(type.use,temp)
}
type.use <- na.omit(type.use)
annotation <- data.frame(
  lineage=type.use$n2,
  stage=type.use$n1
)

rownames(annotation) <- type.use$type

zz <- centermatrix[type.use$type,]
#zz <- scale(zz)
max(zz);min(zz)
zz[zz>1.5] <- 1.5
zz[zz<(-1.5)] <- (-1.5)
zz <- apply(zz, 1, rescale,to=c(-2,2))
pdf("max3_20cluster.pdf",width = 10,height = 10)
pheatmap::pheatmap(t(zz),
                   #cellwidth = 9,cellheight = 9,
                   border=F,annotation_row = annotation,annotation_colors = anno.col,
                   clustering_method = "ward.D2",cluster_rows = F,cluster_cols = T,show_rownames = T,
                   color = colorRampPalette(c("#4575B4" ,"#91BFDB" ,"#E0F3F8", "white" ,"#FEE090" ,"#FC8D59" ,"#D73027"))(100))
dev.off()

type.use$num <- type.use$n2
type.use$num <- gsub("Germ","1",type.use$num)
type.use$num <- gsub("Erythroid","2",type.use$num)
type.use$num <- gsub("Immune","3",type.use$num)
type.use$num <- gsub("Endothelial","4",type.use$num)
type.use$num <- gsub("Neuron","5",type.use$num)
type.use$num <- gsub("Muscle","6",type.use$num)
type.use$num <- gsub("Stromal","7",type.use$num)
type.use$num <- gsub("Epithelial","8",type.use$num)
type.use$num <- gsub("Secretory","9",type.use$num)
type.use$num <- gsub("Proliferating","10",type.use$num)
type.use$num <- gsub("Other","10",type.use$num)

type.use$num <- as.numeric(type.use$num)
type.use <- type.use[order(type.use$num),]

annotation <- data.frame(
  lineage=type.use$n2,
  stage=type.use$n1
)
rownames(annotation) <- type.use$type

zz <- centermatrix[type.use$type,]

max(zz);min(zz)
#zz[zz<0.8] <- 0
#zz[zz>0.8] <- 1
zz[zz>2] <- 2
zz[zz<(-2)] <- (-2)
library(scales)
zz <- apply(zz, 1, rescale,to=c(-2,2))

pheatmap::pheatmap(t(zz),
                   #cellwidth = 9,cellheight = 9,
                   border=F,annotation_row = annotation,annotation_colors = anno.col,
                   clustering_method = "ward.D2",cluster_rows = F,cluster_cols = T,
                   color = colorRampPalette(c("#4575B4" ,"#91BFDB" ,"#E0F3F8", "white" ,"#FEE090" ,"#FC8D59" ,"#D73027"))(100))

zz <- t(zz)
zz <- zz[,c(7,9,6,8,3,10,4,15,13,11,12,1,5,2,14)]

#pdf("./figure/heatmap_cmeans_15module_rescale2_0227.pdf",width = 10,height = 10)
pheatmap::pheatmap(zz,
                   #cellwidth = 9,cellheight = 9,
                   border=F,annotation_row = annotation,annotation_colors = anno.col,
                   clustering_method = "ward.D2",cluster_rows = F,cluster_cols = F,
                   gaps_col = c(11,14),
                   color = colorRampPalette(c("#4575B4" ,"#91BFDB" ,"#E0F3F8", "white" ,"#FEE090" ,"#FC8D59" ,"#D73027"))(100))
#dev.off()

class(colnames(membership))
membership.use <- membership[,c(7,9,6,8,3,10,4,15,13,11,12,1,5,2,14)]
colnames(membership.use) <- c(1:15)
centermatrix.use <- centermatrix[,c(7,9,6,8,3,10,4,15,13,11,12,1,5,2,14)]
colnames(centermatrix.use) <- c(1:15)

zz <- centermatrix.use[type.use$type,]

max(zz);min(zz)
#zz[zz<0.8] <- 0
#zz[zz>0.8] <- 1
zz[zz>2] <- 2
zz[zz<(-2)] <- (-2)
zz <- apply(zz, 1, rescale,to=c(-2,2))
pdf("./figure/heatmap_cmeans_15module_rescale2_600gene_0409.pdf",width = 10,height = 10)
pheatmap::pheatmap(t(zz),
                   #cellwidth = 9,cellheight = 9,
                   border=F,annotation_row = annotation,annotation_colors = anno.col,
                   clustering_method = "ward.D2",cluster_rows = F,cluster_cols = F,
                   gaps_col = c(9,13),#gaps_row = c(5,11,18,25,30,37,44,51,58,62,69),
                   color = colorRampPalette(c("#4575B4" ,"#91BFDB" ,"#E0F3F8", "white" ,"#FEE090" ,"#FC8D59" ,"#D73027"))(100))
dev.off()


###0.2
summaryTF <- NULL
for(i in 1:length(colnames(membership.use))){
  temp <- as.data.frame(membership.use[,i])
  colnames(temp) <- "value"
  temp$TF <- rownames(temp)
  temp <- as.data.frame(temp[order(temp$value,decreasing = T),])
  
  temp <- temp[temp$value>0.2,]
  temp$module <- i
  
  summaryTF <- rbind(summaryTF,temp)
}
table(summaryTF$module)
summaryTF$module <- as.numeric(summaryTF$module)
write.csv(summaryTF,file="./summaryTF_0.2_0409.csv")
#rm(pbmc)
save.image("./Aucell_result_0409.RData")


