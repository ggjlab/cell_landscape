setwd("/media/ggj/ggjlab/RData/agingatlas/MCA2/all/SCENIC/results/")
library(reshape2)
library(data.table)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(Hmisc)
####TF auc
anno <- read.table("./anno_600gene_324810cell_0401.csv",sep=",",header=T,row.names = 1)
#anno_accurate <- read.csv("../../anno/Mouse_accurate_annotation_10tissue_20220328.csv",row.names = 1)
#rownames(anno_accurate) <- anno_accurate$cellID
#setdiff(anno$cellID,anno_accurate$cellID)
#anno_accurate <- anno_accurate[intersect(anno$cellID,anno_accurate$cellID),]
anno_accurate <- anno
auc <- fread("./AUCELL.csv") 
auc <- as.data.frame(auc)  
auc[1:5,1:5]
rownames(auc) <- auc$Cell
auc <- auc[,-1]
auc <- as.data.frame(t(auc))
rownames(anno_accurate) <- anno_accurate$cellID
anno <- anno_accurate
anno <- anno[intersect(colnames(auc),anno$cellID),]

setwd("/media/ggj/ggjlab/RData/agingatlas/MCA2/all/")
library(data.table)
library(ggplot2)
library(reshape2)
tissue <- unique(anno$tissue)
tissue <- tissue[-c(15)]

for(i in 1:length(tissue)){
  file <- paste0("./pathway/GO:BP/use_14tissue/",tissue[i],".regulonAUC.top500.RData")
  load(file)
  aucdata1 <- aucdata
  aucdata1 <- as.data.frame(t(aucdata1))
  
  aucdata1$stage <- colsplit(rownames(aucdata1),tissue[i],names = c("n1","n2"))$n1
  unique(aucdata1$stage)
  aucdata1 <- aucdata1[aucdata1$stage%in%c("Adult","OneYear","EighteenMonths","TwoYears"),]
  aucdata1$stage <- factor(aucdata1$stage,levels=c("Adult","OneYear","EighteenMonths","TwoYears"))
  #ggplot(aucdata1,aes(x=stage,y=down))+geom_boxplot()+theme_bw()
  #ggplot(aucdata1,aes(x=stage,y=up))+geom_boxplot()+theme_bw()
  
  anno1 <- anno[anno$tissue==tissue[i],]
  cellname <- intersect(rownames(aucdata1),anno1$cellID)
  #setdiff(anno$cellID,rownames(aucdata1))
  
  aucdata1 <- aucdata1[cellname,]
  auc1 <- auc[,cellname]
  auc1 <- as.data.frame(t(auc1))
  
  use <- cbind(aucdata1,auc1)
  #ggplot(use,aes(x=stage,y=`Esrrg(+)`))+geom_boxplot()
  
  use <- use[,-3]
  use$down <- as.numeric(use$down)
  use$up <- as.numeric(use$up)
  corr <- rcorr(use,method = "spearman")
  corr <- as.data.frame(corr)
  corr$TF <- rownames(corr)

  save <- paste0("./SCENIC/results/geneset_correlation/",tissue[i],"_correlation.csv")
  write.csv(corr,file=save)
  }


#########
setwd("/media/ggj/ggjlab/RData/agingatlas/MCA2/all/SCENIC/results/geneset_correlation/")
library(reshape2)
file <- list.files("./")
file <- file[-10]
tissue <- colsplit(file,"_",names = c("n1","n2"))$n1

down <- NULL
up <- NULL
for(i in 1:length(file)){
  temp <- read.csv(file[i],row.names = 1)
  temp.down <- temp[temp$down>0.3,]
  temp.down <- temp.down['down']
  temp.down$tissue <- tissue[i]
  temp.down <- temp.down[-1,]
  temp.down$TF <- rownames(temp.down)
  temp.down <- na.omit(temp.down)
  down <- rbind(down,temp.down)
  
  temp.up <- temp[temp$up>0.3,]
  temp.up <- temp.up['up']
  temp.up$tissue <- tissue[i]
  temp.up <- temp.up[-1,]
  temp.up$TF <- rownames(temp.up)
  temp.up <- na.omit(temp.up)
  up <- rbind(up,temp.up)
  
}

up$TF <- gsub("\\(\\+\\)","",up$TF)
down$TF <- gsub("\\(\\+\\)","",down$TF)

up.use <- as.data.frame(table(up$TF))
up.use <- up.use[up.use$Freq>5,]
down.use <- as.data.frame(table(down$TF))
down.use <- down.use[down.use$Freq>4,]

down <- down[down$TF%in%as.character(down.use$Var1),]
up <- up[up$TF%in%as.character(up.use$Var1),]

dir.create("./result")
setwd("./result/")
write.csv(down,file="down_TF_cor_0.3_all.csv",quote = F)
write.csv(up,file="up_TF_cor_0.3_all.csv",quote = F)

###anno
downTF <- data.frame(TF=unique(down$TF),type="TF")
rownames(downTF) <- downTF$TF
downTF <- downTF[as.character(down.use$Var1),]
downTF$num <- down.use$Freq

upTF <- data.frame(TF=unique(up$TF),type="TF")
rownames(upTF) <- upTF$TF
upTF <- upTF[as.character(up.use$Var1),]
upTF$num <- up.use$Freq

tissue <- data.frame(TF=down$tissue,type="tissue")


write.csv(downTF,file="down_TF_cor_0.3_anno.csv",quote = F)
write.csv(upTF,file="up_TF_cor_0.3_anno.csv",quote = F)
write.csv(tissue,file="tissue_cor_0.3_anno.csv",quote = F)
#########################
library(ggplot2)
library(pheatmap)
result <- NULL
for(i in 1:length(file)){
  temp <- read.csv(file[i],row.names = 1)
  temp <- temp[,c(1,2)]
  temp <- na.omit(temp)
  temp <- temp[-c(1,2),]
  
  colnames(temp) <- paste0(tissue[i],"_",colnames(temp))
  result <- merge(result,temp,by="row.names",all=T)
  rownames(result)<-result$Row.names
  result<-result[,-1]
  
}

result <- na.omit(result)

result.up <- result[,grep("up",colnames(result))]
result.up[result.up>0.3] <- 0.3
result.up[result.up<(-0.3)] <- c(-0.3)
result.up <- result.up[as.character(up.use[up.use$Freq>1,]$Var1),]
bk <- c(seq(-0.9,-0.1,by=0.03),seq(0,0.3,by=0.01))
bk
pheatmap(t(result.up),color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),
                              colorRampPalette(colors=c("white","red"))(length(bk)/2)),method="ward.D2")
pheatmap(t(result.up))


result.down <- result[,grep("down",colnames(result))]
result.down[result.down>0.3] <- 0.3
result.down[result.down<(-0.3)] <- c(-0.3)
result.down <- result.down[as.character(down.use[down.use$Freq>1,]$Var1),]
bk <- c(seq(-0.9,-0.1,by=0.03),seq(0,0.3,by=0.01))
bk
pheatmap(t(result.down),color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),
                                colorRampPalette(colors=c("white","red"))(length(bk)/2)),method="ward.D2")
pheatmap(t(result.down))
