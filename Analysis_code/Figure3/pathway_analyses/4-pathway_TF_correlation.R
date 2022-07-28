setwd("/media/ggj/ggjlab/RData/agingatlas/MCA2/all/SCENIC/results/")
library(reshape2)
library(data.table)
library(ggplot2)
library(ggrepel)
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

tissue <- c("Brain","Liver","Lung","Kidney","Heart","Intestine","Testis","Uterus",
            "Pancreas","Stomach","Bladder","Spleen","Thymus","Prostate")
term <- c("GO:0006119")
dir.create("pathway_correlation")
for(i in 1:length(tissue)){
  file <- paste0("../../pathway/GO:BP/use_10tissue/",tissue[i],".regulonAUC.top500.RData")
  load(file)
  #aucdata[1:5,1:5]
  aucdata <- aucdata[term,]
  aucdata <- as.data.frame(t(aucdata))
  aucdata$cellID <- rownames(aucdata)
  anno.temp <- anno[anno$tissue==tissue[i]&anno$stage%in%c("Adult","EighteenMonths","TwoYears","OneYear"),]
  
  #aa <- as.data.frame(intersect(anno.temp$cellID,aucdata$cellID))
  
  aucdata.use <- as.data.frame(aucdata[intersect(anno.temp$cellID,aucdata$cellID),])
  aucdata.use <- aucdata.use[term]
  
  auc.temp <- as.data.frame(t(auc[,anno.temp$cellID])) 
  auc.temp2 <- auc.temp[intersect(anno.temp$cellID,rownames(aucdata.use)),]
  auc.temp2 <- auc.temp2[,colSums(auc.temp2)>0]
  auc.temp.use <- cbind(aucdata.use,auc.temp2)
  correlation <- cor(auc.temp.use,method = "spearman")
  
  
  
  data <- as.data.frame(correlation)[term]
  data$TF <- rownames(data)
  colnames(data)[1] <- "correlation"
  data <- data[order(data$correlation,decreasing = T),]
  data <- data[-1,]
  data$order <- rep(c(1:length(rownames(data))))
  
  data$label <- data$TF
  #data$label <- ifelse(abs(data$correlation)>0.6,data$label,"")
  data$label <- ifelse(data$order%in%c(1:10,c(length(rownames(data))-9):length(rownames(data)))|data$TF=="Pparg(+)",data$label,"")
  ggplot(data,aes(y=correlation,x=order))+geom_point()+theme_bw()+xlab("rank")+
    geom_text_repel(data = data, aes(label = label),
                    size = 5,box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.1, "lines"), 
                    segment.color = "black", 
                    show.legend = FALSE,max.overlaps=10000)
  
  savefile <- paste0("./pathway_correlation/",term,"_",tissue[i],".pdf")
  ggsave(savefile,width = 10,height = 6)
  
  file2 <- paste0("./pathway_correlation/",term,"_",tissue[i],".txt")
  write.table(data,file=file2,sep="\t")
}


# lipid metabolic ------------------------------------------------------------
tissue <- c("Brain","Liver","Lung","Kidney","Heart","Intestine","Testis","Uterus",
            "Pancreas","Stomach")
term <- c("GO:0006629")
dir.create("./pathway_correlation/celltype")
for(i in 1:length(tissue)){
  file <- paste0("../../pathway/GO:BP/use_10tissue/",tissue[i],".regulonAUC.top500.RData")
  load(file)
  #aucdata[1:5,1:5]
  aucdata <- aucdata[term,]
  aucdata <- as.data.frame(t(aucdata))
  aucdata$cellID <- rownames(aucdata)
  anno.temp2 <- anno[anno$tissue==tissue[i],]
  
  celltype <- unique(anno.temp2$celllineage)
  for(j in 1:length(celltype)){
    anno.temp <- anno.temp2[anno.temp2$celllineage==celltype[j],]
    aucdata.use <- as.data.frame(aucdata[intersect(anno.temp$cellID,aucdata$cellID),])
    aucdata.use <- aucdata.use[term]
    
    auc.temp <- as.data.frame(t(auc[,anno.temp$cellID])) 
    auc.temp2 <- auc.temp[intersect(anno.temp$cellID,rownames(aucdata.use)),]
    auc.temp.use <- cbind(aucdata.use,auc.temp2)
    correlation <- cor(auc.temp.use,method = "spearman")
    
    
    
    data <- as.data.frame(correlation)[term]
    data$TF <- rownames(data)
    colnames(data)[1] <- "correlation"
    data <- data[order(data$correlation,decreasing = T),]
    data <- data[-1,]
    data$order <- rep(c(1:length(rownames(data))))
    
    data$label <- data$TF
    #data$label <- ifelse(abs(data$correlation)>0.6,data$label,"")
    data$label <- ifelse(data$order%in%c(1:10,408:417)|data$TF=="Pparg(+)",data$label,"")
    ggplot(data,aes(y=correlation,x=order))+geom_point()+theme_bw()+xlab("rank")+
      geom_text_repel(data = data, aes(label = label),
                      size = 5,box.padding = unit(0.5, "lines"),
                      point.padding = unit(0.1, "lines"), 
                      segment.color = "black", 
                      show.legend = FALSE,max.overlaps=10000)+ggtitle(paste0(tissue[i],"_",celltype[j]))
    
    
    savefile <- paste0("./pathway_correlation/celltype/",term,"_",tissue[i],"_",celltype[j],".pdf")
    ggsave(savefile,width = 10,height = 6)
    
    file2 <- paste0("./pathway_correlation/celltype/",term,"_",tissue[i],"_",celltype[j],".txt")
    write.table(data,file=file2,sep="\t")
    
  }
  #aa <- as.data.frame(intersect(anno.temp$cellID,aucdata$cellID))
  
  
  
}


# oxidative ---------------------------------------------------------------

tissue <- c("Brain","Liver","Lung","Kidney","Heart","Intestine","Testis","Uterus",
            "Pancreas","Stomach")
term <- c("path:mmu00190")
dir.create("pathway_correlation")
for(i in 1:length(tissue)){
  file <- paste0("../../pathway/KEGG/use_10tissue/",tissue[i],".regulonAUC.top500.RData")
  load(file)
  #aucdata[1:5,1:5]
  aucdata <- aucdata[term,]
  aucdata <- as.data.frame(t(aucdata))
  aucdata$cellID <- rownames(aucdata)
  anno.temp <- anno[anno$tissue==tissue[i],]
  
  #aa <- as.data.frame(intersect(anno.temp$cellID,aucdata$cellID))
  
  aucdata.use <- as.data.frame(aucdata[intersect(anno.temp$cellID,aucdata$cellID),])
  aucdata.use <- aucdata.use[term]
  
  auc.temp <- as.data.frame(t(auc[,anno.temp$cellID])) 
  auc.temp2 <- auc.temp[intersect(anno.temp$cellID,rownames(aucdata.use)),]
  auc.temp.use <- cbind(aucdata.use,auc.temp2)
  correlation <- cor(auc.temp.use,method = "spearman")
  
  
  
  data <- as.data.frame(correlation)[term]
  data$TF <- rownames(data)
  colnames(data)[1] <- "correlation"
  data <- data[order(data$correlation,decreasing = T),]
  data <- data[-1,]
  data$order <- rep(c(1:length(rownames(data))))
  
  data$label <- data$TF
  #data$label <- ifelse(abs(data$correlation)>0.6,data$label,"")
  data$label <- ifelse(data$order%in%c(1:10,408:417)|data$TF=="Pparg(+)",data$label,"")
  ggplot(data,aes(y=correlation,x=order))+geom_point()+theme_bw()+xlab("rank")+
    geom_text_repel(data = data, aes(label = label),
                    size = 5,box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.1, "lines"), 
                    segment.color = "black", 
                    show.legend = FALSE,max.overlaps=10000)
  
  savefile <- paste0("./pathway_correlation/",term,"_",tissue[i],".pdf")
  ggsave(savefile,width = 10,height = 6)
  
  file2 <- paste0("./pathway_correlation/",term,"_",tissue[i],".txt")
  write.table(data,file=file2,sep="\t")
}