setwd("/media/ggj/ggjlab/RData/agingatlas/lifespan/Mouse_4stage_2/Kidney/")

library(Seurat)
library(openxlsx)
library(data.table)
##########load AUcell data
anno_auc <- read.table("/media/ggj/ggjlab/RData/agingatlas/MCA2/all/SCENIC/results/anno_600gene_324810cell_0401.csv",sep=",",header=T,row.names = 1)
auc <- fread("/media/ggj/ggjlab/RData/agingatlas/MCA2/all/SCENIC/results/AUCELL.csv") 
auc <- as.data.frame(auc)  
auc[1:5,1:5]
rownames(auc) <- auc$Cell
auc <- auc[,-1]
auc <- as.data.frame(t(auc))
rownames(anno_auc) <- anno_auc$cellID
anno_auc <- anno_auc[colnames(auc),]

##########load Kidney
load("./Kidney_SCT.RData")

marker <- read.xlsx("./Kidney_SCT_marker.xlsx",1)
marker <- na.omit(marker)
marker$celltype <- paste0(marker$celltype,"_",marker$cluster)
new.cluster.ids<- marker$celltype
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc@meta.data$celltype<-Idents(pbmc)

anno <- pbmc@meta.data
anno <- anno[c("celltype","stage")]
anno$cellID <- rownames(anno)
anno$celltype <- as.character(anno$celltype)
unique(anno$celltype)

#####select muscle
celltype.use <- unique(anno$celltype)
celltype.use
celltype.use <- celltype.use[-c(1,2,3,4,6,8,9,10,11,22,27,29)]
celltype.use

anno.use <- anno[anno$celltype%in%celltype.use,] 

####
anno.final <- anno.use[intersect(rownames(anno.use),anno_auc$cellID),]
#anno.final$celltype <- gsub("Neuron_Crabp1 high|mature neuron|mature Neuron","Neuron",anno.final$celltype)
#anno.final$celltype <- gsub("Radial glial cell|Astrocyte|glial cell|Oligodendrocyte","Glial cell",anno.final$celltype)

unique(anno.final$celltype)

anno.final$type<-paste0(anno.final$celltype,"_",anno.final$stage)
test<-as.data.frame(table(anno.final$type))
use<-test[test$Freq>49,]

auc.use <- auc[,anno.final$cellID]
pbmc<-CreateSeuratObject(auc,meta.data = anno.final)
Idents(pbmc)<-pbmc$type
pbmc.use<-subset(pbmc,idents=as.character(use$Var1))
avg.alltiessue.allstage.50cell<-AverageExpression(pbmc.use)

avg.alltiessue.allstage.50cell<-as.data.frame(avg.alltiessue.allstage.50cell$RNA)
avg.alltiessue.allstage.50cell <- na.omit(avg.alltiessue.allstage.50cell)
cor.spearman<-cor(avg.alltiessue.allstage.50cell,method = "spearman")

min(cor.spearman)
lower<-cor.spearman
lower[lower.tri(lower)]<-NA
library(reshape2)
meltdata<-melt(lower)
meltdata<-na.omit(meltdata)
meltdata.use<-meltdata[-which(meltdata$value==1),]
meltdata.final<-meltdata.use[which(meltdata.use$value>=0.9),]
colnames(meltdata.final)<-c("source","target","value")
#meltdata.final$celltype<-colsplit(meltdata.final$source,"_",names = c("time","celltype"))$celltype
#meltdata.final$time<-colsplit(meltdata.final$source,"_",names = c("time","celltype"))$time
meltdata.final$stage<-colsplit(meltdata.final$source,"_",names = c("cluster","tissue","stage"))$stage
meltdata.final$tissue<-colsplit(meltdata.final$source,"_",names = c("cluster","tissue","stage"))$tissue
meltdata.final$cluster<-colsplit(meltdata.final$source,"_",names = c("cluster","tissue","stage"))$cluster

write.csv(meltdata.final,file = "./Kidney_correlation_50min_0.9.csv",row.names = F,quote = F)
