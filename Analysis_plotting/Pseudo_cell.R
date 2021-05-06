library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(DoubletFinder)
library(RColorBrewer)
library(ggthemes)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

load("./Adult-Lung_rm.batch.Rdata")
Stage <- "AdultLung"

message("Pseudocell start")
raw_UMI <- mergedge.rm
Inter<-raw_UMI
raw1<-Inter
raw.10e6<-t(t(raw1)/colSums(raw1))*1000000
Inter<-raw.10e6
Inter<-as.data.frame(Inter)

id<-mergedge.rm.Anno
Inter.id<-cbind(as.character(rownames(id)),as.character(id[,1]))
Inter.id<-as.data.frame(Inter.id)
colnames(Inter.id)<-cbind("Cell_id","Celltype")
rownames(Inter.id)<-Inter.id$Cell_id
Inter.id$Stage <- Stage
Inter.id$Celltype <- as.factor(paste(Stage, Inter.id$Celltype,sep = ""))

pseudocell.size = 20 ## 10 test
new_ids_list = list()
for (i in 1:length(levels(Inter.id$Celltype))) {
  cluster_id = levels(Inter.id$Celltype)[i]
  cluster_cells <- rownames(Inter.id[Inter.id$Celltype == cluster_id,])
  cluster_size <- length(cluster_cells)		
  pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
  pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
  names(pseudo_ids) <- sample(cluster_cells)	
  new_ids_list[[i]] <- pseudo_ids		
}

new_ids <- unlist(new_ids_list)
new_ids <- as.data.frame(new_ids)
new_ids_length <- table(new_ids)

new_colnames <- rownames(new_ids)
all.data<-Inter[,as.character(new_colnames)] ###add
all.data <- t(all.data)###add

new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
                    list(name=new_ids[,1]),FUN=mean)
rownames(new.data)<-new.data$name
new.data<-new.data[,-1]

new_ids_length<-as.matrix(new_ids_length)##
short<-which(new_ids_length<10)##
new_good_ids<-as.matrix(new_ids_length)
if(length(short) > 0 ) {
  new_good_ids<-as.matrix(new_ids_length[-short,])#
}
#
result<-t(new.data)[,rownames(new_good_ids)]
rownames(result)<-rownames(Inter)
cc<-gsub("[_]Cell.*$","",colnames(result))
new.phe<-cbind(colnames(result),Stage,cc) ##
colnames(new.phe)<-c("Sample_ID","Study_ID","Celltype")
OUT <- paste(Stage,"_pse20.Rdata",sep = "")
Pse_Matrix <- result
Pse_pheno<- new.phe
save(Pse_Matrix,Pse_pheno,file=OUT)










