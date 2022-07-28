setwd("/media/ggj/ggjlab/RData/agingatlas/MCA2/all/MAST_time_related_gene/")
library(reshape2)
library(openxlsx)
library(clusterProfiler)
library(org.Mm.eg.db)
file <- list.files("./")
file <- file[-c(10,13)]
upgene <- NULL
downgene <- NULL
for(i in 1:length(file)){
  load(file[i])
  
  temp.up <- de_res[de_res$age.logFC_z>0&de_res$age.H_fdr<1e-6,]
  temp.up$tissue <- colsplit(file[i],"_",names = c("n1","n2"))$n1
  upgene <- rbind(upgene,temp.up)
  
  temp.down <- de_res[de_res$age.logFC_z<0&de_res$age.H_fdr<1e-6,]
  temp.down$tissue <- colsplit(file[i],"_",names = c("n1","n2"))$n1
  downgene <- rbind(downgene,temp.down)
  gc()
}

upgene.summary <- as.data.frame(table(upgene$gene))
downgene.summary <- as.data.frame(table(downgene$gene))
save(upgene,downgene,file="./result/result_14tissue_1e6.RData")

###1e3
file <- file[-c(10,13)]
upgene <- NULL
downgene <- NULL
for(i in 1:length(file)){
  load(file[i])
  
  temp.up <- de_res[de_res$age.logFC_z>0&de_res$age.H_fdr<1e-3,]
  temp.up$tissue <- colsplit(file[i],"_",names = c("n1","n2"))$n1
  upgene <- rbind(upgene,temp.up)
  
  temp.down <- de_res[de_res$age.logFC_z<0&de_res$age.H_fdr<1e-3,]
  temp.down$tissue <- colsplit(file[i],"_",names = c("n1","n2"))$n1
  downgene <- rbind(downgene,temp.down)
  gc()
}

upgene.summary <- as.data.frame(table(upgene$gene))
downgene.summary <- as.data.frame(table(downgene$gene))
save(upgene,downgene,file="./result/result_14tissue_1e3.RData")
###1e10
upgene <- NULL
downgene <- NULL
file
file <- file[-c(10,13)]
for(i in 1:length(file)){
  load(file[i])
  
  temp.up <- de_res[de_res$age.logFC_z>0&de_res$age.H_fdr<1e-10,]
  temp.up$tissue <- colsplit(file[i],"_",names = c("n1","n2"))$n1
  upgene <- rbind(upgene,temp.up)
  
  temp.down <- de_res[de_res$age.logFC_z<0&de_res$age.H_fdr<1e-10,]
  temp.down$tissue <- colsplit(file[i],"_",names = c("n1","n2"))$n1
  downgene <- rbind(downgene,temp.down)
  gc()
}

upgene.summary <- as.data.frame(table(upgene$gene))
downgene.summary <- as.data.frame(table(downgene$gene))
save(upgene,downgene,file="./result/result_14tissue_1e10.RData")
write.csv(downgene,file="./result/result_14tissue_1e10_downgene.csv")
write.csv(upgene,file="./result/result_14tissue_1e10_upgene.csv")
downgene$type <- "Down"
upgene$type <- "Up"
gene <- rbind(downgene,upgene)
write.csv(gene,file="./result/result_14tissue_1e10_gene.csv")
# ###
# dir.create("./result")
# num=6
# symbol=as.character(upgene.summary[upgene.summary$Freq>num,]$Var1)
# write.table(symbol,file="./result/upgene_7tissues_1e10_from14tissues.txt",sep="\t",quote = F,row.names = F,col.names = F)
# symbol=as.character(downgene.summary[downgene.summary$Freq>num,]$Var1)
# write.table(symbol,file="./result/downgene_7tissues_1e10_from14tissue.txt",sep="\t",quote = F,row.names = F,col.names = F)
# eg = bitr(symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
# id = as.character(eg[,2])
# head(id)
# #length(gene)
# length(id)
# ego <- enrichGO(gene = id,
#                 OrgDb = org.Mm.eg.db,
#                 ont = "BP",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff = 0.05,
#                 qvalueCutoff = 0.05,
#                 readable = TRUE)
# View(ego@result)

# up_GO <- ego
# View(down_GO@result)

# load("./result/result_1e10.RData")
# upgene.summary <- as.data.frame(table(upgene$gene))
# downgene.summary <- as.data.frame(table(downgene$gene))
# num=6
# symbol=as.character(upgene.summary[upgene.summary$Freq>num,]$Var1)
# upgene.use <- upgene[upgene$gene%in%symbol,]
# upgene.use <- upgene.use[,c(1,6)]
# upgene.use2 <- dcast(upgene.use,gene~tissue)

#############
# essay.gene.down <- read.xlsx("../essay/elife-62293-supp2 (1).xlsx",1)
# 
# essay.gene.down[essay.gene.down=="TRUE"] <- 1
# essay.gene.down[essay.gene.down=="FALSE"] <- 0
# essay.gene.down.process <- melt(essay.gene.down)
# essay.gene.down.process <- essay.gene.down.process[essay.gene.down.process$value==1,]
# essay.down.summary <- as.data.frame(table(essay.gene.down.process$X1))
