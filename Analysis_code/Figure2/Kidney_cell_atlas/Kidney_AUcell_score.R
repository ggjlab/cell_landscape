setwd("/media/ggj/ggjlab/RData/agingatlas/lifespan/Mouse_4stage_2/Kidney/")

library(Seurat)
library(ggplot2)
library(reshape2)
library(future)
load("./Kidney_SCT.RData")

#grep Epithelial
anno <- read.csv("./Kidney_anno.csv",row.names = 1)
Epi_anno <- anno[anno$lineage=="Epithelial",]
Epi <- pbmc[,Epi_anno$cellID]
DimPlot(Epi)
DimPlot(Epi,split.by = "stage")
table(Epi$stage)
stage <- as.character(unique(Epi$stage))
Idents(Epi) <- "stage"
Epi.use <- subset(Epi,idents = stage[3:10])

DimPlot(Epi.use)

###
library(future)
plan("multicore", workers = 48)
stage.marker <- FindAllMarkers(object =Epi.use, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold  = 0.25)
write.csv(stage.marker,file="./stage_DEG.csv")

####
stage.marker <- read.csv("./stage_DEG.csv",row.names = 1)
stage.marker.use <- stage.marker[stage.marker$p_val_adj<0.05,]
table(stage.marker.use$cluster,stage.marker.use$avg_log2FC>1)

###GO
library(clusterProfiler)
library(org.Mm.eg.db)
final <- NULL
for(i in 1:length(unique(stage.marker.use$cluster))){
  temp <-stage.marker.use[stage.marker.use$cluster==unique(stage.marker.use$cluster)[i],]
  
  gene <- temp[temp$avg_log2FC>=1&temp$p_val_adj<0.05,]$gene
  
  eg = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  id = as.character(eg[,2])
  head(id)
  length(gene)
  length(id)
  
  ##GO
  ego <- enrichGO(gene = id,
                  OrgDb = org.Mm.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  use<- ego@result
  use$type <- unique(stage.marker.use$cluster)[i]
  
  final <- rbind(final,use)
}
final <- final[final$p.adjust<0.05,]
write.csv(final,"GO_result.csv")


####select GO term
GO_term <- c("GO:0008380",
             "GO:0006397",
             "GO:0022613",
             "GO:0042254",
             "GO:0006520",
             "GO:0001822",
             "GO:0072001",
             "GO:0007346",
             "GO:0061440",
             "GO:0055064",
             "GO:0055083",
             "GO:0044282",
             "GO:0006066",
             "GO:0006631",
             "GO:0006790",
             "GO:0010038",
             "GO:0006959",
             "GO:0055078",
             "GO:0001912",
             "GO:0001910",
             "GO:0002705",
             "GO:0001906",
             "GO:0072044",
             "GO:0097205",
             "GO:0032816",
             "GO:0003014",
             "GO:0019882",
             "GO:0001822",
             "GO:0072017",
             "GO:0072070",
             "GO:0003096"
)
library(GSEABase)
library(data.table)
library(AUCell)
geneSets <- getGmt("/media/ggj/ggjlab/RData/agingatlas/MCA2/all/pathway/GO:BP/mmusculus.GO:BP.name.filter.by.MCA2.2y.genes.large2.gmt")

num <- c(which(names(geneSets)%in%GO_term))
geneSets.use <- geneSets[num]

GO_anno <- sapply(GO_term, function(x){
  temp <- geneSets[[x]]
  name <- data.frame(anno=temp@shortDescription,ID=temp@setName)
  return(name)
})

GO_anno <- t(GO_anno)
GO_anno <- GO_anno[!duplicated(GO_anno),]
GO_anno <- as.data.frame(GO_anno)
exprMat <- as.matrix(GetAssayData(Epi.use@assays$RNA,slot="data"))
nCores=4
aucellRankings=AUCell_buildRankings(exprMat, nCores=nCores, plotStats=FALSE)
# abline(v=500, col="skyblue3", lwd=5, lty=3) # aucellRankings@nGenesDetected["5%"]

regulonAUC=AUCell_calcAUC(geneSets.use, aucellRankings, aucMaxRank=500, nCores=nCores) # aucellRankings@nGenesDetected["5%"]
regulonMatix=getAUC(regulonAUC)
regulonMatix[1:5,1:5]
dim(regulonMatix)

use <- melt(regulonMatix)
anno <- data.frame(stage=Epi.use$stage,cell=colnames(Epi.use))
final <- merge(use,anno,by.x='cells',by.y='cell')
colnames(final)[2] <- 'genesets'
final2 <- aggregate(final,by=list(final$genesets,final$stage),mean)
final2 <- final2[,c(1,2,5)]

final3 <- dcast(final2,Group.1~Group.2)
rownames(final3) <- final3$Group.1
final3 <- final3[,-1]
final3 <- scale(t(final3))

library(scales)
final3 <- apply(final3, 1, rescale,to=c(-1,1))
GO_anno <- GO_anno[rownames(final3),]
rownames(final3) <- GO_anno$anno
colnames(final3) <- c("E14.5","Neonatal","10d","3w","6w~8w",
                      "12m","18m","24m")

annotation <- data.frame(name=colnames(final3),stage=colnames(final3))
rownames(annotation) <- annotation$name
annotation <- annotation['stage']
library(RColorBrewer)
col_flg<-colorRampPalette(brewer.pal(11,"RdBu"))(8)
col_flg <- rev(col_flg)
anno.col <- list(
  
  stage=c(
    E14.5 = col_flg[1],
    Neonatal = col_flg[2],
    '10d' = col_flg[3],
    '3w' = col_flg[4],
    '6w~8w' = col_flg[5],
    '12m' = col_flg[6],
    '18m' = col_flg[7],
    '24m' = col_flg[8]
    
  )
  
)
coul <- colorRampPalette(brewer.pal(9, "YlGn"))(100)
pdf("./Kidney_Auc_pheatmap.pdf",width=8,height =10)
pheatmap::pheatmap(final3,cluster_cols  = F,border_color = "white",
                   clustering_method = "ward.D2",annotation_col = annotation,show_colnames = F,
                   cellwidth = 15,cellheight = 15,annotation_colors = anno.col,
                   #color = colorRampPalette(c("#4575B4" ,"#91BFDB" ,"#E0F3F8", "white" ,"#FEE090" ,"#FC8D59" ,"#D73027"))(100))
color=coul)
dev.off()
save(regulonAUC,regulonMatix,file="./AUCell_result.RData")
