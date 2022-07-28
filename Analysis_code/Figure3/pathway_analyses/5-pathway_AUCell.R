# calculate the activity scores of gene lists

library(AUCell)
library(GSEABase)
library(data.table)

rm(list=ls())
gc()

gmtFile <- 'mmusculus.KEGG.name.filter.by.MCA2.2y.genes.large2.gmt' #gmt file is downloaded from g:Profiler(https://biit.cs.ut.ee/gprofiler/)
geneSets <- getGmt(gmtFile)

input <- './'
utissue <- c('Brain', 'Heart', 'Intestine', 'Kidney', 'Liver', 'Lung', 'Pancreas', 'Stomach', 'Testis', 'Uterus', 'Bladder', 'Spleen', 'Thymus', 'Prostate')
output <- './KEGG/'

nCores=1

for (i in 1:length(utissue)) {
  ttissue=utissue[i]
  tfilename=paste0(input,ttissue,"/",ttissue,'.exprs.norm.csv') #normalized data processed by Seurat
  afilename=paste0(input,ttissue,"/",ttissue,'.anno.csv') #annotation file includes stage and tissue
  nfile<-read.csv(tfilename,row.names = 1)
  nanno<-read.csv(afilename,row.names = 1)
  toutput=paste0(output,ttissue)
  dir.create(toutput,showWarning=F)
  setwd(toutput)
  
  stage=nanno[,'stage']
  ustage=as.character(unique(stage))
  aucdata=data.frame()
  aucanno=vector()
  for (n in 1:length(ustage)) {
    tstage=ustage[n]
    spos=which(stage %in% tstage)
    exprMat=as.matrix(nfile[,spos])
    fanno=nanno[spos,]
    aucellRankings=AUCell_buildRankings(exprMat, nCores=nCores, plotStats=FALSE)
    regulonAUC=AUCell_calcAUC(geneSets, aucellRankings, aucMaxRank=500, nCores=nCores) 
    regulonMatix=getAUC(regulonAUC)
    save(regulonAUC,file=paste0(tstage,'.',ttissue,'.regulonAUC.top500.RData'))
    fwrite(as.data.frame(regulonMatix), paste0(tstage,'.',ttissue,'.regulonAUC.top500.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
    write.table(fanno,paste0(tstage,'.',ttissue,'.regulonAUC.top500.anno.txt'),col.names=NA,sep='\t',quote=F)
    
    tmp=merge(aucdata,regulonMatix,by='row.names',all=T)
    tmp[is.na(tmp)]=0
    body=tmp[,2:ncol(tmp)]
    rownames(body)=tmp[,1]
    aucdata=body;
    aucanno=rbind(aucanno,fanno)
    
    rm(exprMat, aucellRankings, regulonAUC); gc();
    print(paste0(ttissue,': ',tstage,' finished'));
  }
  setwd(output)
  save(aucdata,file=paste0(ttissue,'.regulonAUC.top500.RData'))
  fwrite(as.data.frame(aucdata),paste0(ttissue,'.regulonAUC.top500.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
  write.table(aucanno,paste0(ttissue,'.regulonAUC.top500.anno.txt'),col.names=NA,sep='\t',quote=F)
  
  rm(nfile,aucdata,aucanno,tmp,body); gc();
}


# plot the pathway

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(reshape2)

anno.dir <- sort(list.files(path = ".", pattern = "anno.txt"))
auc.dir <- sort(list.files(path = ".", pattern = ".RData"))

fig.total <- list()
for (i in 1:length(anno.dir)) {
  setwd("./")
  load(auc.dir[i])
  tissue <- strsplit(anno.dir[i], "\\.")[[1]][1]
  anno <- read.table(anno.dir[i], sep = "\t", row.names = 1)
  path <- rownames(aucdata)
  
  # analyze oxidative phosphorylation pathway
  select.data <- aucdata["path:mmu00190", ]
  
  # calculate mean AUC(>0) for each cell
  result <- as.data.frame(apply(select.data, 2, sum))
  colnames(result) <- "sum"
  result$count <- 1
  result$avg <- result$sum / result$count
  
  # merge by stage
  result$stage <- anno[rownames(result), "stage"]
  
  res <- aggregate(result$avg, by = list(result$stage), FUN = mean)
  colnames(res) <- c("stage", "value")
  
  # set stage order
  stage.order <- c("Fetal", "Neonatal", "TenDays", "ThreeWeeks",
                   "Adult", "OneYear", "EighteenMonths", "TwoYears")
  
  # stage to num
  stage.num <- c(1:8)
  names(stage.num) <- stage.order
  res$stage <- factor(res$stage, levels = stage.order)
  res$num <- stage.num[res$stage]
  
  #plot figure
  fig <- ggscatter(res, x = "num", y = "value", 
                   add = "reg.line", 
                   add.params = list(color = "#E69F00", fill = "lightgray"),
                   conf.int = TRUE,
                   cor.coef = TRUE,
                   cor.coeff.args =
                     list(method = "pearson",
                          label.x.npc = "left",
                          label.y.npc = "top"
                     ),
                   cor.coef.size = 3,
  ) +
    scale_x_discrete(limits = factor(c(1:8)), labels = stage.order) +
    labs(x = "Stage", y = "Activity", title = tissue) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      title = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"))
  
  fig
  
  # pdf
  setwd("./results/oxidative_phosphorylation/")
  out.dir <- paste(tissue, ".pdf", sep = "")
  ggsave(fig, filename = out.dir, width = 4, height = 4)
  
  fig.total[[i]] <- fig
  
}

ggfig.total[['nrow']] <- 2
fig.total[['ncol']] <- 5

pdf('./all.tissues.KEGG:mmu00190.pdf', width = 25, height = 10)
do.call('grid.arrange', fig.total)
dev.off()
