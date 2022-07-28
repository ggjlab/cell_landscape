# calculate the activity scores of gene lists
setwd("/media/ggj/ggjlab/RData/agingatlas/MCA2/all/pathway/GO:BP/")
library(AUCell)
library(GSEABase)
library(data.table)

rm(list=ls())
gc()

#gmtFile <- 'mmusculus.GO:BP.name.filter.by.MCA2.2y.genes.large2.gmt' #gmt file is downloaded from g:Profiler(https://biit.cs.ut.ee/gprofiler/)
#geneSets <- getGmt(gmtFile)

input <- '/media/ggj/ggjlab/RData/agingatlas/MCA2/all/gene/'
output <- './use_14tissue/'

nCores=8
load("../../MAST_time_related_gene/result/result_14tissue_1e10.RData")
utissue <- unique(upgene$tissue)
for (i in 1:length(utissue)) {
  temp.up <- upgene[upgene$tissue==unique(upgene$tissue)[i],]
  temp.up <- temp.up$gene
  temp.down <- downgene[downgene$tissue==unique(upgene$tissue)[i],]
  temp.down <- temp.down$gene
  names(temp.up) <- rep("up",times=length(temp.up))
  names(temp.down) <- rep("down",times=length(temp.down))
  temp <- list(up=temp.up,
               down=temp.down)
  lst=lapply(names(temp), function(x){
    GeneSet(temp[[x]],setName=x)
  })
  geneSets <- GeneSetCollection(lst)
  
  setwd("/media/ggj/ggjlab/RData/agingatlas/MCA2/all/pathway/GO:BP/")
  ttissue=utissue[i]
  tfilename=paste0(input,ttissue,"/",ttissue,'.exprs.norm.csv') #normalized data processed by Seurat
  afilename=paste0(input,ttissue,"/",ttissue,'.anno.csv') #annotation file includes stage and tissue
  nfile<-fread(tfilename)
  nfile <- as.data.frame(nfile)
  rownames(nfile) <- nfile$V1
  nfile <- nfile[,-1]
  
  nanno<-read.csv(afilename,row.names = 1)
  toutput=paste0(output,ttissue)
  dir.create(toutput,showWarning=T)
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
  setwd("../../")
  setwd(output)
  save(aucdata,file=paste0(ttissue,'.regulonAUC.top500.RData'))
  fwrite(as.data.frame(aucdata),paste0(ttissue,'.regulonAUC.top500.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
  write.table(aucanno,paste0(ttissue,'.regulonAUC.top500.anno.txt'),col.names=NA,sep='\t',quote=F)
  
  rm(nfile,aucdata,aucanno,tmp,body); gc();
}


