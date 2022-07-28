# average several cells for different tissue cells from different stages
rm(list=ls())
library(data.table)
library(Seurat)
print('Begin load expression matrix.')
exprsinput='/media/ggj/ggjlab/RData/agingatlas/MCA2/all/'
setwd(exprsinput)
#filename='rawdata.rdata' 
filename='./mca_10stages_20211215_sct.rdata' 
load(filename)
print('All expression data loaded successfully.')
adata=data

annoinput='/media/ggj/ggjlab/RData/agingatlas/MCA2/all/';
setwd(annoinput)
anno=read.csv('./anno/Mouse_annotation_20220120.csv' ,sep="\t",header=T, row.names=1, check.names=F,stringsAsFactors = F)
colnames(anno)[6] <- "leiden"
anno$leiden<-paste0("cluster",anno$leiden)
anno$cellID<-rownames(anno)
#write.csv(anno,file="./cell_anno_fromh5ad.csv")
sample1=colnames(adata)
sample2=rownames(anno)

data=adata[,sample2]
rm(adata)
gc()
# tissue=anno[,'tissue']
# # utissue=as.character(unique(tissue))
# unique(tissue)
# utissue=as.character(unique(tissue))
# utissue=c("Bladder","Brain","Heart","Intestine","Kidney","Liver","Lung","Pancreas","Spleen","Stomach","Testis","Thymus","Uterus","Prostate")
# #utissue=c("Bladder","Brain","Heart","Intestine","Kidney","Liver","Lung","Pancreas","Spleen","Stomach","Testis","Thymus","Uterus")
utissue<-unique(anno$leiden)
tissue<-anno[,'leiden']
print('Begin to aggregate the data.')

output='/media/ggj/ggjlab/RData/agingatlas/MCA2/all/pseudocell_cluster/'

### dge and anno file for each tissue and each stage 
for (i in 1:length(utissue)) {
  avg=vector()
  avg_anno=vector()
  tmp=utissue[i]
  toutput=paste0(output,tmp,'/')
  dir.create(toutput,showWarning=F)
  setwd(toutput)
  
  pos=which(tissue %in% tmp)
  ttdata=data[,pos]
  ttdata<-GetAssayData(ttdata@assays$RNA,slot="counts")
  min<-floor(length(colnames(ttdata))*0.001)
  pbmc<-CreateSeuratObject(ttdata,min.cells = min)
  ttdata<-as.matrix(GetAssayData(pbmc@assays$RNA,slot="counts"))
  tanno=anno[pos,]
  savefilename=paste0(utissue[i],".RData")
  annofilename=paste0(utissue[i],"_anno.csv")
  dgefilename=paste0(utissue[i],"_dge.csv")
  tdata=as.data.frame(ttdata)
  ct=tanno[,c('leiden','Tissue','Stage')]
  #save(ct,tdata,file=savefilename)
  fwrite(ct,file=annofilename,row.names = T)
  fwrite(tdata,file=dgefilename,row.names = T)
  message(utissue[i]," has finished!")
}


##############compute pseudocell in python################



# combine all pseudo cells together 
library(data.table)

rm(list=ls()); gc()
input='/media/ggj/ggjlab/RData/agingatlas/MCA2/all/pseudocell_cluster/'
output=paste0(input,'allTissues/')
dir.create(output,showWarning=F)
setwd(output)
#atissue=c('Brain', 'Heart', 'Intestine', 'Kidney', 'Liver', 'Lung', 'Pancreas', 'Stomach', 'Testis', 'Uterus', 'Bladder', 'Spleen', 'Thymus', 'Prostate');

favg=vector()
favg_anno=vector()
setwd("../")
temp_pseudocell<-list.files(pattern = "*dge_new.csv")
temp_pseudocell_anno<-list.files(pattern = "*anno_new.csv")


#dge_pseudocell
all_tissue_stage_pseudocell<-NULL
for(i in 1:length(temp_pseudocell)){
  tempvalue<-read.csv(temp_pseudocell[i],check.names = F,row.names = 1,stringsAsFactors = F)
  tempvalue<-t(tempvalue)
  all_tissue_stage_pseudocell<-merge(tempvalue,all_tissue_stage_pseudocell,by="row.names",all=T)
  rownames(all_tissue_stage_pseudocell)<-all_tissue_stage_pseudocell$Row.names
  all_tissue_stage_pseudocell<-all_tissue_stage_pseudocell[,-1]
  all_tissue_stage_pseudocell[is.na(all_tissue_stage_pseudocell)]<-0
  message(temp_pseudocell[i]," has finished!")
}
write.csv(all_tissue_stage_pseudocell,file="allTissues/all_tissue_stage_pseudocell100.csv")

#anno_pseudocell
all_tissue_stage_pseudocell_anno<-NULL
for(i in 1:length(temp_pseudocell_anno)){
  tempvalue<-read.csv(temp_pseudocell_anno[i],check.names = F,row.names = 1,stringsAsFactors = F)
  all_tissue_stage_pseudocell_anno<-rbind(all_tissue_stage_pseudocell_anno,tempvalue)
  message(temp_pseudocell[i]," has finished!")
}
all_tissue_stage_pseudocell_anno$cellID<-rownames(all_tissue_stage_pseudocell_anno)
write.csv(all_tissue_stage_pseudocell_anno,file="allTissues/all_tissue_stage_pseudocell100_anno.csv")

###only stage>4 tissue
atissue=c('Brain', 'Heart', 'Intestine', 'Kidney', 'Liver', 'Lung', 'Pancreas', 'Stomach', 'Testis', 'Uterus', 'Bladder', 'Spleen', 'Thymus', 'Prostate');
tissue=c(atissue)
input='/media/ggj/ggjlab/RData/agingatlas/MCA2/all/pseudo_result/'
output=paste0(input,'allTissues/')
dir.create(output,showWarning=F)
setwd(output)

favg=vector()
favg_anno=vector()
setwd("../tissue_use/")
#dge_pseudocell
temp_pseudocell<-list.files(pattern = "*dge_new.csv")
temp_pseudocell_anno<-list.files(pattern = "*anno_new.csv")
all_tissue_stage_pseudocell<-NULL
for(i in 1:length(temp_pseudocell)){
  tempvalue<-read.csv(temp_pseudocell[i],check.names = F,row.names = 1,stringsAsFactors = F)
  tempvalue<-t(tempvalue)
  all_tissue_stage_pseudocell<-merge(tempvalue,all_tissue_stage_pseudocell,by="row.names",all=T)
  rownames(all_tissue_stage_pseudocell)<-all_tissue_stage_pseudocell$Row.names
  all_tissue_stage_pseudocell<-all_tissue_stage_pseudocell[,-1]
  all_tissue_stage_pseudocell[is.na(all_tissue_stage_pseudocell)]<-0
  message(temp_pseudocell[i]," has finished!")
}
write.csv(all_tissue_stage_pseudocell,file="../allTissues/all_tissue_stage_pseudocell100.csv")
#anno_pseudocell
all_tissue_stage_pseudocell_anno<-NULL
for(i in 1:length(temp_pseudocell_anno)){
  tempvalue<-read.csv(temp_pseudocell_anno[i],check.names = F,row.names = 1,stringsAsFactors = F)
  all_tissue_stage_pseudocell_anno<-rbind(all_tissue_stage_pseudocell_anno,tempvalue)
  message(temp_pseudocell[i]," has finished!")
}
all_tissue_stage_pseudocell_anno$cellID<-rownames(all_tissue_stage_pseudocell_anno)
write.csv(all_tissue_stage_pseudocell_anno,file="../allTissues/all_tissue_stage_pseudocell100_anno.csv")
##

# for (i in 1:length(utissue)) {
#   avg=vector()
#   avg_anno=vector()
#   tmp=utissue[i]
#   toutput=paste0(output,tmp,'/')
#   dir.create(toutput,showWarning=F)
#   setwd(toutput)
#   
#   pos=which(tissue %in% tmp)
#   ttdata=data[,pos]
#   ttdata<-GetAssayData(ttdata@assays$RNA,slot="counts")
#   min<-floor(length(colnames(ttdata))*0.001)
#   pbmc<-CreateSeuratObject(ttdata,min.cells = min)
#   ttdata<-as.matrix(GetAssayData(pbmc@assays$RNA,slot="counts"))
#   dgefilename=paste0(utissue[i],"_allstage","_dge.csv")
#   write.csv(ttdata,file=dgefilename,row.names = T)
#  
#   message(utissue[i]," has finished!")
# }


