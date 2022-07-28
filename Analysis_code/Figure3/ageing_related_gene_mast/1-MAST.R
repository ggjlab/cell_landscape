setwd("/media/ggj/ggjlab/RData/agingatlas/MCA2/all/")
library(MAST)
library(readxl)
library(Seurat)
library(data.table)
dir.create("MAST_time_related_gene")
anno<-read.table("./anno/Mouse_accurate_annotation_14tissue_20220401.csv",sep=",",row.names = 1,header = T)

tissue <- unique(anno$tissue)
#options(mc.cores=8)
for(i in c(10:length(tissue))){
  file <- paste0("./gene/",tissue[i],"/",tissue[i],".exprs.norm.csv")
  cpm <- fread(file)
  cpm <- as.data.frame(cpm)
  cpm[1:5,1:5]
  rownames(cpm) <- cpm$V1
  cpm <- cpm[,-1]
  
  temp.anno <- anno[anno$tissue==tissue[i]&anno$stage%in%c("Adult","OneYear","EighteenMonths","TwoYears"),]
  cpm <- cpm[,temp.anno$cellID]
  
  stage_to_num <- c(1:4)
  names(stage_to_num) <- c("Adult","OneYear","EighteenMonths","TwoYears")
  
  #######MAST
  print(table(temp.anno$stage))
  #cells <- a@meta.data[a$nFeature_RNA>=500,]
  #print(table(cells$Stage))
  #a <- a[,rownames(cells)]
  age_num = as.numeric(stage_to_num[temp.anno$stage])  
  #log_counts <- a@assays$RNA@data
  log_counts <- as.matrix(cpm)
  gc()
  log_counts[log_counts < 0.1] <- 0
  log_counts <- as(log_counts, "dgCMatrix")
  
  fData = data.frame(names=rownames(log_counts))
  rownames(fData) = rownames(log_counts)
  
  
  
  cData = data.frame(cond=age_num)
  rownames(cData) = colnames(log_counts)
  
  s = as.matrix(log_counts)
  obj <- FromMatrix(s, cData, fData)
  colData(obj)$cngeneson <- scale(colSums(assay(obj)>0))
  zlmCond <- zlm(formula = as.formula("~cond+cngeneson"), obj)
  
  summaryCond <- summary(zlmCond, doLRT="cond")
  summaryDt <- summaryCond$datatable
  dt1 = summaryDt[contrast=="cond" & component=="H", .(primerid, `Pr(>Chisq)`)]
  dt2 = summaryDt[contrast=="cond" & component=="logFC", .(primerid, coef, z)]
  de_res = merge(dt1, dt2, by="primerid")
  colnames(de_res) <- c("gene", "age.H_p", "age.logFC", 'age.logFC_z')
  de_res$age.H_fdr <- p.adjust(de_res$age.H_p, "fdr")
  de_res = de_res[order(-de_res$age.logFC),]
  de_res <- de_res[de_res$age.H_fdr<0.01&de_res$age.logFC!="NaN",]
  
  file1 <- paste0("./MAST_time_related_gene/",tissue[i],"_result.RData")
  save(de_res,summaryCond,summaryDt,zlmCond,file=file1)
  
  message(tissue[i]," finished!")
  gc()
  }


# 4 aging tissue ----------------------------------------------------------

tissue <- c("Bladder","Prostate","Spleen","Thymus")
for(i in c(1:length(tissue))){
  file <- paste0("./gene/",tissue[i],"/",tissue[i],".exprs.norm.csv")
  cpm <- fread(file)
  cpm <- as.data.frame(cpm)
  cpm[1:5,1:5]
  rownames(cpm) <- cpm$V1
  cpm <- cpm[,-1]
  
  temp.anno <- anno[anno$tissue==tissue[i]&anno$stage%in%c("Adult","OneYear","EighteenMonths","TwoYears"),]
  cpm <- cpm[,temp.anno$cellID]
  
  stage_to_num <- c(1:4)
  names(stage_to_num) <- c("Adult","OneYear","EighteenMonths","TwoYears")
  
  #######MAST
  print(table(temp.anno$stage))
  #cells <- a@meta.data[a$nFeature_RNA>=500,]
  #print(table(cells$Stage))
  #a <- a[,rownames(cells)]
  age_num = as.numeric(stage_to_num[temp.anno$stage])  
  #log_counts <- a@assays$RNA@data
  log_counts <- as.matrix(cpm)
  gc()
  log_counts[log_counts < 0.1] <- 0
  log_counts <- as(log_counts, "dgCMatrix")
  
  fData = data.frame(names=rownames(log_counts))
  rownames(fData) = rownames(log_counts)
  
  
  
  cData = data.frame(cond=age_num)
  rownames(cData) = colnames(log_counts)
  
  s = as.matrix(log_counts)
  obj <- FromMatrix(s, cData, fData)
  colData(obj)$cngeneson <- scale(colSums(assay(obj)>0))
  zlmCond <- zlm(formula = as.formula("~cond+cngeneson"), obj)
  
  summaryCond <- summary(zlmCond, doLRT="cond")
  summaryDt <- summaryCond$datatable
  dt1 = summaryDt[contrast=="cond" & component=="H", .(primerid, `Pr(>Chisq)`)]
  dt2 = summaryDt[contrast=="cond" & component=="logFC", .(primerid, coef, z)]
  de_res = merge(dt1, dt2, by="primerid")
  colnames(de_res) <- c("gene", "age.H_p", "age.logFC", 'age.logFC_z')
  de_res$age.H_fdr <- p.adjust(de_res$age.H_p, "fdr")
  de_res = de_res[order(-de_res$age.logFC),]
  de_res <- de_res[de_res$age.H_fdr<0.01&de_res$age.logFC!="NaN",]
  
  file1 <- paste0("./MAST_time_related_gene/",tissue[i],"_result.RData")
  save(de_res,summaryCond,summaryDt,zlmCond,file=file1)
  
  message(tissue[i]," finished!")
  gc()
}

# # structure ---------------------------------------------------------------
# anno <- anno[anno$celllineage%in%c("Epithelial","Endothelial","Stromal"),]
# table(anno$celllineage)
# tissue <- c("Bladder","Spleen","Thymus","Prostate")
# options(mc.cores=4)
# for(i in c(4:length(tissue))){
#   file <- paste0("./gene/",tissue[i],"/",tissue[i],".exprs.norm.csv")
#   cpm <- fread(file)
#   cpm <- as.data.frame(cpm)
#   cpm[1:5,1:5]
#   rownames(cpm) <- cpm$V1
#   cpm <- cpm[,-1]
  
#   temp.anno <- anno[anno$tissue==tissue[i]&anno$stage%in%c("Adult","OneYear","EighteenMonths","TwoYears"),]
#   cpm <- cpm[,temp.anno$cellID]
  
#   stage_to_num <- c(1:4)
#   names(stage_to_num) <- c("Adult","OneYear","EighteenMonths","TwoYears")
  
#   celltype <- unique(anno$celllineage)
#   for(j in 1:length(celltype)){
#     anno.use <- temp.anno[temp.anno$celllineage==celltype[j],]
#     if(length(rownames(anno.use))<100|length(unique(anno.use$stage))==1){
#       next
#     }
#     else{
#       cpm.use <- cpm[,anno.use$cellID]
      
#       #######MAST
#       print(table(anno.use$stage))
#       #cells <- a@meta.data[a$nFeature_RNA>=500,]
#       #print(table(cells$Stage))
#       #a <- a[,rownames(cells)]
#       age_num = as.numeric(stage_to_num[anno.use$stage])  
#       #log_counts <- a@assays$RNA@data
#       log_counts <- as.matrix(cpm.use)
#       gc()
#       log_counts[log_counts < 0.1] <- 0
#       log_counts <- as(log_counts, "dgCMatrix")
      
#       fData = data.frame(names=rownames(log_counts))
#       rownames(fData) = rownames(log_counts)
      
      
      
#       cData = data.frame(cond=age_num)
#       rownames(cData) = colnames(log_counts)
      
#       s = as.matrix(log_counts)
#       obj <- FromMatrix(s, cData, fData)
#       colData(obj)$cngeneson <- scale(colSums(assay(obj)>0))
#       zlmCond <- zlm(formula = as.formula("~cond+cngeneson"), obj)
      
#       summaryCond <- summary(zlmCond, doLRT="cond")
#       summaryDt <- summaryCond$datatable
#       dt1 = summaryDt[contrast=="cond" & component=="H", .(primerid, `Pr(>Chisq)`)]
#       dt2 = summaryDt[contrast=="cond" & component=="logFC", .(primerid, coef, z)]
#       de_res = merge(dt1, dt2, by="primerid")
#       colnames(de_res) <- c("gene", "age.H_p", "age.logFC", 'age.logFC_z')
#       de_res$age.H_fdr <- p.adjust(de_res$age.H_p, "fdr")
#       de_res = de_res[order(-de_res$age.logFC),]
#       de_res <- de_res[de_res$age.H_fdr<0.05&de_res$age.logFC!="NaN",]
      
#       file1 <- paste0("./MAST_time_related_gene/structure/",tissue[i],"_",celltype[j],"_result.RData")
#       save(de_res,summaryCond,summaryDt,zlmCond,file=file1)
      
#       message(tissue[i],"_",celltype[j]," finished!")
#       gc()
#     }
   
#   }
  
# }

