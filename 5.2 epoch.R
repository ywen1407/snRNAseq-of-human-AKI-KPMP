library(Seurat)
library(SeuratDisk)
library(dplyr)
library(epoch)

rm(list=ls())
gc()
memory.limit(size = 9999999999)

Convert("C:/Users/kastu/Downloads/KPMP/snRNA/pt_01042022_4110hvg.h5ad", dest = "h5seurat", overwrite = TRUE)

pt<-LoadH5Seurat("~/pt_01042022_4110hvg.h5Seurat")
TF<-utils_loadObject("~/hsTFs.rda")
expDat<-as.matrix(pt[["RNA"]]@data)
sampTab<-pt@meta.data
expDat2<-expDat[rowSums(expDat)>0,]
dim(expDat2)
dgenes<-rownames(expDat2)
grnDF <- reconstructGRN(expDat2, TF, dgenes, method="pearson", zThresh=3)
grnDF<-grnDF[,c(2,1,3,4)]
grnDF<-grnDF[,c(1,2,3)]
colnames(grnDF)<-c("TF","target","importance")
write.csv(grnDF,"~/grnDF_epoch_MI_pt_0104_5301hvg.csv",row.names = F)
