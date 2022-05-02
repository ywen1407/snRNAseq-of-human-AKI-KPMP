library(ggplot2,verbose=F)
library(Seurat,verbose=F)
library(dplyr,verbose=F)
library(SeuratData,verbose=F)
library(SeuratDisk,verbose=F)
library(tidyverse,verbose=F)
library("data.table",verbose=F)
library(scmap)
library(SingleCellExperiment,verbose=F)
rm(list=ls())
gc()

memory.limit(size=99999999)
#read database and create merged label
human<- LoadH5Seurat("~/human4integration_01042022.h5seurat")
ident.list<-unique(human$orig.ident)
human_hvg<-c()
for (i in 1:length(ident.list)){
  temp<-subset(human,subset=orig.ident==ident.list[i])
  temp <- NormalizeData(temp, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
  temp <- FindVariableFeatures(temp, selection.method = "vst", nfeatures =500,verbose = F)
  human_hvg<-c(human_hvg,temp@assays$RNA@var.features)
}
human_hvg<-unique(human_hvg)
length(human_hvg)
human@assays$RNA@var.features<-human_hvg
saveRDS(human,"~/human4scmap_0104.rds")

mouse <- LoadH5Seurat("~/mouse4integration_01042022.h5seurat")
mouse
ident.list<-unique(mouse$orig.ident)
mouse_hvg<-c()
for (i in 1:length(ident.list)){
  temp<-subset(mouse,subset=orig.ident==ident.list[i])
  temp <- NormalizeData(temp, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
  temp <- FindVariableFeatures(temp, selection.method = "vst", nfeatures =500,verbose = F)
  mouse_hvg<-c(mouse_hvg,temp@assays$RNA@var.features)
}
mouse_hvg<-unique(mouse_hvg)
length(mouse_hvg)
mouse@assays$RNA@var.features<-mouse_hvg
saveRDS(mouse,"~/mouse4scmap_0104.rds")

human<-read_rds("~/human4scmap_0104.rds")
mouse<-read_rds("~/mouse4scmap_0104.rds")
#human: case_type,0=AKI, 1=healthy control
#mouse: Group, Control 4hours  12hours 2days   14days  6weeks
human$celltype_merge<-plyr::mapvalues(human$nuctype2,from=unique(human$nuctype2),
                                      to=c("DCT","PT_Healthy","PT_Injury","Fib","TAL","Lyjmphocyte","PC",
                                           "EC","ICA","PT_Healthy","Pericyte","Myeloid","PT_Injury","ICB","POD",
                                           "PEC",
                                           "DTL_ATL","PT_Healthy"))
human$celltype_merge2<-human$nuctype3
mouse$celltype_merge<-plyr::mapvalues(mouse$celltype2,from=unique(mouse$celltype2),
                                      to=c("PT_Healthy","PT_Healthy","Fib","EC","PT_Injury",
                                           "DCT","TAL","TAL","DTL_ATL","ICB","CNT","DCT","ICA","TAL","Myeloid","PC",
                                           "Pericyte","PC","Lymphocyte","POD","PEC","URO","TAL","EC",
                                           "PT_Injury","PT_Injury","PT_Injury","PT_Injury"))
mouse$celltype_merge2<-plyr::mapvalues(mouse$celltype2,from=unique(mouse$celltype2),
                                       to=c("PT","PT","Fib","EC","PT",
                                            "DCT","TAL","TAL","DTL_ATL","ICB","CNT","DCT","ICA","TAL","Myeloid","PC",
                                            "Pericyte","PC","Lymphocyte","POD","PEC","URO","TAL","EC",
                                            "PT","PT","PT","PT"))

#convert to singelcellexperiment
#repeat following code for different subgroups described in our manuscript. 
human2<-human
mouse2<-mouse
human.sce <- as.SingleCellExperiment(human2)
mouse.sce <- as.SingleCellExperiment(mouse2)

rowData(human.sce)$feature_symbol <- rownames(human.sce)
rowData(mouse.sce)$feature_symbol <- rownames(mouse.sce)

rowData(human.sce)$scmap_features<-FALSE
rowData(mouse.sce)$scmap_features<-FALSE
rowData(human.sce)[rowData(human.sce)$feature_symbol %in% human2@assays$RNA@var.features,'scmap_features']<-TRUE
rowData(mouse.sce)[rowData(mouse.sce)$feature_symbol %in% mouse2@assays$RNA@var.features,'scmap_features']<-TRUE

human.sce <- indexCluster(human.sce,cluster_col="celltype_merge2")
mouse.sce <- indexCluster(mouse.sce,cluster_col="celltype_merge2")

set.seed(30729)
scmapCluster_results <- scmapCluster(
  projection = mouse.sce, 
  index_list = list(
    ownprojection = metadata(human.sce)$scmap_cluster_index
  ),threshold = 0
)
plot(
  getSankey(
    colData(mouse.sce)$celltype_merge2, 
    scmapCluster_results$scmap_cluster_labs[,'ownprojection'],
    plot_height = 400
  )
)
zoomout<-as.data.frame(cbind(as.character(colData(mouse.sce)$celltype_merge2),scmapCluster_results$scmap_cluster_labs))
write.csv(zoomout,"~/global.csv")

#PT subclusters
human<-read_rds("~/human4scmap_0104.rds")
mouse<-read_rds("~/mouse4scmap_0104.rds")

#rename celltype_merge based on annotation described in our paper. 
human$celltype_merge<-plyr::mapvalues(human$nuctype2,from=unique(human$nuctype2),
                                      to=c("DCT","PT_Healthy","PT_FR","Fib","TAL","Lyjmphocyte","PC",
                                           "EC","ICA","PT_Healthy","Pericyte","Myeloid","PT_Injury","ICB","POD",
                                           "PEC",
                                           "DTL_ATL","PT_Healthy"))
human$celltype_merge2<-human$nuctype3
mouse$celltype_merge<-plyr::mapvalues(mouse$celltype2,from=unique(mouse$celltype2),
                                      to=c("PT_Healthy","PT_Healthy","Fib","EC","PT_FR",
                                           "DCT","TAL","TAL","DTL_ATL","ICB","CNT","DCT","ICA","TAL","Myeloid","PC",
                                           "Pericyte","PC","Lymphocyte","POD","PEC","URO","TAL","EC",
                                           "PT_Injury","PT_Injury","PT_Injury","PT_Injury"))
mouse$celltype_merge2<-plyr::mapvalues(mouse$celltype2,from=unique(mouse$celltype2),
                                       to=c("PT","PT","Fib","EC","PT",
                                            "DCT","TAL","TAL","DTL_ATL","ICB","CNT","DCT","ICA","TAL","Myeloid","PC",
                                            "Pericyte","PC","Lymphocyte","POD","PEC","URO","TAL","EC",
                                            "PT","PT","PT","PT"))
humanpt<-subset(human,subset=celltype_merge2=="PT")
mousept<-subset(mouse,subset=celltype_merge2=="PT")

humanpt <- NormalizeData(humanpt, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
humanpt<-FindVariableFeatures(humanpt, selection.method = "vst", nfeatures = 1500,verbose=F)
mousept <- NormalizeData(mousept, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
mousept<-FindVariableFeatures(mousept, selection.method = "vst", nfeatures = 1500,verbose=F)

#repeat following code as described in our paper
humanpt2<-humanpt
mousept2<-mousept


humanpt.sce <- as.SingleCellExperiment(humanpt2)
mousept.sce <- as.SingleCellExperiment(mousept2)

rowData(humanpt.sce)$feature_symbol <- rownames(humanpt.sce)
rowData(mousept.sce)$feature_symbol <- rownames(mousept.sce)

rowData(humanpt.sce)$scmap_features<-FALSE
rowData(mousept.sce)$scmap_features<-FALSE
rowData(humanpt.sce)[rowData(humanpt.sce)$feature_symbol %in% humanpt2@assays$RNA@var.features,'scmap_features']<-TRUE
rowData(mousept.sce)[rowData(mousept.sce)$feature_symbol %in% mousept2@assays$RNA@var.features,'scmap_features']<-TRUE

humanpt.sce <- indexCluster(humanpt.sce,cluster_col="celltype_merge")
mousept.sce <- indexCluster(mousept.sce,cluster_col="celltype_merge")

set.seed(30729)
scmapCluster_results <- scmapCluster(
  projection = mousept.sce, 
  index_list = list(
    ownprojection = metadata(humanpt.sce)$scmap_cluster_index
  ),threshold = 0
)
plot(
  getSankey(
    colData(mousept.sce)$celltype_merge, 
    scmapCluster_results$scmap_cluster_labs[,'ownprojection'],
    plot_height = 400
  )
)
zoomin<-as.data.frame(cbind(as.character(colData(mousept.sce)$celltype_merge),scmapCluster_results$scmap_cluster_labs))
write.csv(zoomin,"C:/Users/kastu/Downloads/pt_all.csv")

#validation in human
set.seed(307289)

human$replicate <- sample(c("rep1", "rep2","rep3","rep4"), size = ncol(human), replace = TRUE)
human1<-subset(human,subset=replicate!="rep4")
human2<-subset(human,subset=replicate=="rep4")

human1.sce <- as.SingleCellExperiment(human1)
human2.sce <- as.SingleCellExperiment(human2)

rowData(human1.sce)$feature_symbol <- rownames(human1.sce)
rowData(human2.sce)$feature_symbol <- rownames(human2.sce)

rowData(human1.sce)$scmap_features<-FALSE
rowData(human2.sce)$scmap_features<-FALSE
rowData(human1.sce)[rowData(human1.sce)$feature_symbol %in% human1@assays$RNA@var.features,'scmap_features']<-TRUE
rowData(human2.sce)[rowData(human2.sce)$feature_symbol %in% human2@assays$RNA@var.features,'scmap_features']<-TRUE

human1.sce <- indexCluster(human1.sce,cluster_col="celltype_merge2")
human2.sce <- indexCluster(human2.sce,cluster_col="celltype_merge2")

set.seed(30729)
scmapCluster_results <- scmapCluster(
  projection = human2.sce, 
  index_list = list(
    ownprojection = metadata(human1.sce)$scmap_cluster_index
  ),threshold = 0
)
plot(
  getSankey(
    colData(human2.sce)$celltype_merge2, 
    scmapCluster_results$scmap_cluster_labs[,'ownprojection'],
    plot_height = 400
  )
)
zoomout<-as.data.frame(cbind(as.character(colData(human2.sce)$celltype_merge2),scmapCluster_results$scmap_cluster_labs))
write.csv(zoomout,"~/global_val.csv")

humanpt<-subset(human,subset=celltype_merge2=="PT")
humanpt <- NormalizeData(humanpt, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
humanpt<-FindVariableFeatures(humanpt, selection.method = "vst", nfeatures = 1500,verbose=F)
set.seed(307289)
humanpt$replicate <- sample(c("rep1", "rep2","rep3","rep4"), size = ncol(humanpt), replace = TRUE)
humanpt1<-subset(humanpt,subset=replicate!="rep4")
humanpt2<-subset(humanpt,subset=replicate=="rep4")

humanpt1.sce <- as.SingleCellExperiment(humanpt1)
humanpt2.sce <- as.SingleCellExperiment(humanpt2)

rowData(humanpt1.sce)$feature_symbol <- rownames(humanpt1.sce)
rowData(humanpt2.sce)$feature_symbol <- rownames(humanpt2.sce)

rowData(humanpt1.sce)$scmap_features<-FALSE
rowData(humanpt2.sce)$scmap_features<-FALSE
rowData(humanpt1.sce)[rowData(humanpt1.sce)$feature_symbol %in% humanpt1@assays$RNA@var.features,'scmap_features']<-TRUE
rowData(humanpt2.sce)[rowData(humanpt2.sce)$feature_symbol %in% humanpt2@assays$RNA@var.features,'scmap_features']<-TRUE

humanpt1.sce <- indexCluster(humanpt1.sce,cluster_col="celltype_merge")
humanpt2.sce <- indexCluster(humanpt2.sce,cluster_col="celltype_merge")

set.seed(30729)
scmapCluster_results <- scmapCluster(
  projection = humanpt2.sce, 
  index_list = list(
    ownprojection = metadata(humanpt1.sce)$scmap_cluster_index
  ),threshold = 0
)
plot(
  getSankey(
    colData(humanpt2.sce)$celltype_merge, 
    scmapCluster_results$scmap_cluster_labs[,'ownprojection'],
    plot_height = 400
  )
)
zoomin<-as.data.frame(cbind(as.character(colData(humanpt2.sce)$celltype_merge),scmapCluster_results$scmap_cluster_labs))
write.csv(zoomin,"C:/Users/kastu/Downloads/pt_val.csv")
