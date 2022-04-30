library(RColorBrewer)
library(ggplot2,verbose=F)
library(Seurat,verbose=F)
library(dplyr,verbose=F)
library(SeuratData,verbose=F)
library(SeuratDisk,verbose=F)
library(tidyverse,verbose=F)
library("data.table",verbose=F)
library("viridis")
library(MetBrewer)

#filter out low quality nuclei and nuclei with high Mt content
allcase[['percent.mt']]<-PercentageFeatureSet(allcase,pattern="^MT-")
allcase<-subset(allcase,subset=nFeature_RNA>500 &nFeature_RNA<5000 & nCount_RNA<10000 &nCount_RNA>500 & percent.mt<2)

#get HVG from each dataset and unionize
ident.list<-unique(allcase$orig.ident)
varlist<-c()
for (i in 1:length(ident.list)){
  temp<-subset(allcase,subset=orig.ident==ident.list[i])
  temp <- NormalizeData(temp, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
  temp <- FindVariableFeatures(temp, selection.method = "vst", nfeatures =500,verbose = F)
  varlist<-c(varlist,temp@assays$RNA@var.features)
}
varlist<-unique(varlist)
length(varlist)
allcase@assays$RNA@var.features<-varlist

#preprocess
allcase <- NormalizeData(allcase, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
allcase<- ScaleData(allcase,verbose=F)
allcase<- RunPCA(allcase,verbose = F)
ElbowPlot(allcase,ndims=30)

allcase <- allcase %>%FindNeighbors(dims=1:13,verbose=F)%>%FindClusters(resolution=0.6,verbose=F)%>%identity()
allcase <- RunUMAP(allcase, dims = 1:13,verbose=F)

feature=c("SLC5A12","SLC22A6","HAVCR1","VCAM1","TACSTD2",
          "UMOD","SLC12A1","SLC12A3","AQP2","SCNN1G",
          "SLC26A7","SLC4A9","CRB2","CLDN1","NPHS2",
          'ACTA2',"DCN","MYH11", "EMCN","FLT1","CD163","IL7R","MYH2")
DotPlot(allcase, features = feature,scale=T) + RotatedAxis()+scale_color_distiller(direction=1)+xlab("Marker Gene")+ylab("Nuclei Type")

#remove cluster 16 as this is skeletal muscle
allcase<-subset(allcase,subset=RNA_snn_res.0.6!=20 & RNA_snn_res.0.6!= 16)

#redo above preprocessing
ident.list<-unique(allcase$orig.ident)
varlist<-c()
for (i in 1:length(ident.list)){
  temp<-subset(allcase,subset=orig.ident==ident.list[i])
  temp <- NormalizeData(temp, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
  temp <- FindVariableFeatures(temp, selection.method = "vst", nfeatures =500,verbose = F)
  varlist<-c(varlist,temp@assays$RNA@var.features)
}
varlist<-unique(varlist)
allcase@assays$RNA@var.features<-varlist

allcase <- NormalizeData(allcase, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
allcase<- ScaleData(allcase,verbose=F)
allcase<- RunPCA(allcase,verbose = F)

#run below section to visualize Number of PC to use to get best clustering results
for (i in c(9:20)){
  allcase <- allcase %>%FindNeighbors(dims=1:i,verbose=F)%>%FindClusters(resolution=0.6,verbose=F)%>%identity()
  allcase <- RunUMAP(allcase, dims = 1:i,verbose=F)
  DimPlot(allcase, label=TRUE, pt.size = .1,label.size = 3,repel=T)
  ggsave(paste0("/content/drive/MyDrive/KPMP/snRNA/snRNA012022/",i,"_nuclei_dim_.png"),plot = last_plot(),scale = 1, width = 5, height = 5, dpi = 600)
  feature=c("SLC5A12","SLC22A6","HAVCR1","VCAM1","TACSTD2",
            "UMOD","SLC12A1","SLC12A3","AQP2","SCNN1G",
            "SLC26A7","SLC4A9","CRB2","CLDN1","NPHS2",
            'ACTA2',"DCN","MYH11", "EMCN","FLT1","CD163","IL7R")
  DotPlot(allcase, features = feature,scale=T) + RotatedAxis()
  ggsave(paste0(i,"_nuclei_dot_.png"),plot = last_plot(),scale = 1, width = 8, height = 5, dpi = 600)
}

#seems PC=17 works well
allcase <- allcase %>%FindNeighbors(dims=1:17,verbose=F)%>%FindClusters(resolution=0.6,verbose=F)%>%identity()
allcase <- RunUMAP(allcase, dims = 1:17,verbose=F)

feature=c("SLC5A12","SLC22A6","HAVCR1","VCAM1","TACSTD2",
          "UMOD","SLC12A1","SLC12A3","AQP2","SCNN1G",
          "SLC26A7","SLC4A9","CRB2","CLDN1","NPHS2",
          'ACTA2',"DCN","MYH11", "EMCN","FLT1","CD163","IL7R")
DotPlot(allcase, features = feature,scale=T) + RotatedAxis()+scale_color_distiller(direction=1)+xlab("Marker Gene")+ylab("Nuclei Type")

cur.ident<-c(10,2,4,7,1,
             18,0,5,19,3,11,9,17,
             21,15,
             8,20,22,14,6,16,13,12
)
new.ident<-c("PT.1","PT.2","PT.3","PT.4","PT.5",
             "DTL_ATL","TAL","DCT.1","DCT.2","PC.1","PC.2","ICA","ICB",
             "PEC","POD",
             "Fib.1","Fib.2","Fib.3","Pericyte","EC.1","EC.2","Myeloid","Lymphocyte"
)
allcase$nuctype1<-plyr::mapvalues(allcase@active.ident,from=cur.ident,to=new.ident)
allcase$nuctype1<-factor(allcase$nuctype1,levels=new.ident)

#less refined cluster
cur.ident<-c(10,2,4,7,1,
             18,0,5,19,3,11,9,17,
             21,15,
             8,20,22,14,6,16,13,12
)
new.ident<-c("PT.1","PT.2","PT.3","PT.4","PT.5",
             "DTL_ATL","TAL","DCT","DCT","PC","PC","ICA","ICB",
             "PEC","POD",
             "Fib","Fib","Fib","Pericyte","EC","EC","Myeloid","Lymphocyte"
)
allcase$nuctype2<-plyr::mapvalues(allcase@active.ident,from=cur.ident,to=new.ident)
allcase$nuctype2<-factor(allcase$nuctype2,levels=c("PT.1","PT.2","PT.3","PT.4","PT.5",
                                                   "DTL_ATL","TAL","DCT","PC","ICA","ICB",
                                                   "PEC","POD",
                                                   "Fib","Pericyte","EC","Myeloid","Lymphocyte"
))
cur.ident<-c(10,2,4,7,1,
             18,0,5,19,3,11,9,17,
             21,15,
             8,20,22,14,6,16,13,12
)
new.ident<-c("PT","PT","PT","PT","PT",
             "DTL_ATL","TAL","DCT","DCT","PC","PC","ICA","ICB",
             "PEC","POD",
             "Fib","Fib","Fib","Pericyte","EC","EC","Myeloid","Lymphocyte"
)
allcase$nuctype3<-plyr::mapvalues(allcase@active.ident,from=cur.ident,to=new.ident)
allcase$nuctype3<-factor(allcase$nuctype3,levels=c("PT",
                                                   "DTL_ATL","TAL","DCT","PC","ICA","ICB",
                                                   "PEC","POD",
                                                   "Fib","Pericyte","EC","Myeloid","Lymphocyte"
))

#to get Figure 1A, 1B in the manuscript
DimPlot(allcase, label=TRUE, pt.size = .1,label.size = 3,repel=T,group.by="nuctype1")+scale_color_manual(value="Hiroshige",n=18)
DotPlot(allcase, features = feature,scale=T) + RotatedAxis()+scale_color_distiller(direction=1)+xlab("Marker Gene")+ylab("Nuclei Type")

#PT subcluster
pt<-subset(allcase,subset=nuctype3=="PT")

#Figure 1C. 
DotPlot(pt, features = c("SLC5A12","SLC22A6","PDZK1IP1","LRP2","SPP1","SOX4","CD24","HAVCR1","VCAM1"),scale=F) + RotatedAxis()+scale_color_distiller(direction=1)+xlab("Marker Gene")+ylab("Nuclei Type")


#to get Figure 1D in the manuscript (PT composition)
pt$nuctype1<-factor(pt$nuctype1,levels=c("PT.1","PT.2","PT.3","PT.4","PT.5"))
pt$nuctype1<-as.character(pt$nuctype1)
pt$orig.ident2<-as.character(pt$orig.ident2)
pt_table<-as.data.frame(cbind(pt$nuctype1,pt$orig.ident2))
colnames(pt_table)<-c("cluster","orig.ident2")
pt_table$count<-1
pt_table$orig.ident2<-as.factor(pt_table$orig.ident2)
pt_table$orig.ident2<-factor(pt_table$orig.ident2,levels=c("AKI #1", "AKI #2", "AKI #3", "AKI #4", "AKI #5", "AKI #6",
                                                           "Health Control #1", "Health Control #2","Health Control #3","Health Control #4",
                                                           "Health Control #5","Health Control #6","Health Control #7"))
ggplot(pt_table, aes(x = orig.ident2, y=count,fill = cluster)) +
  geom_bar(position = "fill",stat="identity")+
  xlab("Participant ID")+
  ylab("Proportion")+
  scale_fill_viridis(discrete = T,direction=-1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_rect(fill = "white", colour = "grey50"))
ggsave("pt_composition.png",plot = last_plot(),scale = 1, width = 6, height = 5, dpi = 600)


