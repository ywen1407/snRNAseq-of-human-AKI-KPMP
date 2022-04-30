library(Seurat,verbose=F)
library(dplyr,verbose=F)
library(SeuratDisk,verbose=F)
library(celda,verbose=F)
library(ggplot2)
library(sctransform)


memory.limit(99999999999999999)

setwd()
count_3210205<-Read10X(data.dir = "32-10205") #no contamination
AKI_32_10205<- CreateSeuratObject(counts = count_3210205, project = "32-10205")

count_3410050<-Read10X(data.dir = "34-10050") 
AKI_34_10050<- CreateSeuratObject(counts = count_3410050, project = "34-10050")

count_3410187<-Read10X(data.dir = "34-10187") #small contamination by UMOD, will use without decont for now.
AKI_34_10187<- CreateSeuratObject(counts = count_3410187, project = "34-10187")

count_3410209<-Read10X(data.dir = "34-10209") 
AKI_34_10209<- CreateSeuratObject(counts = count_3410209, project = "34-10209")

count_3010044<-Read10X(data.dir = "30-10044") #contaminated by UMOD, no PT, will not use for now. 
count_3410184<-Read10X(data.dir = "34-10184") #full contamination with UMOD, unable to fully eliminate, will not use for now. 

additional<-merge(AKI_32_10205,y=c(AKI_34_10050,AKI_34_10187,AKI_34_10209),
                  add.cell.ids = c("32-10205","34-10050","34-10187","34-10209"))
additional[["percent.mt"]] <- PercentageFeatureSet(additional, pattern = "^MT-")
additional<-subset(additional,subset=nFeature_RNA>=500&nFeature_RNA<=5000&nCount_RNA>500&nCount_RNA<10000&percent.mt<=2)
dim(additional)
saveRDS(additional,"additional_4case_02152022.RDS")

#read previously precessed non-covid AKI, healthy reference file
noncovid<-readRDS("allcase_01032022_alt_processed.RDS")
noncovid<-subset(noncovid,subset=orig.ident!="18-142") #remove this healthy reference as significant batch effect noted from previous analysis
dim(noncovid)

#merge previously processed non-covid AKI+ reference with additional AKI cases
noncovid2<-merge(noncovid,y=additional)
saveRDS(noncovid2,"noncovid_AKI_02152022_filtered.RDS")

#read COVID AKI cases shared by Blue Lake. 
covid<-LoadH5Seurat("COVID_bluelake.h5Seurat")
dim(covid)
covid<-subset(covid,subset=nFeature_RNA>=500&nFeature_RNA<=5000&nCount_RNA>500&nCount_RNA<10000&percent.mt<=2)
dim(covid)
colnames(covid@meta.data)
saveRDS(covid,"covid_filtered.RDS")


#read previously merged non-covid AKI+ reference cases. 
#then merge non-covid and covid together
setwd("//homer.win.ad.jhu.edu/Parikh Lab/Fellows/Yumeng/KPMP_snRNA")

noncovid<-readRDS("noncovid_AKI_02152022_filtered.RDS")
covid<-readRDS("covid_filtered.RDS")
allcase<-merge(noncovid,y=covid)
saveRDS(allcase,"allcase_combined_02152022_filtered.RDS")


allcase[["percent.mt"]] <- PercentageFeatureSet(allcase, pattern = "^MT-")
dim(allcase)
allcase<-subset(allcase,subset=nFeature_RNA>=500&nFeature_RNA<=5000&nCount_RNA>500&nCount_RNA<10000&percent.mt<=2)
dim(allcase)
allcase <- NormalizeData(allcase, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
allcase <- FindVariableFeatures(allcase, selection.method = "vst", nfeatures =1800,verbose = F)
allcase<- ScaleData(allcase,verbose=F)

allcase<- RunPCA(allcase,verbose = F)
allcase <- allcase %>%FindNeighbors(dims=1:15,verbose=F)%>%FindClusters(resolution=0.6,verbose=F)%>%identity()
allcase <- RunUMAP(allcase, dims = 1:15,verbose=F)
feature=c("SLC5A12","SLC22A6","HAVCR1","LCN2","VCAM1","TACSTD2",
          "UMOD","SLC12A1","SLC12A3","AQP2","SCNN1G",
          "SLC26A7","SLC4A9","CRB2","CLDN1","NPHS1","NPHS2","EMCN","FLT1",
          'ACTA2',"DCN","MYH11","C1QC","S100A8","IL7R","NKG7","MS4A1","MZB1")
DotPlot(allcase, features = feature,scale=F) + RotatedAxis()
RidgePlot(allcase,features = "UMOD")
