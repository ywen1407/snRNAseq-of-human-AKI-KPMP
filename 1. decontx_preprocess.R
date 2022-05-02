library(scater)
library(ggplot2)
library(loomR)
library(Seurat)
library(dplyr)
library(patchwork)
library(SeuratData)
library(SeuratDisk)
library(tidyverse)
library(celda)
library(limma)

rm(list=ls())
gc()
memory.limit(size = 999999)

###read and deconx for AKI samples
setwd()
AKI.1_counts <-Read10X(data.dir = "~/30-10034")
results<-decontX(AKI.1_counts,delta=c(10,100),estimateDelta=F)
AKI_3010034 <- CreateSeuratObject(counts = results$decontXcounts, project = "30-10034", min.cells = 0, min.features = 0)
dim(AKI_3010034)

AKI.2_counts <-Read10X(data.dir = "~/32-2")
AKI_322 <- CreateSeuratObject(counts = AKI.2_counts, project = "32-2", min.cells = 0, min.features = 0)
dim(AKI_322)

AKI.3_counts <-Read10X(data.dir = "~/32-10003")
results<-decontX(AKI.3_counts,delta=c(10,75),estimateDelta=F)
AKI_3210003 <- CreateSeuratObject(counts = results$decontXcounts, project = "32-10003", min.cells = 0, min.features = 0)
#contamination indicated with CUBN and SLC8A1 in all clusters
dim(AKI_3210003)

AKI.4_counts <-Read10X(data.dir = "~/32-10034")
results<-decontX(AKI.4_counts,delta=c(10,50),estimateDelta=F)
AKI_3210034 <- CreateSeuratObject(counts = results$decontXcounts, project = "32-10034", min.cells = 0, min.features = 0)
dim(AKI_3210034)

AKI.5_counts <-Read10X(data.dir = "~/33-10005")
results<-decontX(AKI.5_counts,delta=c(10,100),estimateDelta=F)
AKI_3310005 <- CreateSeuratObject(counts = results$decontXcounts, project = "33-10005", min.cells = 0, min.features = 0)

AKI.6_counts <-Read10X(data.dir = "~/33-10006")
results<-decontX(AKI.6_counts,delta=c(10,10),estimateDelta=F)
AKI_3310006 <- CreateSeuratObject(counts = results$decontXcounts, project = "33-10006", min.cells = 0, min.features = 0)
#contamination by umod
dim(AKI_3310006)

####read and deconx for health control#############################################################################################
health.1_counts <-Read10X(data.dir = "~/18-142/")
results<-decontX(health.1_counts,delta=c(10,50),estimateDelta=F)
health_18142 <- CreateSeuratObject(counts = results$decontXcounts, project = "18-142", min.cells = 0, min.features = 0)
dim(health_18142)

health.4_counts <-Read10X(data.dir = "~/3535/")
results<-decontX(health.4_counts,delta=c(10,75),estimateDelta=F)
health_3535 <- CreateSeuratObject(counts = results$decontXcounts, project = "3535", min.cells = 0, min.features = 0)

health.5_counts <-Read10X(data.dir = "~/KRP446/")
results2<-decontX(health.5_counts,delta=c(10,10),estimateDelta=F)
health_KRP446 <- CreateSeuratObject(counts = results2$decontXcounts, project = "KRP446", min.cells = 0, min.features = 0)

health.6_counts <-Read10X(data.dir = "~/KRP460/")
health_KRP460 <- CreateSeuratObject(counts = health.6_counts, project = "KRP460", min.cells = 0, min.features = 0)

health.7_counts <-Read10X(data.dir = "~/KRP462/")
results<-decontX(health.7_counts,delta=c(10,25),estimateDelta=F)
health_KRP462 <- CreateSeuratObject(counts = results$decontXcounts, project = "KRP462", min.cells = 0, min.features = 0)

health.8_counts <-Read10X(data.dir = "~/3504/")
results<-decontX(health.8_counts,delta=c(10,125),estimateDelta=F)
health_3504 <- CreateSeuratObject(counts = results$decontXcounts, project = "3504", min.cells = 0, min.features = 0)

health.9_counts <-Read10X(data.dir = "~/KRP461/")
health_KRP461 <- CreateSeuratObject(counts = health.9_counts, project = "KRP461", min.cells = 0, min.features = 0)

####now merge all together
allcase <- merge(AKI_3010034,y=c(AKI_322,AKI_3210003,AKI_3210034,AKI_3310005,AKI_3310006,health_18142,
                                 health_3535,health_KRP446,health_KRP460,health_KRP462,health_3504,health_KRP461),
                 add.cell.ids=c("30-10034","32-2","32-10003","32-10034","33-10005","33-10006","18-142","3535","KRP446","KRP460","KRP462","3504","KRP461"),
                 project="AKI_Health")
allcase$orig.ident<-factor(x=allcase$orig.ident,levels = c("30-10034","32-2","32-10003","32-10034","33-10005","33-10006","18-142",
                                                           "3535","KRP446","KRP460","KRP462","3504","KRP461"))
allcase$case_type<-plyr::mapvalues(allcase$orig.ident,from = unique(allcase$orig.ident),
                                   to=c("AKI","AKI","AKI","AKI","AKI","AKI","Healthy_Control",
                                        "Healthy_Control","Healthy_Control","Healthy_Control","Healthy_Control","Healthy_Control","Healthy_Control"))

saveRDS(allcase,file="allcase_01042022_raw.rds")
