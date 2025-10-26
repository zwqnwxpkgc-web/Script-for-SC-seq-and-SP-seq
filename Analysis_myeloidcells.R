library(installr)
library(infercnv)
library(mindr)
library(multtest)
library(Seurat)
library(dplyr)
library(tidyverse)
library(harmony)
library(devtools)

OC01_ml <- subset(OC01,idents = c('Myeloid cells'))
OC01_ml<-ScaleData(OC01_ml,features=rownames(OC01_ml))
OC01_ml<-RunHarmony(OC01_ml,group.by.vars = 'orig.ident')

OC01_ml<-JackStraw(OC01_ml,num.replicate = 200)
OC01_ml<-ScoreJackStraw(OC01_ml,dims=1:15)
JackStrawPlot(OC01_ml,dims=1:15)
ElbowPlot(OC01_ml,reduction = 'harmony')

OC01_ml<-FindNeighbors(OC01_ml,dims = 1:4,reduction = 'harmony')
OC01_ml<-FindClusters(OC01_ml,resolution = 0.1,reduction='harmony')
table(OC01_ml@meta.data$seurat_clusters)
unique(OC01_ml$seurat_clusters)
OC01_ml<-RunUMAP(OC01_ml,dims=1:4,reduction = 'harmony')
OC01_ml<-RunTSNE(OC01_ml,dims=1:4,reduction='harmony')
DimPlot(OC01_ml,reduction="umap")#,split.by = 'orig.ident')
DimPlot(OC01_ml,reduction='tsne')


cellnumber<-table(OC01_ml$orig.ident,OC01_ml$seurat_clusters)
cellnumber<-as.data.frame(cellnumber)
write.csv(cellnumber,'myeloid_number.csv')
genes_to_check = c('CD68','HLA-DRA','HLA-DRB1','ITGAM','SELL','CXCR4','FCGR3B','ARG1','ITGAE','THBD','CD1A','CD1C','MMP9')
DotPlot(OC01_ml,features = unique(genes_to_check ),scale.max=50,cols=c('#fff7ac','#E6550D')) + RotatedAxis()
DimPlot(OC01_ml,reduction = 'umap',label = T)
OC01_ml<-RenameIdents(OC01_ml,'0'='MDSCs','1'='TAM','2'='TAM')
saveRDS(OC01_ml,'./Myeloid/Myeloidcells.rds')

###注释###

new.cluster<-as.data.frame(fread('./Myeloid/ml_annot.csv'))
new.cluster<-new.cluster$`Cell type`

levels(OC01_ml)
unique(OC01_ml$seurat_clusters)
names(new.cluster)<-unique(OC01_ml$seurat_clusters)
OC01_ml<-RenameIdents(OC01_ml,new.cluster)
OC01_ml.markers<-FindAllMarkers(OC01_ml,only.pos=TRUE,min.pct = 0.25,logfc.threshold = 0.25)
write.csv(OC01_ml.markers,'OC01_ml_markers.csv')

top10_3<-OC01_ml.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
DoHeatmap(OC01_ml,features = top10_3$gene)

DotPlot(OC01_ml, features = c('FPR1','FPR2','FPR3'),scale.max = 20,cols=c('#fff7ac','#E6550D')) + RotatedAxis()

####以FPR1/2/3表达值高低进行分群分析####
FPRexp<-FetchData(OC01_ml,vars =c('FPR1','FPR2','FPR3'))
mean_FPR1<-mean(FPRexp$FPR1)
mean_FPR2<-mean(FPRexp$FPR2)
mean_FPR3<-mean(FPRexp$FPR3)
OC01_ml_FPR1<-subset(OC01_ml,FPR1>mean_FPR1)
OC01_ml_FPR2<-subset(OC01_ml,FPR2>mean_FPR2)
OC01_ml_FPR2_low<-subset(OC01_ml,FPR2<mean_FPR2)
OC01_ml_FPR3<-subset(OC01_ml,FPR3>mean_FPR3)
saveRDS(OC01_ml_FPR1,file='FPR1high_myeloidcells.rds')
saveRDS(OC01_ml_FPR2,file='FPR2high_myeloidcells.rds')
saveRDS(OC01_ml_FPR3,file='FPR3high_myeloidcells.rds')
OC01_ml_FPR1<-readRDS(file='FPR1high_myeloidcells.rds')
OC01_ml_FPR2<-readRDS(file='FPR2high_myeloidcells.rds')
OC01_ml_FPR3<-readRDS(file='FPR3high_myeloidcells.rds')


cellnumber<-table(Idents(OC01_ml_FPR2),OC01_ml_FPR2$orig.ident)
cellnumber<-table(Idents(OC01_ml_FPR1),OC01_ml_FPR1$orig.ident)
cellnumber<-table(Idents(OC01_ml_FPR3),OC01_ml_FPR3$orig.ident)
cellnumber<-table(Idents(OC01_ml_FPR2_low),OC01_ml_FPR2_low$orig.ident)
cellnumber<-table(Idents(OC01_ml),OC01_ml$orig.ident)
write.csv(cellnumber,'FPR2high_myeloid_cellnumber.CSV')
write.csv(cellnumber,'FPR1high_myeloid_cellnumber.CSV')
write.csv(cellnumber,'FPR3high_myeloid_cellnumber.CSV')
write.csv(cellnumber,'FPR2low_myeloid_cellnumber.CSV')
write.csv(cellnumber,'myeloid_cellnumber.CSV')
DimPlot(OC01_ml_FPR2,reduction = 'umap',label = T)

FPRexp$cell<-rownames(FPRexp)
cellgroup<-as.data.frame(OC01_ml$group)
cellgroup$cell<-rownames(cellgroup)
celltype<-as.data.frame(OC01_ml@active.ident)
celltype$cell<-rownames(celltype)
FPRexp=merge(FPRexp,celltype,by='cell',all=F)
FPRexp=merge(FPRexp,cellgroup,by='cell',all=F)
write.csv(FPRexp,'FPR123exp_myeloid_cells_group_celltype.csv')


DotPlot(OC01_ml_FPR2,features=c('CD68','HLA-DRA','HLA-DRB1','ITGAM','SELL','CXCR4','FCGR3B',
                                'ARG1','ITGAE','THBD','CD1A','CD1C','MMP9'),scale = T,scale.max = 50,cols=c('#fff7ac','#E6550D'))+RotatedAxis()
DotPlot(OC01_ml_2,features=c('FPR1','FPR2','FPR3'),scale = T)



###注释

new.cluster<-c('MDSCs','TAM')
new.cluster<-new.cluster$`Cell type`
levels(OC01_ml)
names(new.cluster)<-levels(OC01_ml)
OC01_ml<-RenameIdents(OC01_ml,new.cluster)
cellnumber<-table(OC01_ml$group,OC01_ml@active.ident)

OC01_ml_.markers<-FindAllMarkers(OC01_ml,only.pos=TRUE,min.pct = 0.25,logfc.threshold = 0.25)
top10_ml<-OC01_ml_.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
DoHeatmap(OC01_ml,features = top10_ml$gene)
DoHeatmap(OC01_ml,features = top10_ml$gene)+scale_fill_gradientn(colors = c('#9ec1d4',"#fff7ac",'#e87651'))


