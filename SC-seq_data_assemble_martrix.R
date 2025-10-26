library(installr)
library(infercnv)
library(mindr)
library(multtest)
library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(harmony)
library(devtools)
##安装所需的包，需注意mindr已被CRAN移除，可手动下载压缩包手动加载##

setwd('D:/SCseq/Matrix/')

rm(list = ls())
unique(OC01$orig.ident)
##自动读取10X数据，三个文件##

PT1.data<-Read10X(data.dir='PT1')
PT2.data<-Read10X(data.dir='PT2')
PT3.data<-Read10X(data.dir='PT3')
PT4.data<-Read10X(data.dir='PT4')
PT5.data<-Read10X(data.dir='PT5')
PT6.data<-Read10X(data.dir='PT6')
PT7.data<-Read10X(data.dir='PT7')
PT8.data<-Read10X(data.dir='PT8')
PT9.data<-Read10X(data.dir='PT9')
PT10.data<-Read10X(data.dir='PT10')
PT11.data<-Read10X(data.dir='PT11')
PT12.data<-Read10X(data.dir='PT12')
PT13.data<-Read10X(data.dir='PT13')
PT14.data<-Read10X(data.dir='PT14')
PT15.data<-Read10X(data.dir='PT15')
PT16.data<-Read10X(data.dir='PT16')

##创建seurat对象##

PT1<-CreateSeuratObject(counts=PT1.data,project='PT1')
PT2<-CreateSeuratObject(counts=PT2.data,project='PT2')
PT3<-CreateSeuratObject(counts=PT3.data,project='PT3')
PT4<-CreateSeuratObject(counts=PT4.data,project='PT4')
PT5<-CreateSeuratObject(counts=PT5.data,project='PT5')
PT6<-CreateSeuratObject(counts=PT6.data,project='PT6')
PT7<-CreateSeuratObject(counts=PT7.data,project='PT7')
PT8<-CreateSeuratObject(counts=PT8.data,project='PT8')
PT9<-CreateSeuratObject(counts=PT9.data,project='PT9')
PT10<-CreateSeuratObject(counts=PT10.data,project='PT10')
PT11<-CreateSeuratObject(counts=PT11.data,project='PT11')
PT12<-CreateSeuratObject(counts=PT12.data,project='PT12')
PT13<-CreateSeuratObject(counts=PT13.data,project='PT13')
PT14<-CreateSeuratObject(counts=PT14.data,project='PT14')
PT15<-CreateSeuratObject(counts=PT15.data,project='PT15')
PT16<-CreateSeuratObject(counts=PT16.data,project='PT16')

###merge###

OC01list=list(PT1,PT2,PT3,PT4,PT5,PT6,PT7,PT8,PT9,P10,PT11,PT12,PT13,PT14,PT15,PT16)
OC01<-merge(x=OC01list[[1]],y=OC01list[-1])
unique(Idents(OC01))
saveRDS(OC01,file='OC_16samples.rds') ####保存一下
OC01<-readRDS(file='OC_16samples.rds')####读取数据

all.genes=rownames(OC01)
all.genes
OC01<-JoinLayers(OC01)
OC01<-NormalizeData(OC01) %>% FindVariableFeatures() %>% ScaleData(features = all.genes) %>% RunPCA(verbose = T)
OC01<-RunHarmony(OC01,group.by.vars = 'orig.ident')
unique(OC01$orig.ident)

##计算线粒体RNA含量##
OC01[["percent.mt"]]<-PercentageFeatureSet(OC01,pattern = "^MT-")
head(OC01@meta.data,5)
VlnPlot(OC01,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)

install.packages("patchwork")


plot1<-FeatureScatter(OC01,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2<-FeatureScatter(OC01,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
CombinePlots(plots=list(plot1,plot2))

###初步质控，过滤RNA过多或过少的细胞##

OC01<-subset(OC01,subset=nFeature_RNA>200&nFeature_RNA<2500&percent.mt<10)

ncol(OC01.data)

ncol(as.data.frame((OC01[["RNA"]]@counts)))

###寻找高频变化基因，FOR PCA###

top10<-head(VariableFeatures(OC01),10)

plot1<-VariableFeaturePlot(OC01)
plot2<-LabelPoints(plot=plot1,points = top10,repel=TRUE)
plot1
plot2

OC01<-JackStraw(OC01,num.replicate = 100)
OC01<-ScoreJackStraw(OC01,dims=1:20)
JackStrawPlot(OC01,dims=1:20)
ElbowPlot(OC01,reduction = 'harmony')

OC01<-FindNeighbors(OC01,dims = 1:15,reduction = 'harmony')
OC01<-FindClusters(OC01,resolution = 0.2,reduction='harmony')
table(OC01@meta.data$seurat_clusters)
OC01<-RunUMAP(OC01,dims=1:11,reduction = 'harmony')
DimPlot(OC01,reduction="umap",label = T,label.size = 5)
DimPlot(OC01,reduction="umap",label = T,split.by = 'orig.ident',label.size = 5)
DimPlot(OC01,reduction="umap",label = F,label.size = 5,group.by = 'orig.ident')
OC01<-RunTSNE(OC01,dims = 1:11,reduction = 'harmony')
DimPlot(OC01,reduction="tsne", label = TRUE)

OC01<-RunPCA(OC01,features = VariableFeatures(object = OC01))
print(OC01[["pca"]],dims=1:5,nfeatures=5)
VizDimLoadings(OC01,dims=1:2,reduction = "pca")
DimPlot(OC01,reduction="pca")

DimHeatmap(OC01,dims=1:12,cells=500,balanced=TRUE)

####注释细胞类型#####
library(data.table)
new.cluster<-as.data.frame(fread('annot.csv'))
new.cluster<-new.cluster$`Cell type`

levels(OC01)

names(new.cluster)<-levels(OC01)
OC01<-RenameIdents(OC01,new.cluster)

DotPlot(OC01,features = c('CD3D','CD4','TPSAB1','KIT','CNTNAP2','CD19','CD79A','TAGLN','RGS5','PECAM1','VWF','DCN','COL1A1','LYZ','CD14','KRT14','KRT5','CD8A'),
        col.min = 0.5,cols=c('#fff7ac','#E6550D'))+RotatedAxis()###Dotplot_celltypes_all


cluster5.markers<-FindMarkers(OC01,ident.1=5,ident.2 = c(0,3),min.pct = 0.25)

OC01.markers<-FindAllMarkers(OC01,only.pos=TRUE,min.pct = 0.25,logfc.threshold = 0.25)
OC01.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
top10_2<-OC01.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
DoHeatmap(OC01,features = top10_2$gene,)+scale_fill_gradientn(colors = c('#9ec1d4',"#fff7ac",'#e87651'))

###将不同来源样本(orig.ident)进行分组###
all_groups<-unique(OC01$orig.ident)

group_to_merge<-c('PT8','PT9','PT10','PT11','PT12','PT13','PT14','PT15','PT16')
merged_group<-'late_OSCC'

group_to_merge2<-c('PT1','PT2','PT3','PT4','PT5','PT6','PT7')
merged_group2<-'Early_OSCC'

OC01$group<-OC01$orig.ident

unique(OC01$group)

OC01$group<-ifelse(OC01$group %in% group_to_merge2, merged_group2,OC01$group)
OC01$group<-ifelse(OC01$group %in% group_to_merge, merged_group,OC01$group)

cols<-c('#2b83ba','#abdda4')

DimPlot(OC01,reduction="umap",label = T,cols=cols,group.by = 'group')

clusterCols <- c("#843C39", "#8C6D31", "#E6550D", "#3182BD", "#54990F","#BD9E39", "#E7BA52", "#31A354", "#E41A1C",'#DA70D6')

DimPlot(OC01,reduction="umap",label = T,cols=clusterCols)


