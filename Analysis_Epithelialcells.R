
###降维聚类####
unique(OC01_KRT$seurat_cluster)
OC01_KRT<-subset(OC01,idents = c('Epithelial cells'))
OC01_KRT<-RunHarmony(OC01_KRT,group.by.vars = 'orig.ident')
OC01_KRT<-JackStraw(OC01_KRT,num.replicate = 100)
OC01_KRT<-ScoreJackStraw(OC01_KRT,dims=1:20)
JackStrawPlot(OC01_KRT,dims=1:20)
ElbowPlot(OC01_KRT,reduction = 'harmony')

OC01_KRT<-FindNeighbors(OC01_KRT,dims = 1:3,reduction = 'harmony')
OC01_KRT<-FindClusters(OC01_KRT,resolution = 0.05,reduction='harmony')
table(OC01_KRT@meta.data$seurat_clusters)
OC01_KRT<-RunUMAP(OC01_KRT,dims=1:3,reduction = 'harmony')
OC01_KRT<-RunTSNE(OC01_KRT,dims = 1:3,reduction = 'harmony')

DimPlot(OC01_KRT,reduction="umap",label = T)
DimPlot(OC01_KRT,reduction="tsne", label = T,split.by = 'group',cols = c('#C6307C','#4991C1'))
DimPlot(OC01_KRT,reduction="tsne",label = T,cols = c('#C6307C','#4991C1'))

DimPlot(OC01_KRT,reduction="tsne", label = T)
DotPlot(OC01_KRT,features =c('ANXA1'))


###寻找特异表达基因###
OC01_KRT.markers<-FindAllMarkers(OC01_KRT,only.pos=T,min.pct = 0.25,logfc.threshold = 0.25)
top10<-OC01_KRT.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
DoHeatmap(OC01_KRT,features = top10$gene)+scale_fill_gradientn(colors = c('#9ec1d4',"#fff7ac",'#e87651'))
write.csv(OC01_KRT.markers,'./NEW/cellmarker-Epi.csv')

#######差异表达分析####
OC01_KRT_group_DEG<-FindAllMarkers(OC01_KRT,group.by = 'group',only.pos=F,min.pct = 0.25,logfc.threshold = 0.25)
write.csv(OC01_KRT_group_DEG,'./NEW/DEG/DEG_group_EarlyvsLate.csv')
OC01_KRT_cluster_DEG<-FindAllMarkers(OC01_KRT,only.pos=F,min.pct = 0.25,logfc.threshold = 0.25)
write.csv(OC01_KRT_cluster_DEG,'./NEW/DEG/DEG_clusters_highvslow.csv')

ANXA1high_Epi<-subset(OC01_KRT,idents=c('ANXA1high_Epi'))
ANXA1high_Epi_group<-FindAllMarkers(ANXA1high_Epi,group.by = 'group',only.pos=F,min.pct = 0.25,logfc.threshold = 0.25)
write.csv(ANXA1high_Epi_group,'./NEW/DEG/DEG_ANXA1high_EarlyvsLate.csv')
ANXA1low_Epi<-subset(OC01_KRT,idents=c('ANXA1low_Epi'))
ANXA1low_Epi_group<-FindAllMarkers(ANXA1low_Epi,group.by = 'group',only.pos=F,min.pct = 0.25,logfc.threshold = 0.25)
write.csv(ANXA1low_Epi_group,'./NEW/DEG/DEG_ANXA1low_EarlyvsLate.csv')


###关注ANXA1的表达情况###
ANXA1exp<-FetchData(OC01_KRT,vars = 'ANXA1')
ANXA1exp$cell<-rownames(ANXA1exp)
cellgroup<-as.data.frame(OC01_KRT$group)
cellgroup$cell<-rownames(cellgroup)
celltype<-as.data.frame(OC01_KRT@active.ident)
celltype$cell<-rownames(celltype)
ANXA1exp=merge(ANXA1exp,celltype,by='cell',all=F)
ANXA1exp=merge(ANXA1exp,cellgroup,by='cell',all=F)
write.csv(ANXA1exp,'./NEW/ANXA1exp_Epithelial_cells_group_celltype.csv')

DotPlot(OC01_KRT,features=c('ANXA1'),cols=c('#fff7ac','#E6550D'))+RotatedAxis()

###注释，根据ANXA1表达高低区分

library(data.table)
new.cluster<-as.data.frame(fread('Epi_annot-2.csv'))
new.cluster<-new.cluster$`Cell type`

levels(OC01_KRT)
??RenameIdents
names(new.cluster)<-levels(OC01_KRT)
OC01_KRT<-RenameIdents(OC01_KRT,new.cluster)


cellnumber<-table(OC01_KRT$orig.ident,OC01_KRT@active.ident)
cellnumber<-as.data.frame(cellnumber)
write.csv(cellnumber,'./NEW/cellnumber_clusters.csv')
####cytoTRACE####
library(CytoTRACE)
expr_epi<-as.matrix(OC01_KRT[['RNA']]@counts)

expr_epi<-as.matrix(GetAssayData(OC01_KRT,assay = 'RNA',layer = 'counts'))

expr_epi<-expr_epi[apply(expr_epi>0,1,sum)>=5,]

expr_epi<as.matrix(expr_epi)
results<-CytoTRACE(expr_epi,ncores = 1)

phenot<-OC01_KRT@active.ident
phenot<-as.character(phenot)
names(phenot)<-rownames(OC01_KRT@meta.data)

emb<-OC01_KRT@reductions[['tsne']]@cell.embeddings
plotCytoTRACE(results,phenotype = phenot,emb = emb,outputDir = './CytoTRACE-2/')

###GSVA
library(GSVA)
library(data.table)
rm(OC_early)

gene.expr<-LayerData(OC01_KRT,assay = 'RNA',layer = 'data')
gene.expr<as.matrix(gene.expr)
dim(gene.expr)
head(gene.expr)
gene.expr<-gene.expr[which(rowSums(data) > 0),] 
is.na(gene.expr)
na.omit(gene.expr)
dim(gene.expr)
??GSVA::gsva

mygeneset<-as.data.frame(fread('Stemness_signature.csv'))
gsva.result <- GSVA::gsva(gene.expr, mygeneset,kcdf='Gaussian')
GSVA_results_2<-t(gsva.result)

GSVA_results_2<-cbind(GSVA_results_2,as.data.frame(OC01_KRT@active.ident))
GSVA_results_2<-cbind(GSVA_results_2,as.data.frame(OC01_KRT$orig.ident))
GSVA_results_2<-cbind(GSVA_results_2,as.data.frame(OC01_KRT$group))

write.csv(GSVA_results_2,'./NEW/GSVA_results_Epi_stemness_origident.csv')
####GSEA###
###比较早期晚期上皮癌细胞的差异
Epi_ANXA1high<-subset(OC01_KRT,idents=c('ANXA1high_Epi'))
Epi_ANXA1low<-subset(OC01_KRT,idents=c('ANXA1low_Epi'))
ANXA1high_Epi_DEG_group<-FindAllMarkers(Epi_ANXA1high,group.by = 'group',only.pos=F,min.pct = 0.25,logfc.threshold = 0.25)
ANXA1low_Epi_DEG_group<-FindAllMarkers(Epi_ANXA1low,group.by = 'group',only.pos=F,min.pct = 0.25,logfc.threshold = 0.25)
Epi_DEG_cluster<-FindAllMarkers(OC01_KRT,only.pos=F,min.pct = 0.25,logfc.threshold = 0.25)
Epi_DEG_group<-FindAllMarkers(OC01_KRT,group.by = 'group',only.pos=F,min.pct = 0.25,logfc.threshold = 0.25)
write.csv(Epi_DEG_group,'./NEW/Epi_DEG_group.csv')

####GSEA####
library(clusterProfiler)
library(data.table)
library(enrichplot)
BiocManager::install('clusterProfiler')
install.packages('clusterProfiler')
install.packages('GOSemSim',version='2.30.2')



diffgene<-as.data.frame(fread('./NEW/Epi_DEG_group_latevsearly.csv'))
rownames(diffgene)<-diffgene$V1
diffgene<-diffgene[,-1]


alldiif<-diffgene[order(diffgene$avg_log2FC,decreasing = T),]

gene_fc=alldiif$avg_log2FC

names(gene_fc)<-rownames(alldiif)

head(gene_fc)

immo<-clusterProfiler::read.gmt('c7.all.v2023.2.Hs.symbols.gmt')
immo<-clusterProfiler::read.gmt('c2.all.v2023.2.Hs.symbols.gmt')
immo<-clusterProfiler::read.gmt('h.all.v2023.2.Hs.symbols.gmt')
immo<-clusterProfiler::read.gmt('GALIE_TUMOR_STEMNESS_GENES.v2023.1.Hs.gmt')
h.all.v2023.2.Hs.symbols.gmt

immo$term<-sub('^[^_]*_','',immo$term)

gsea_1<-clusterProfiler::GSEA(gene_fc,TERM2GENE = immo,verbose = T,pvalueCutoff = 1)

g1<-as.data.frame(gsea_1)


g1<-subset(g1,qvalue<0.25)

write.csv(g1,'./NEW/Hallmark_Epi_DEG_late_vs_early.csv')







