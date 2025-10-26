
###Infercnv###
####安装环境####
library(data.table)
library(infercnv)
library(readr)
##############

OC01_KRT_CNV <- subset(OC01,idents = c('Epithelial cells','CD8 T cells','CD4 T cells','Pericytes'))

unique(Idents(OC01))
gene_order<-read_tsv('hg38_gencode_v27.txt',col_names=c('gene','chr','start','end'))
gene.expr <-  GetAssayData(OC01,assay='RNA',layer = 'data')
anno<-as.data.frame(Idents(OC01))
OC01<-JoinLayers(OC01)
unique(Idents(OC01))

gene.expr<-gene.expr[rownames(gene.expr)%in%gene_order$gene,]

infercnv_obj=CreateInfercnvObject(raw_counts_matrix = gene.expr,
                                  annotations_file = anno,
                                  delim = '\t',
                                  gene_order_file = 'hg38_gencode_v27.txt',
                                  ref_group_names = c('T cells','B cells'))

infercnv_obj=infercnv::run(infercnv_obj,cutoff = 0.1,
                           out_dir = 'InferCNV/',
                           cluster_by_groups=T,
                           denoise=F,
                           write_expr_matrix=T,
                           HMM=F,
                           min_cells_per_gene = 10,
)
save(infercnv_obj,file='./inferCNV/infercnv_obj.Rdata')

data<-read.table('./inferCNV/infercnv.observations.txt',header = T)
expr=data%>%as.matrix()
expr.scale=scale(t(expr))
tmp1=sweep(expr.scale,2,apply(expr.scale,2,min),'-')
tmp2=apply(expr.scale,2,max) - apply(expr.scale,2,min)
expr_1=t(2*sweep(tmp1,2,tmp2,"/")-1)

cnv_score=as.data.frame(colSums(expr_1*expr_1))

colnames(cnv_score)='cnv_score'
cnv_score=rownames_to_column(cnv_score,var = 'cell')
gene_dat

gene_data2$cell<-rownames(gene_data2)
colnames(gene_data2)[1]='cluster'
test=merge(cnv_score,gene_data2,by='cell',all=F)

test2=merge(test,celltype,by='cell',all=F)

test3=merge(test2,cellidents,by='cell',all=F)

cellidents<-as.data.frame(OC01$orig.ident)

cellidents$cell<-rownames(cellidents)

ggplot2::ggplot(test,aes(x=cluster,y=cnv_score))+
  geom_violin(aes(fill=cluster),cex=1.2)+
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(color='black'),
        legend.position = 'none')
write.csv(test3,'cnv_score.csv')
