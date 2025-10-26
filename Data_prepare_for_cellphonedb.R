####提取细胞类型矩阵
###以PT1为例###
unique(OC01$orig.ident)
PT1<-subset(OC01, orig.ident == 'PT1')
celltype<-as.data.frame(PT1@active.ident)
name<-rownames(celltype)
rownames(celltype)<-NULL
celltype<-cbind(name,celltype)
colnames(celltype)<-c('Cell','cell_type')
write.table(celltype, file="./cellphonedb/separate/CellphoneDB_celltype_PT1.txt", row.names=T, col.names=NA,sep='\t',quote = F)

####提取细胞表达矩阵

expr_S<-GetAssayData(PT1,assay='RNA',layer='counts')
count_norm_S <- apply(expr_S, 2, function(x) (x/sum(x))*10000)
write.table(count_norm_S, file="./cellphonedb/separate/CellphoneDB_expr_normalized_PT1.txt", row.names=T, col.names=NA,sep='\t',quote = F)


