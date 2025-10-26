library(Seurat)
library(tidyverse)
library(Matrix)
library(spacexr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(readr)
library(config)
library(png)
library(patchwork)
library(SingleR)
library(celldex)
library(data.table)
library(devtools)
library(SeuratDisk)
library(hdf5r)
library(zellkonverter)
library(SingleCellExperiment)
library(reticulate)
OC.img=Seurat::Read10X_Image('./spatial',image.name = 'tissue_lowres_image.png',assay = 'Spatial')
OC=Seurat::Load10X_Spatial(data.dir = './',
                           filename = 'filtered_feature_bc_matrix.h5',
                           assay = 'Spatial',
                           slice = 'Slice1',
                           image=OC.img)
OC01<-readRDS(file='OC_16samples.rds')

SCearly<-subset(OC01,group==c('Early_OSCC'))
SClate<-subset(OC01,group==c('late_OSCC'))


OC01<-JoinLayers(OC01)

sc_counts<-GetAssayData(SCearly,assay='RNA',layer='counts')###以SCearly为例###

cellType<-as.data.frame(Idents(SCearly))

unique(Idents(SCearly))

cell_types = cellType$`Idents(SCearly)`; names(cell_types) <- rownames(cellType) # create cell_types named list
cell_types = as.factor(cell_types) # convert to factor data type
head(cell_types)

###修改空转对象中barcodes名使其与包含于coords中的barcode
SP_counts<-GetAssayData(OC,assay = 'Spatial',layer = 'counts')
###########
coords<-as.data.frame(fread('tissue_positions.csv'))
rownames(coords)<-coords$barcode
coords<-coords[,-1]
head(coords)
sc_nUMI = colSums(sc_counts)
sp_nUMI <- colSums(SP_counts)
########准备输入数据
reference = Reference(sc_counts, cell_types, sc_nUMI) #构建参考数据集
puck <- SpatialRNA(coords, SP_counts, sp_nUMI) #构建空间数据集
#执行RCTD
myRCTD <- create.RCTD(puck, reference)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet') #run.RCTD

#查找marker基因

get_marker_data <- function(cell_type_names, cell_type_means, gene_list) {
  
  marker_means = cell_type_means[gene_list,]
  
  marker_norm = marker_means / rowSums(marker_means)
  
  marker_data = as.data.frame(cell_type_names[max.col(marker_means)])
  
  marker_data$max_epr <- apply(cell_type_means[gene_list,],1,max)
  
  colnames(marker_data) = c('cell_type','max_epr')
  
  rownames(marker_data) = gene_list
  
  marker_data$log_fc <- 0
  
  epsilon <- 1e-9
  
  for(cell_type in unique(marker_data$cell_type)) {
    
    cur_genes <- gene_list[marker_data$cell_type == cell_type]
    
    other_mean = rowMeans(cell_type_means[cur_genes,cell_type_names != cell_type])
    
    marker_data$log_fc[marker_data$cell_type == cell_type] <- log(epsilon + cell_type_means[cur_genes,cell_type]) - log(epsilon + other_mean)
    
  }
  
  return(marker_data)
  
}

cell_type_info_restr = myRCTD@cell_type_info$info

de_genes <- get_de_genes(cell_type_info_restr, puck, fc_thresh = 3, expr_thresh = .0001, MIN_OBS = 3)

marker_data_de = get_marker_data(cell_type_info_restr[[2]], cell_type_info_restr[[1]], de_genes)


#构建绘图数据

results <- myRCTD@results

results_df <- results$results_df

barcodes = rownames(results_df[results_df$spot_class != 'reject' & puck@nUMI >= 1,])

my_table = puck@coords[barcodes,]

my_table$class = results_df[barcodes,]$first_type

#绘图

cal_zoom_rate <- function(width, height){
  
  std_width = 1000
  
  std_height = std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
  
  if (std_width / std_height > width / height){
    
    scale = width / std_width
    
  }
  
  else{
    
    scale = height / std_height
    
  }
  
  return(scale)
  
}

png <- readPNG('tissue_lowres_image.png')

zoom_scale = cal_zoom_rate(dim(png)[2], dim(png)[1])

my_table = my_table %>% mutate(across(c(x,y), ~.x*zoom_scale))

col = c('#F56867','#FEB915','#C798EE','#9BE86E','#7495D3','#D1D1D1','#6D1A9C','#15821E','#3A84E6','#997273','#787878','#DB4C6C','#9E7A7A','#554236','#AF5F3C','#93796C','#F9BD3F','#DAB370')

p = ggplot(my_table, aes(x = x, y = dim(png)[1] - y))+
  background_image(png)+
  geom_point(shape = 16, size = 1.8, aes(color = class))+
  
  coord_cartesian(xlim = c(0, dim(png)[2]), y = c(0, dim(png)[1]), expand = FALSE)+
  
  scale_color_manual(values = col)+
  
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  
  guides(color = guide_legend(override.aes = list(size = 2.5, alphe = 0.1)))

p = ggplot(my_table, aes(x = x, y = dim(png)[1] - y))+
  
  background_image(png)+
  
  geom_point(shape = 16, size = 0.1, aes(color = class))+
  
  coord_cartesian(xlim = c(-600, dim(png)[2]), y = c(-100, dim(png)[1]), expand = FALSE)+
  
  scale_color_manual(values = col)+
  
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  
  guides(color = guide_legend(override.aes = list(size = 5, alphe = 0.1)))

ggsave(p, file = 'spatial_all.png', width=10, height=7, dpi = 300)

###绘制概率分布图####
unique(myRCTD@cell_type_info)
barcodes <- colnames(myRCTD@spatialRNA@counts)
weights <- myRCTD@results$weights
norm_weights <- normalize_weights(weights)         
head(norm_weights)
levels(OC01)

p <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'Myeloid cells'],ylimit = c(0,0.2) , 
                          title ='plot of Myeloid cells', size=1, alpha=0.8) 
p <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'CD8 T cells'],ylimit = c(0,0.05) , 
                          title ='plot of CD8 T cells', size=2, alpha=0.8) 
p <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'CD4 T cells'],ylimit = c(0,0.2) , 
                          title ='plot of CD4 T cells', size=2, alpha=0.8)
p <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'Fibroblasts'],ylimit = c(0,0.2) , 
                          title ='plot of Fibroblasts', size=2, alpha=0.8)
p <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'Endothelial cells'],ylimit = c(0,0.2) , 
                          title ='plot of Endothelial cells', size=2, alpha=0.8)
p <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'Mast cells'],ylimit = c(0,0.2) , 
                          title ='plot of Mast cells', size=2, alpha=0.8)
p <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'Pericytes'],ylimit = c(0,0.2) , 
                          title ='plot of Pericytes', size=2, alpha=0.8)
p <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'B cells'],ylimit = c(0,0.2) , 
                          title ='plot of B cells', size=2, alpha=0.8)
p <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'Nerve-like cells'],ylimit = c(0,0.2) , 
                          title ='plot of Nerve-like cells', size=2, alpha=0.8)
p <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,'Epithelial cells'],ylimit = c(0,1) , 
                          title ='plot of Epithelial cells', size=1, alpha=0.8)

#############绘制饼图#########
devtools::install_github('JEFworks-Lab/STdeconvolve')
library(STdeconvolve)
library(ggplot2)
library(ggsci)
barcodes <- colnames(myRCTD@spatialRNA@counts)
weights <- myRCTD@results$weights
norm_weights <- normalize_weights(weights) 
clusterCols <- c("#843C39", "#8C6D31", "#E6550D", "#3182BD", "#54990F","#BD9E39", "#E7BA52", "#31A354", "#E41A1C")
clusterCols <- c("#843C39","#31A354", "#3182BD", "#E6550D", "#E7BA52","#8C6D31", "#E41A1C","#BD9E39", "#54990F",'black')
clusterCols <- c("#843C39", "#8C6D31", "#E6550D", "#3182BD", "#54990F","#BD9E39", "#E7BA52", "#31A354", "#E41A1C",'#DA70D6')

m <- as.matrix(norm_weights)
p <- coords
colnames(p)<-c('x','y')

plt <- vizAllTopics(theta = m,
                    pos = p,
                    topicOrder=seq(ncol(m)),
                    topicCols=clusterCols,
                    groups = NA,
                    group_cols = NA,
                    r = 0.7, # size of scatterpies; adjust depending on the coordinates of the pixels
                    lwd = 0.01,
                    showLegend = TRUE,
                    plotTitle = "scatterpies")

plt <- plt + ggplot2::guides(fill=ggplot2::guide_legend(ncol=2))
plt

ggsave(paste0(savedir, "Spaital_scatterpies.png"), width=12, height=6, plot=plt, bg="white")



