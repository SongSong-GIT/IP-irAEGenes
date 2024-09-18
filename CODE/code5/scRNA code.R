
#####################################################################
# 加载必要的库（如果尚未安装，请先安装）
setwd("E:\\iraes-pathway工作总结\\图\\单细胞\\SKCM")
# 4 download from URL
packageurl <- "https://cran.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_4.1.4.tar.gz" 
install.packages(packageurl, repos=NULL, type="source")
packageurl <- "https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_4.4.0.tar.gz" 
install.packages(packageurl, repos=NULL, type="source")

library(Seurat)
#ls("package:DoubletFinder")   查看包里面的函数
# 解压并加载原始计数数据
count_data <- read.table(gzfile("E:\\软件\\GSE72056_melanoma_single_cell_revised_v2.txt.gz"), header = TRUE)
count_data=count_data[-c(1:3),]
exp=aggregate(.~Cell,mean,data=count_data)
head(exp)
rownames(exp) <- exp$Cell
exp <- exp[,-1]

# 创建 Seurat 对象
sce <- CreateSeuratObject(counts =exp)
#####################################
mito_genes=rownames(sce)[grep("^MT-", rownames(sce))] 
mito_genes #13个线粒体基因
sce=PercentageFeatureSet(sce, "^MT-", col.name = "percent_mito")
fivenum(sce@meta.data$percent_mito)
#计算核糖体基因比例
ribo_genes=rownames(sce)[grep("^Rp[sl]", rownames(sce),ignore.case = T)]
ribo_genes
sce.all=PercentageFeatureSet(sce.all, "^RP[SL]", col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)
#计算红血细胞基因比例
rownames(sce)[grep("^Hb[^(p)]", rownames(sce),ignore.case = T)]
sce=PercentageFeatureSet(sce, "^HB[^(P)]", col.name = "percent_hb")
fivenum(sce@meta.data$percent_hb)
#可视化细胞的上述比例情况
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_hb")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce, group.by = "orig.ident", features = feats, pt.size = 0.01, ncol = 2) + 
  NoLegend()
p1
library(ggplot2)
#ggsave(filename="Vlnplot1.pdf",plot=p1)
ggsave(filename="Vlnplot1.png",plot=p1)
feats <- c("percent_mito", "percent_hb")
p2=VlnPlot(sce, group.by = "orig.ident", features = feats, pt.size = 0.01, ncol = 3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2	
# ggsave(filename="Vlnplot2.pdf",plot=p2)
ggsave(filename="Vlnplot2.png",plot=p2)

p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
ggsave(filename="Scatterplot.pdf",plot=p3)
#根据上述指标，过滤低质量细胞/基因
#过滤指标1:最少表达基因数的细胞&最少表达细胞数的基因
selected_c <- WhichCells(sce, expression = nFeature_RNA <5000)
selected_f <- rownames(sce)[Matrix::rowSums(sce@assays$RNA@counts >0) > 3]

sce.all.filt <- subset(sce, features = selected_f, cells = selected_c)

dim(sce) 
dim(sce.all.filt) 
#  可以看到，主要是过滤了基因，其次才是细胞

#过滤指标2:线粒体/核糖体基因比例(根据上面的violin图)
dim(sce.all.filt)
selected_mito <- WhichCells(sce.all.filt, expression = percent_mito < 15)
selected_hb <- WhichCells(sce.all.filt, expression = percent_hb < 0.05)
length(selected_hb)
length(selected_mito)


sce.all.filt <- subset(sce.all.filt, cells = selected_mito)
sce.all.filt <- subset(sce.all.filt, cells = selected_hb)
dim(sce.all.filt)

table(sce.all.filt$orig.ident) 
#可视化过滤后的情况
feats <- c("nFeature_RNA", "nCount_RNA")
p1_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 2) + 
  NoLegend()
ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered)

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + 
  NoLegend()
ggsave(filename="Vlnplot2_filtered.pdf",plot=p2_filtered)

#过滤指标3:过滤特定基因
# Filter MALAT1 管家基因
sce.all.filt <- sce.all.filt[!grepl("MALAT1", rownames(sce.all.filt),ignore.case = T), ]
# Filter Mitocondrial 线粒体基因
sce.all.filt <- sce.all.filt[!grepl("^MT-", rownames(sce.all.filt),ignore.case = T), ]
# 当然，还可以过滤更多
da=readRDS("E:\\iraes-pathway工作总结\\图\\单细胞\\SKCM\\GSE244983\\da.rds")
dim(sce.all.filt)
saveRDS(sce.all.filt,"da.rds")
######################################
da=sce.all.filt
da <- NormalizeData(da, normalization.method = "LogNormalize")
da <- FindVariableFeatures(da, selection.method = "vst", nfeatures = 2000)
da <- ScaleData(da)
da <- RunPCA(da, npcs = 20, features = VariableFeatures(da),verbose = FALSE)
da <- FindNeighbors(da, reduction = "pca", dims = 1:20)
da <- FindClusters(da, resolution = 0.3) 
#da <- FindClusters(da, resolution = 1)####### 可以自己调整
da <- RunUMAP(da, reduction = "pca", dims = 1:20)
DimPlot(da,group.by = "RNA_snn_res.0.3",label = T) 
#########################################################
cellAnnotation <- function(obj, markerList, assay = 'SCT', slot = 'data') {
  markergene <- unique(do.call(c, markerList))
  
  Idents(obj) <- 'seurat_clusters'
  cluster.averages <- AverageExpression(obj, assays = assay, slot = slot, features = markergene, return.seurat = TRUE)
  if (assay == 'SCT') {
    scale.data <- cluster.averages@assays$SCT@data
  } else {
    scale.data <- cluster.averages@assays$RNA@data
  }
  print(scale.data)
  
  cell_score <- sapply(names(markerList), function(x) {
    tmp <- scale.data[rownames(scale.data) %in% markerList[[x]], ]
    if (is.matrix(tmp)) {
      if (nrow(tmp) >= 2) {
        res <- apply(tmp, 2, max)
        return(res)
      } else {
        return(rep(-2, ncol(tmp)))
      }
    } else {
      return(tmp)
    }
  })
  
  print(cell_score)
  
  celltypeMap <- apply(cell_score, 1, function(x) {
    colnames(cell_score)[which(x == max(x))]
  }, simplify = TRUE)
  
  obj@meta.data$cellType_auto <- plyr::mapvalues(x = as.character(obj@active.ident), from = names(celltypeMap), to = celltypeMap)
  
  return(obj)
}

B.cells <- c("MS4A1","CD19","CD1C","CD79A")
T.cells <- c("CD3D", "CD3E", "CD8A")
CD4T <-c("CD4","CD3")
Endothelial.cells <- c("CD31","VWF") 
DC.cell<-c("CD1C","CD11C")
Fibroblasts <- c("COL1A1","DCN")
Mac<-c("CD68","CD163")
NK.cells <- c("CD56", "GZMB") 
markergenelist <- list(T_cells=T.cells,DCs=DC.cell,B_cells=B.cells,Fibroblasts=Fibroblasts,
                       NK_cells=NK.cells,CD4T=CD4T,Mac=Mac,Endothelial=Endothelial.cells )
SKCM_obj <- cellAnnotation(obj=da,assay = "RNA",markerList=markergenelist)

DimPlot(SKCM_obj, reduction = "umap",
        group.by = "cellType_auto",label = T)
saveRDS()
##################         注释           ###############################
# B.cells <- c("MS4A1","CD19","CD1C","CD79A")
# T.cells <- c("CD3D", "CD3E", "CD8A")
# CD4T <-c("CD4","CD3")
# Endothelial.cells <- c("CD31","VWF") 
# Fibroblasts <- c("COL1A1","DCN")
# Mac<-c("CD68","CD163")
# NK.cells <- c("CD56", "GZMB") 
# Mon <- c("CD14","S100A9") 


genes_to_check = c('MS4A1','CD19','CD1C','CD79A','CD3D', 'CD3E', 'CD8A','CD4','CD3',
                   'CD31','VWF','COL1A1','DCN','CD68','CD163','CD56', 'GZMB','CD14','S100A9')
p_all_markers <- DotPlot(da, features = genes_to_check,
                         assay='RNA'  )  + coord_flip()

p_all_markers
ggsave(plot=p_all_markers,
       filename="check_all_marker_by_seurat_cluster.pdf",width = 12)
#######################   更换level         ###################################
library(dplyr)
library(stringr)
da=readRDS("./da_celltype.rds")
seurat_clusters <- da@meta.data$seurat_clusters
res3 <- da@meta.data$RNA_snn_res.0.3
new_df <- data.frame(seurat_clusters = as.character(seurat_clusters), res3 = res3)
#new_df$levels <- paste(new_df$T_celltype, new_df$seurat_clusters, sep = "_")
#da <- RenameIdents(da, "NK_cells"="T_cells")
new_df=distinct(new_df)
da=RenameIdents(da,setNames(new_df$res3,new_df$seurat_clusters))
levels(da)
table(da@meta.data$celltype)
##################################################
celltype=data.frame(ClusterID=0:10,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0,1,6),2]='T_cells'  
#celltype[celltype$ClusterID %in% 6,2]='NK_cells' 
celltype[celltype$ClusterID %in% 2,2]='B_cells'
celltype[celltype$ClusterID %in% c( 7,10 ),2]='Endothelial'
celltype[celltype$ClusterID %in% c(5,9),2]='Macrophage'  
celltype[celltype$ClusterID %in% c(3,4,8),2]='Fibroblasts'  
head(celltype)
table(celltype$celltype)
da@meta.data$celltype_old = "NA"
for(i in 1:nrow(celltype)){
 da@meta.data[which(da@meta.data$RNA_snn_res.0.3 == celltype$ClusterID[i]),'celltype_old'] <- celltype$celltype[i]}
table(da@meta.data$celltype_old)
head(da)
da@meta.data$celltype=da@meta.data$celltype_old
DimPlot(da, reduction = "umap",
        group.by = "celltype",label = T)
setwd("E:\\iraes-pathway工作总结\\图\\单细胞\\SKCM\\GSE72056")
saveRDS(da,"da_cluster.rds")
################################  差异   ###################################################
seurat_clusters <- da@meta.data$RNA_snn_res.0.3
celltype <- da@meta.data$celltype_old
new_df <- data.frame(seurat_clusters = as.character(seurat_clusters), celltype= celltype)
#new_df$levels <- paste(new_df$T_celltype, new_df$seurat_clusters, sep = "_")
new_df=distinct(new_df)
da=RenameIdents(da,setNames(new_df$seurat_clusters,new_df$celltype))
levels(da)
DefaultAssay(da) <- "RNA"
da.marker <- FindAllMarkers(da,only.pos = TRUE,min.pct = 0.5,logfc.threshold = 0.5)
head(da.marker)
write.csv(da.marker,"da.marker.csv")

gene=intersect(da.marker$gene,geneset$symbol)


################################    ssGSEA             ###########################################################
setwd("E:\\iraes-pathway工作总结\\图\\单细胞\\SKCM\\GSE72056")
#BiocManager::install("GSVA")
library(ggplot2)
library(tinyarray)
library(GSVA)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggpubr)

#输入分组信息，设置好样本对应的状态
da=readRDS("./da_celltype.rds")
mygroup <- data.frame(group = colnames(da), celltype= da$celltype)
row.names(mygroup) <- mygroup[,1]
mygroup <- mygroup[,-1]   #使mygroup最终为character 图4

#输入表达矩阵 不要取log（标化后十几内可以），不可以有负值和缺失值---所有基因全表达矩阵
exp=da@assays[["RNA"]]@counts
exp<- as.matrix(exp)    #注意将表达谱的data.frame转化为matrix，否则后续的gsva分析会报错

#####免疫细胞丰度计算
#导入免疫细胞基因集
geneset =rio::import("E:\\iraes-pathway工作总结\\图\\生存\\总的\\SKCM_sur_all.csv") #存储在文件夹中，要加载，加载后转化为list
##注意行名:转换时将行名转换为list名
geneset=geneset[,1]
genelist <- as.list(as.data.frame(gene)) ##需要转置后才能将行名转换为为list的名
# 使用na.omit()函数删除NA值
genelist <- lapply(genelist, na.omit)
lapply(genelist[1:3], head)   #转为list----成图3

#进行gsva分析
re <- gsva(exp, genelist, method="ssgsea",
           mx.diff=FALSE, verbose=FALSE)        #注意表达谱exp载入后需转化为matrix，前面已转换



#免疫细胞浸润分组绘图
draw_boxplot(re,mygroup)
#单基因相关性分析---目的基因与免疫细胞相关性
mygene <- c("PENK","GCGR",'RBP7','INHA')  #定义你的目的基因
nc = t(rbind(re,exp[mygene,]))  ;#将你的目的基因匹配到表达矩阵---行名匹配--注意大小写
m = rcorr(nc)$r[1:nrow(re),(ncol(nc)-length(mygene)+1):ncol(nc)]

##计算p值
p = rcorr(nc)$P[1:nrow(re),(ncol(nc)-length(mygene)+1):ncol(nc)]
head(p)


tmp <- matrix(case_when(as.vector(p) < 0.01 ~ "**",
                        as.vector(p) < 0.05 ~ "*",
                        TRUE ~ ""), nrow = nrow(p))

##绘制热图
p1 <- pheatmap(t(m),
               display_numbers =t(tmp),
               angle_col =45,
               color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
               border_color = "white",
               cellwidth = 20, 
               cellheight = 20,
               width = 7, 
               height=9.1,
               treeheight_col = 0,
               treeheight_row = 0)
#################################  AddModuleScore和AUcell      ###################################################
library(Seurat) 
library(tidyverse)
library(msigdbr) #提供基因集
install.packages("E:\\软件\\AUCell_1.26.0.zip")
BiocManager::install("AUCell")
library(AUCell)
#读取数据
#load("sce.anno.RData")
#载入颜色
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
colour1=c("#F08080","#87CEEB","#40E0D0","#DC143C","#FF6347")
#设置Assay
DefaultAssay(da) <- "RNA"
# 手动输入基因向量，并转为list形式 
#随意选择其中一条通路，转为list
WNT_features <- genelist
sce2 <- AddModuleScore(da,
                       features = WNT_features,
                       ctrl = 100,
                       name = "WNT_features")
head(sce2@meta.data)
#这里就得到了基因集评分结果，但是注意列名为 WNT_features1
colnames(sce2@meta.data)[10] <- 'WNT_Score'
####     AUCelll
cells_rankings <- AUCell_buildRankings(sce2@assays$RNA@data,splitByBlocks=TRUE) 

cells_AUC <- AUCell_calcAUC(genelist, cells_rankings, 
                            aucMaxRank=nrow(cells_rankings)*0.1)
geneSet <- "geneset"
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
#添加至metadata中
sce2$AUCell <- AUCell_auc
head(sce2@meta.data)
########   可视化
VlnPlot(sce2,features = 'WNT_Score', 
        pt.size = 0,group.by = "celltype")

ggboxplot(sce2@meta.data, x="celltype", y="AUCell", width = 0.6, 
          color = "black",#轮廓颜色
          fill="celltype",#填充
          palette = colour1,
          #palette =c("#E7B800", "#00AFBB"),#分组着色
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=1, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right") #图例
library(ggrepel)
#提取umap坐标数据
umap <- data.frame(sce2@meta.data, sce2@reductions$umap@cell.embeddings)
head(umap)

#1）计算每个celltype的median坐标位置
cell_type_med <- umap %>%
  group_by(celltype) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
#2）使用ggrepel包geom_label_repel 添加注释

ggplot(umap, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(colour = AUCell)) +
  ggrepel::geom_label_repel(
    aes(label = celltype),
    data = cell_type_med,
    fontface = "bold",
    point.padding = unit(0.5, "lines")
  ) +
  scale_colour_gradient(low = "#3F84B4", high = "#E27067")
saveRDS(sce2,"sce2_AUCell.rds")

###########################       拟时序       #############################################
library(monocle)
library(tidyverse)
library(Matrix)
library(monocle)
library(dplyr)
library(RColorBrewer)
display.brewer.all() #显示所有调色板

pbmc_assays <- pbmc@assays
HSMM_sample_sheet <- pbmc@meta.data #保存细胞信息
HSMM_gene_annotation <- pbmc_assays[["RNA"]]@meta.features #保存基因信息
HSMM_expre_matrix<- pbmc_assays[["RNA"]]@counts #保存矩阵信息

#========== 输入文件，生成cellDataSet对象，将Monocle把单细胞表达数据存放在CellDataSet对象里
pd <- new('AnnotatedDataFrame',data=HSMM_sample_sheet) # 细胞信息
fd <- new('AnnotatedDataFrame',data=HSMM_gene_annotation) # 基因信息
test=newCellDataSet(as(as.matrix(expr_matrix),'sparseMatrix'),phenoData = pd,featureData = fd)
#大数据集使用稀疏矩阵，节省内存，加快运算
test <- estimateSizeFactors(test) 
test <- estimateDispersions(test)
test=detectGenes(test,min_expr = 0.1) #计算每个基因在多少细胞中表达

marker_gene=import("E:\\iraes-pathway工作总结\\图\\生存\\总的\\SKCM_sur_all.csv")
test_ordering_genes=unique(marker_gene$symbol)
test=setOrderingFilter(HSMM,ordering_genes = test_ordering_genes) 
#指明哪些基因用于后续的聚类/排序
test=reduceDimension(test,reduction_method = "DDRTree",max_components = 2, norm_method = 'log',residualModelFormulaStr = "~num_genes_expressed") 
#residualModelFormulaStr减少其他因素的影响，比如不同样本、不同批次
test=orderCells(test)

plot_cell_trajectory(test,color_by = "celltype")
ggsave("celltype.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
plot_cell_trajectory(test,color_by = "State")
ggsave("State.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
plot_cell_trajectory(test,color_by = "Pseudotime")
ggsave("Pseudotime.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
plot_cell_trajectory(test,color_by = "celltype")+facet_wrap(~celltype,nrow=1)
ggsave("celltypeb.pdf",device = "pdf",width = 21,height = 9,units = c("cm"))

###########################
expressed_genes=row.names(subset(fData(test),num_cells_expressed>=10)) #在部分基因里面找
pseudotime_de <- differentialGeneTest(test[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]
states_de <- differentialGeneTest(test[expressed_genes,],
                                  fullModelFormulaStr = "~State")
states_de <- states_de[order(states_de$qval), ]

saveRDS(test, file = "test_monocle.rds")
write.table(pseudotime_de, file = "pseudotime_de.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
write.table(states_de, file = "states_de.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

BEAM_res=BEAM(test,branch_point = 1,cores = 1,progenitor_method="duplicate")
#会返回每个基因的显著性，显著的基因就是那些随不同branch变化的基因
#这一步很慢
BEAM_res=BEAM_res[,c("id","pval","qval")]
saveRDS(BEAM_res, file = "BEAM_res.rds")
BEAM_res=readRDS("./BEAM_res.rds")
BEAM_gene=BEAM_res[BEAM_res$id%in%intersect(BEAM_res$id,marker_gene$symbol),]
top100=BEAM_res[BEAM_res$id%in%top50,]
top100=top100[1:100,]
top50=c("FGD3","FCRL3","PARP15","CYTIP","CD3G","SP140","IL2RG","LCP1","CD48","PTPRCAP","GZMK","CD96","LCK","CD3D","SLA","HCST","EMB",
        "TYR","GPNMB","GZMB","CXCL13","CTLA4","SAMD3","TARP","CXCR3","GZMA","CD8A","CST7","TOX","CD27","LSP1","TRAT1","GZMM","IFNG")
BEAM_gene=rbind(BEAM_gene,top100)
tmp1=plot_genes_branched_heatmap(test[subset(BEAM_gene,qval < 1e-4)$id,],
                                 branch_point = 1,
                                 num_clusters = 4, #这些基因被分成几个group
                                 cores = 1,
                                 #branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 #hmcols = NULL, #默认值
                                #hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 #branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                 use_gene_short_name = F,
                                 show_rownames = T,
                                 return_heatmap = T #是否返回一些重要信息
)

pdf("branched_heatmap.pdf",width = 5,height = 6)
tmp1$ph_res
dev.off()

##################    细胞通讯                ##############################

library(CellChat)
library(tidyverse)
library(Seurat)
library(rio)
library(ggplot2)
options(stringsAsFactors = FALSE)


###### 保存图片的function
ps=function(filename,plot=FALSE,w=12,h=6){
  if (is.object(plot)) {
    print(plot)
  }
  plot <- recordPlot()
  pdf(file = filename, onefile = T, width = w, height = h)
  replayPlot(plot)
  dev.off()
}

###### 输入注释好的单细胞数据
setwd("E:\\iraes-pathway工作总结\\图\\单细胞\\SKCM\\GSE72056")
data=readRDS("./da_celltype.rds")

############################  所有时间一起  ####################################
setwd("./细胞通讯/")
if(T){

  ###### 创建cellchat对象
  cellchat=createCellChat(object = data,group.by = "celltype")
  summary(cellchat)
  str(cellchat)
  levels(cellchat@idents)
  ## 查看每个细胞类型有多少细胞
  groupsize=table(cellchat@idents)
  
  ###### 导入配体受体数据库
  #cellchatdb=CellChatDB.mouse # --------------------------------------------------> 小鼠的配体受体库
   cellchatdb=CellChatDB.human # -----------------------------------------------> 人类的配体受体库
  str(cellchatdb)
  # 包含"interaction","complex","cofactor","geneInfo"  
  colnames(cellchatdb$interaction)
  cellchatdb$interaction[1:4,1:4]
  head(cellchatdb$complex)
  head(cellchatdb$cofactor)
  head(cellchatdb$geneInfo)
  showDatabaseCategory(cellchatdb)
  
  ###### 从侧面来刻画细胞通讯相互作用(提取细胞和细胞之间的关系)
  unique(cellchatdb$interaction$annotation)
  cellchatdb.use=subsetDB(cellchatdb,search="Secreted Signaling")
  cellchat@DB=cellchatdb.use
  
  ###### 预处理
  cellchat=subsetData(cellchat)
  ## 识别每个细胞群中高表达的配体受体(findmarkers)
  cellchat=identifyOverExpressedGenes(cellchat) 
  ## 结果在cellchat@LR$LRsig
  cellchat=identifyOverExpressedInteractions(cellchat) 
  ## 找到配受体后映射到PPI上对@data.signaling表达之矫正，结果在@data.project
  #cellchat=projectData(cellchat,PPI.mouse)  #------------------------------------> 小鼠的PPI
   cellchat=projectData(cellchat,PPI.human)  #----------------------------------> 人类的PPI
  
  
  ###### 推断细胞通讯网络
  #### 推断配体受体水平细胞通讯网络，结果在@net下面
  ## 如果不想用上一步PPI网络矫正的结果，可以设置raw.use = T
  cellchat=computeCommunProb(cellchat,raw.use = F,population.size = T)
  cellchat=filterCommunication(cellchat,min.cells = 10)
  df.net=subsetCommunication(cellchat)
  ## 计算配受体对相互作用的通信概率推断信号通路水平上的通信概率
  export(df.net,"net_lr.csv")
  
  #### 推断信号通路水平的细胞通讯网络，结果在@netP下面
  ## 我们可以计算链路的数量或者汇总通信概率来计算细胞间的聚合通信网络
  cellchat=computeCommunProbPathway(cellchat)
  df.netp=subsetCommunication(cellchat,slot.name = "netP")
  export(df.netp,"df.netp.csv")
  
  
  ###### 细胞互作关系展示
  #### 所有细胞群总体展示：细胞互作数量和强度统计分析
  ## 统计细胞和细胞直接按通信的数量（有多少个配体受体对）和强度（概率）
  cellchat=aggregateNet(cellchat)
  ## 计算每种细胞各有多少个
  groupsize=as.numeric(table(cellchat@idents))
  ## 可视化
  par(mfrow=c(1,2),xpd=T)
  netVisual_circle(cellchat@net$count,
                   vertex.weight = groupsize,
                   weight.scale = T,
                   label.edge = F,
                   title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight,
                   vertex.weight = groupsize,
                   weight.scale = T,
                   label.edge = F,
                   title.name = "Interaction weights/strength")
  ps("所有细胞群互作数量和强度统计可视化.pdf")
  
  ## 检查每种细胞发出的信号
  ## 如果出不来图使用dev.new()
  cellchat=updateCellChat(cellchat)
  # require(devtools)
  ## 自己原本的igraph版本
  # install_version("igraph", version = "1.5.1", repos = "http://cran.us.r-project.org")
  # install_version("igraph", version = "1.4.0", repos = "http://cran.us.r-project.org")
  ## 展示通讯个数
  mat=cellchat@net$count
  par(mfrow=c(3,3),xpd=T)
  for(i in 1:nrow(mat)){  #------------------------------------------------------> igraph1.3.5
    # i=2
    print(i)
    mat2=matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
    mat2[i,]=mat[i,]
    netVisual_circle(mat2,
                     vertex.weight = groupsize,
                     weight.scale = T,
                     # arrow.width = 0.2,
                     # arrow.size = 0.1,
                     edge.weight.max = max(mat),
                     title.name = rownames(mat)[i])
  }
  ps("检查每种细胞发出的信号数目的可视化.pdf",w=12,h = 15)
  ## 展示通讯强度
  mat=cellchat@net$weight
  par(mfrow=c(3,3),xpd=T)
  for(i in 1:nrow(mat)){  #------------------------------------------------------> igraph1.3.5
    # i=2
    print(i)
    mat2=matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
    mat2[i,]=mat[i,]
    netVisual_circle(mat2,
                     vertex.weight = groupsize,
                     weight.scale = T,
                     # arrow.width = 0.2,
                     # arrow.size = 0.1,
                     edge.weight.max = max(mat),
                     title.name = rownames(mat)[i])
  }
  ps("检查每种细胞发出的信号强度的可视化.pdf",w=12,h = 15)
  
  ## 这是分开出图
  if(F){
    # ####  单个信号通路或配体受体介导的细胞互作可视化
    # ## 查看有哪些信号通路
    # cellchat@netP$pathways
    # # [1] "GALECTIN" "IL2"      "MIF"      "CCL"      "PARs"     "CD40"    
    # # [7] "FASLG"    "CSF" 
    # ## 选择一个信号通路例如IL2
    # pathways.show=c("IL2")
    # levels(cellchat@idents)
    # # [1] "CD4+CD8+"         "CD4+Exhausted"    "CD4+Naive T(Tfh)"
    # # [4] "CD4+Treg"         "CD8+Cytotoxic"    "CD8+Naive"       
    # # [7] "CD8+Tcm"          "Trm" 
    # ## 可视化层次图
    # vertex.receiver=c(2:5)
    # par(mfrow=c(1,1),xpd=T)
    # netVisual_aggregate(cellchat,
    #                     signaling = pathways.show,
    #                     vertex.receiver = vertex.receiver,
    #                     layout = "hierarchy")
    # ps(paste0("层次图可视化",pathways.show,"通路.pdf"),w = 8,h = 6)
    # ## 可视化网络图
    # vertex.receiver=c(1:7)
    # par(mfrow=c(1,1),xpd=T)
    # netVisual_aggregate(cellchat,
    #                     signaling = pathways.show,
    #                     vertex.receiver = vertex.receiver,
    #                     layout = "circle")
    # ps(paste0("网络图可视化",pathways.show,"通路.pdf"),w = 6,h = 6)
    # ## 可视化和弦图
    # vertex.receiver=c(1:7)
    # par(mfrow=c(1,1),xpd=T)
    # netVisual_aggregate(cellchat,
    #                     signaling = pathways.show,
    #                     layout = "chord")
    # ps(paste0("和弦图可视化",pathways.show,"通路.pdf"),w = 6,h = 6)
    # ## 可视化热图
    # vertex.receiver=c(1:7)
    # par(mfrow=c(1,1),xpd=T)
    # netVisual_heatmap(cellchat,
    #                   signaling = pathways.show,
    #                   color.heatmap = "Reds")
    # ps(paste0("热图可视化",pathways.show,"通路.pdf"),w = 8,h = 8)
    
    # #### 配体-受体层级的可视化（计算各个配体-受体对对信号通路的贡献）
    # ## 计算配体受体对选定信号通路的贡献值（在这里就是查看哪条通路对（pathways.show）贡献最大）
    # pathways.show="GALECTIN"
    # netAnalysis_contribution(cellchat,signaling = pathways.show) -> gg
    # ggsave(paste0("配受体对对",pathways.show,"贡献值.pdf"),height = 5,width = 4)
    # ## 提取对 pathway.show 有贡献的所有配受体对
    # pairLR.ccl=extractEnrichedLR(cellchat,
    #                              signaling = pathways.show,
    #                              geneLR.return = FALSE)
    # ## 提取通路中贡献最大的配体受体对来展示（也可以选择其他的配体受体对）
    # LR.show=pairLR.ccl[1,]
    # vertex.receiver=c(2:5)
    # ## 层级图展示
    # par(mfrow=c(1,1),xpd=T)
    # netVisual_individual(cellchat,
    #                      signaling = pathways.show,
    #                      show.legend = pathways.show,
    #                      pairLR.use = LR.show,
    #                      vertex.receiver = vertex.receiver,
    #                      layout = "hierarchy")
    # ps(paste0("层级图可视化",pathways.show,"通路贡献最大配受体对",LR.show,".pdf"),w = 8,h = 6)
    # ## 网络图展示
    # par(mfrow=c(1,1),xpd=T)
    # netVisual_individual(cellchat,
    #                      signaling = pathways.show,
    #                      show.legend = pathways.show,
    #                      pairLR.use = LR.show,
    #                      layout = "circle")
    # ps(paste0("网络图可视化",pathways.show,"通路贡献最大配受体对",LR.show,".pdf"),w = 6,h = 6)
    # ## 和弦图
    # par(mfrow=c(1,1),xpd=T)
    # netVisual_individual(cellchat,
    #                      signaling = pathways.show,
    #                      pairLR.use = LR.show,
    #                      layout = "chord")
    # ps(paste0("和弦图可视化",pathways.show,"通路贡献最大配受体对",LR.show,".pdf"),w = 6,h = 6)
  }
  
  #### 批量保存每个信号通路的互作结果
  pathways.show.all=cellchat@netP$pathways
  dir.create("all_pathways_com_circle")
  path=getwd()
  setwd("all_pathways_com_circle")
  for(i in 1:length(pathways.show.all)){
    ## 绘制网络图
    netVisual(cellchat,
              signaling = pathways.show.all[i],
              out.format = c("pdf"),
              vertex.receiver = vertex.receiver,
              layout = "circle")   #-------------------------------------------->报错
    gg=netAnalysis_contribution(cellchat,signaling = pathways.show.all[i])
    ggsave(filename = paste0(".//",
                             pathways.show.all[i],
                             "_L-R_contribution.pdf"),
           plot = gg,width = 5,height = 2.5,)
    
    pathways.show=pathways.show.all[i]
    #######  通路在细胞互作中的作用
    if(T){
      vertex.receiver=c(2:4)
      par(mfrow=c(1,1),xpd=T)
      netVisual_aggregate(cellchat,
                          signaling = pathways.show,
                          vertex.receiver = vertex.receiver,
                          layout = "hierarchy")
      ps(paste0(pathways.show,"层次图可视化","通路.pdf"),w = 8,h = 6)
      ## 可视化网络图
      vertex.receiver=c(1:5)
      par(mfrow=c(1,1),xpd=T)
      netVisual_aggregate(cellchat,
                          signaling = pathways.show,
                          vertex.receiver = vertex.receiver,
                          layout = "circle")
      ps(paste0(pathways.show,"网络图可视化","通路.pdf"),w = 6,h = 6)
      ## 可视化和弦图
      vertex.receiver=c(1:5)
      par(mfrow=c(1,1),xpd=T)
      netVisual_aggregate(cellchat,
                          signaling = pathways.show,
                          layout = "chord")
      ps(paste0(pathways.show,"和弦图可视化","通路.pdf"),w = 6,h = 6)
      ## 可视化热图
      vertex.receiver=c(1:4)
      par(mfrow=c(1,1),xpd=T)
      netVisual_heatmap(cellchat,
                        signaling = pathways.show,
                        color.heatmap = "Reds")
      ps(paste0(pathways.show,"热图可视化","通路.pdf"),w = 8,h = 8)
    }
    
    ######  通路中最大贡献的配受体对
    if(T){
      ## 提取对 pathway.show 有贡献的所有配受体对
      pairLR.ccl=extractEnrichedLR(cellchat,
                                   signaling = pathways.show,
                                   geneLR.return = FALSE)
      ## 提取通路中贡献最大的配体受体对来展示（也可以选择其他的配体受体对）
      LR.show=pairLR.ccl[1,]
      
      ## 层级图展示
      vertex.receiver=c(2:4)
      par(mfrow=c(1,1),xpd=T)
      netVisual_individual(cellchat,
                           signaling = pathways.show,
                           show.legend = pathways.show,
                           pairLR.use = LR.show,
                           vertex.receiver = vertex.receiver,
                           layout = "hierarchy")
      ps(paste0(pathways.show,"层级图可视化","通路贡献最大配受体对",LR.show,".pdf"),w = 8,h = 6)
      ## 网络图展示
      vertex.receiver=c(1:5)
      par(mfrow=c(1,1),xpd=T)
      netVisual_individual(cellchat,
                           signaling = pathways.show,
                           show.legend = pathways.show,
                           pairLR.use = LR.show,
                           layout = "circle")
      ps(paste0(pathways.show,"网络图可视化","通路贡献最大配受体对",LR.show,".pdf"),w = 6,h = 6)
      ## 和弦图
      vertex.receiver=c(1:5)
      par(mfrow=c(1,1),xpd=T)
      netVisual_individual(cellchat,
                           signaling = pathways.show,
                           pairLR.use = LR.show,
                           layout = "chord")
      ps(paste0(pathways.show,"和弦图可视化","通路贡献最大配受体对",LR.show,".pdf"),w = 6,h = 6)
    }
    
    print(i)
  }
  setwd(path)
  
  #### 多个配体受体介导的细胞互作关系可视化
  ## 气泡图（全部配受体）
  gg=netVisual_bubble(cellchat,
                      sources.use = c(1:5),
                      targets.use = c(1:5),
                      remove.isolate = FALSE)
  ggsave("所有细胞配受体对.pdf",width = 12,height = 5,plot = gg)
  ## 气泡图（指定信号通路或者配体受体）
  ## signaling参数后面跟的是cellchat@netP$pathways的结果
  netVisual_bubble(cellchat,
                   sources.use = c(3,4,5,8),
                   targets.use = c(5,6,7),
                   signaling = c("CCL","IL2","CD40"),
                   remove.isolate = F)
  ## 气泡图（指定信号通路或者配体受体并指定细胞）
  pairLR.use=extractEnrichedLR(cellchat,signaling = c("CCL","IL2","CD40"))
  netVisual_bubble(cellchat,
                   sources.use = c(1:8),
                   targets.use = c(1:8),
                   pairLR.use = pairLR.use,
                   remove.isolate = T)
  ## 参与某条信号通路（cellchat@netP$pathways）的所有基因在细胞群众的表达情况
  plotGeneExpression(cellchat,signaling = "GALECTIN")
  plotGeneExpression(cellchat,signaling = "CD40",type = "dot")
}


################################  两组之间比较  #################################
setwd("E:\\cellchat-test\\sham_Clp48h")
if(T){
  ###### 所有时间预处理流程和两组比较流程不一样，为啥？？
  if(F){
    cellchat=createCellChat(object = data,group.by = "T_celltype")
    cellchat@DB=CellChatDB.mouse
    cellchat=subsetData(cellchat)
    cellchat=subsetData(cellchat)
    cellchat=identifyOverExpressedGenes(cellchat) 
    cellchat=identifyOverExpressedInteractions(cellchat) 
    cellchat=projectData(cellchat,PPI.mouse)  #------------------------------------> 小鼠的PPI
    cellchat=computeCommunProb(cellchat,raw.use = F,population.size = T) #raw.use为T不需要矫正
    cellchat=filterCommunication(cellchat,min.cells = 10)
    df.net=subsetCommunication(cellchat)
    export(df.net,"net_lr.csv")
    cellchat=computeCommunProbPathway(cellchat)
    df.netp=subsetCommunication(cellchat,slot.name = "netP")
    export(df.netp,"df.netp.csv")
  }
  
  sham=subset(data,sample=="Sham")
  clp48h=subset(data,sample=="Clp48h")
  
  cellchat=createCellChat(object = sham,group.by = "T_celltype")
  cellchat@DB=CellChatDB.mouse
  cellchat=subsetData(cellchat)
  cellchat=identifyOverExpressedGenes(cellchat)
  cellchat=identifyOverExpressedInteractions(cellchat)
  cellchat=computeCommunProb(cellchat,raw.use = T,population.size = T)
  cellchat=aggregateNet(cellchat)
  cellchat=netAnalysis_computeCentrality(cellchat,slot.name = "netP")
  sham=cellchat
  export(sham,"sham.cellchat.rds")
  
  cellchat=createCellChat(object = clp48h,group.by = "T_celltype")
  cellchat@DB=CellChatDB.mouse
  cellchat=subsetData(cellchat)
  cellchat=identifyOverExpressedGenes(cellchat)
  cellchat=identifyOverExpressedInteractions(cellchat)
  cellchat=computeCommunProb(cellchat,raw.use = T,population.size = T)
  cellchat=aggregateNet(cellchat)
  cellchat=netAnalysis_computeCentrality(cellchat,slot.name = "netP")
  clp48h=cellchat
  export(clp48h,"clp48h.cellchat.rds")
  
  
  
  
  
}










































cellAnnotation <- function(obj, markerList, assay = 'SCT', slot = 'data') {
  markergene <- unique(do.call(c, markerList))
  
  Idents(obj) <- 'seurat_clusters'
  cluster.averages <- AverageExpression(obj, assays = assay, slot = slot, features = markergene, return.seurat = TRUE)
  if (assay == 'SCT') {
    scale.data <- cluster.averages@assays$SCT@data
  } else {
    scale.data <- cluster.averages@assays$RNA@data
  }
  print(scale.data)
  
  cell_score <- sapply(names(markerList), function(x) {
    tmp <- scale.data[rownames(scale.data) %in% markerList[[x]], ]
    if (is.matrix(tmp)) {
      if (nrow(tmp) >= 2) {
        res <- apply(tmp, 2, max)
        return(res)
      } else {
        return(rep(-2, ncol(tmp)))
      }
    } else {
      return(tmp)
    }
  })
  
  print(cell_score)
  
  celltypeMap <- apply(cell_score, 1, function(x) {
    colnames(cell_score)[which(x == max(x))]
  }, simplify = TRUE)
  
  obj@meta.data$cellType_auto <- plyr::mapvalues(x = as.character(obj@active.ident), from = names(celltypeMap), to = celltypeMap)
  
  return(obj)
}

B.cells <- c("MS4A1","CD19","CD1C","CD79A")
T.cells <- c("CD3D", "CD3E", "CD8A")
CD4T <-c("CD4","CD3")
Endothelial.cells <- c("CD31","VWF") 
DC.cell<-c("CD1C","CD11C")
Fibroblasts <- c("COL1A1","DCN")
Mac<-c("CD68","CD163")
NK.cells <- c("CD56", "GZMB") 
markergenelist <- list(T_cells=T.cells,DCs=DC.cell,B_cells=B.cells,Fibroblasts=Fibroblasts,
                       NK_cells=NK.cells,CD4T=CD4T,Mac=Mac,Endothelial=Endothelial.cells )
SKCM_obj <- cellAnnotation(obj=da,assay = "RNA",markerList=markergenelist)

DimPlot(SKCM_obj, reduction = "umap",
        group.by = "cellType_auto",label = T)

remove.packages("Matrix")

remotes::install_version("Matrix", version = "1.6-1.1")
devtools::install_github("cole-trapnell-lab/HSMMSingleCell")
BiocManager::install("SingleR")
BiocManager::install("HSMMSingleCell")
BiocManager::install("RBGL")
devtools::install_github('junjunlab/scRNAtoolVis')

install.packages("E:\\软件\\SingleR_2.6.0.zip")
install.packages("E:\\软件\\monocle_2.32.0.zip")
install.packages("E:\\软件\\biocViews_1.72.0.zip")
install.packages("E:\\软件\\RBGL_1.80.0.zip")
install_github("LTLA/SingleR")
library(monocle)
install.packages("devtools")
library(devtools)
library("SingleR")
BiocManager::install("beachmat")
install.packages("beachmat")




VlnPlot(data, 
        features = c("ANXA1","TNFRSF9"),
        pt.size = 0,
        ncol = 2) 

