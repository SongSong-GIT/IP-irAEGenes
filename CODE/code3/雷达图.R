#加载包

library(AnnotationDbi)
library(org.Hs.eg.db)#基因注释包
library(clusterProfiler)#富集包
library(dplyr)
library(ggplot2)#画图包

setwd("E:\\iraes-pathway工作总结\\图")
########################寻找  hub_gene   #############################
library(progress)
library(RISmed)

dir = 'E:\\iraes-pathway工作总结\\imm\\result\\esetData'
file_list = list.files(path = dir, pattern = "*.csv",recursive = TRUE,full.names = TRUE)  
result_all <- data.frame()
for (f in 1:20) {
  exp_data <- read.csv(file_list[f], row.names = 1)
  
  if (f == 1) {
    k <- as.data.frame(rownames(exp_data))
    colnames(k) <- "gene"
    result_all <- rbind(result_all, k)
    result_all[,"count"] <- 1
  } else {
    k <- as.data.frame(rownames(exp_data))
    colnames(k) <- "gene"
    
    for (i in 1:nrow(k)) {
      gene_name <- k[i, "gene"]
      if (gene_name %in% result_all[, "gene"]) {
        # If the gene is already in the result_all, increment its count
        result_all[result_all$gene == gene_name, "count"] <- result_all[result_all$gene == gene_name, "count"] + 1
      } else {
        # If the gene is not in the result_all, add it with count = 1
        result_all <- rbind(result_all, data.frame(gene = gene_name, count = 1))
      }
    }
  }
}
hub_gene=result_all
#rm(k,i,gene_name,f,cyto_data,result_all,hub_data_bySymbol,exp_data)
library(progress)
library(RISmed)
hub_gene[,"pubmed"]=0
pb <- progress_bar$new(total = nrow(hub_gene))
system.time({
  for (i  in 1:nrow(hub_gene)) {
    pb$tick()
    search_topic <- paste("(",hub_gene[i,1],")")
    Sys.sleep(0.1)
    search_query <- EUtilsSummary(search_topic,db="pubmed")
    #retmax=10000,datetype='pdat', mindate=2000, maxdate=2023)
    hub_gene$pubmed[i]<-search_query@count
  }
})
hub=hub_gene[which(hub_gene$count>14),]
#hub=import("./hub_gene.txt")
write.table(result_all,"all_gene.txt",sep = "\t",row.names = F,quote = F)
write.table(hub,"hub_gene.txt",sep = "\t",row.names = F,quote = F)
######################   富集        ####################################
setwd("E:\\iraes-pathway工作总结\\irAEs\\图3")
#2、基因id转换，kegg和go富集用的ID类型是ENTREZID）
diff<-import("./all_gene.txt")
gene.df <- bitr(diff$gene,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)#TCGA数据框如果没有进行基因注释，那么fromType应该是Ensembl，各种ID之间可以互相转换,toType可以是一个字符串，也可以是一个向量，看自己需求                     
gene <- gene.df$ENTREZID

#1、KEGG富集
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism= "human",qvalueCutoff = 0.05,pvalueCutoff=0.05)

#2、可视化
###柱状图
hh <- as.data.frame(kk)#自己记得保存结果哈！
rownames(hh) <- 1:nrow(hh)
write.csv(hh,"kegg.result.csv")
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
hh=hh[1:20,]
ggplot(hh,aes(y=order,x=Count,fill=p.adjust))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  #coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "red",high ="blue" )+#颜色自己可以换
  labs(title = "KEGG Pathways Enrichment",
       x = "Gene numbers", 
       y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()
ggsave("./Down.KEGG.bar.pdf")
##########
library(dplyr)

df=hh[,c("ID","Description","p.adjust","Count","GeneRatio")]
df=df[order(df$Count,decreasing = T),]
df=df[1:20,]
df$Description=str_to_sentence(df$Description)
df$Description=factor(df$Description,df$Description)

## 气泡图
ggplot(df,aes(x=Count, y=Description)) + 
  geom_point(aes(size=Count,color=-1*log10(p.adjust)))+  #点的大小根据Count数变化，颜色根据P值变化
  scale_colour_gradient(low="blue", high="red")+  #设置图例
  labs(
    color=expression(-log[10](p.adjust)),
    size="Num of Genes",
    x="Num of Genes")+
  theme_bw()+  #设置背景
  theme(
    axis.text.y = element_text(size = rel(1)),
    axis.title.x = element_text(size=rel(1.0)),
    axis.title.y = element_blank())
ggsave("./Down.KEGG.point.pdf")

################################################################

hub_gene=import("./hub_gene.txt")
active=import("E:\\iraes-pathway工作总结\\irAEs\\图2\\Active.txt")
mild=import("E:\\iraes-pathway工作总结\\irAEs\\图2\\Mild.txt")
silent=import("E:\\iraes-pathway工作总结\\irAEs\\图2\\Silent.txt")
############           active         ######
cancer = unique(active$type)
cancer_active = list()

for (i in 1:length(cancer)) {
  cancer_type_data = active[active$type == cancer[i],]
  cancer_active[[i]] = cancer_type_data
}

names(cancer_active) = cancer
intersect(cancer_active[["BLCA"]][["gene"]],hub_gene$gene)
intersection_counts = list()

for (cancer_type in names(cancer_active)) {
  intersect_genes = intersect(cancer_active[[cancer_type]]$gene, hub_gene$gene)
  intersection_counts[[cancer_type]] = length(intersect_genes)
}

# Convert the list to a data frame for easier viewing if needed
intersection_counts_df = data.frame(
  CancerType = names(intersection_counts),
  IntersectionCount = unlist(intersection_counts)
)
a=intersection_counts_df
############           mild         ######
cancer = unique(mild$type)
cancer_mild = list()

for (i in 1:length(cancer)) {
  cancer_type_data = mild[mild$type == cancer[i],]
  cancer_mild[[i]] = cancer_type_data
}

names(cancer_mild) = cancer
intersect(cancer_mild[["BLCA"]][["gene"]],hub_gene$gene)

intersection_counts = list()

for (cancer_type in names(cancer_mild)) {
  intersect_genes = intersect(cancer_mild[[cancer_type]]$gene, hub_gene$gene)
  intersection_counts[[cancer_type]] = length(intersect_genes)
}

# Convert the list to a data frame for easier viewing if needed
intersection_counts_df = data.frame(
  CancerType = names(intersection_counts),
  IntersectionCount = unlist(intersection_counts)
)
b=intersection_counts_df
############           silent        ######
cancer = unique(silent$type)
cancer_silent = list()

for (i in 1:length(cancer)) {
  cancer_type_data = silent[silent$type == cancer[i],]
  cancer_silent[[i]] = cancer_type_data
}

names(cancer_silent) = cancer
intersect(cancer_silent[["BLCA"]][["gene"]],hub_gene$gene)
intersection_counts = list()

for (cancer_type in names(cancer_silent)) {
  intersect_genes = intersect(cancer_silent[[cancer_type]]$gene, hub_gene$gene)
  intersection_counts[[cancer_type]] = length(intersect_genes)
}

# Convert the list to a data frame for easier viewing if needed
intersection_counts_df = data.frame(
  CancerType = names(intersection_counts),
  IntersectionCount = unlist(intersection_counts)
)
c=intersection_counts_df

a=as.data.frame(t(a))
a=a[-1,]
a$group="active"

b=as.data.frame(t(b))
b=b[-1,]
b$group="mild"

c=as.data.frame(t(c))
c=c[-1,]
c$group="silent"

dt=rbind(a,b,c)

dt$BLCA=as.numeric(dt$BLCA)
dt$BRCA=as.numeric(dt$BRCA)
dt$CESC=as.numeric(dt$CESC)
dt$CHOL=as.numeric(dt$CHOL)
dt$COAD=as.numeric(dt$COAD)
dt$ESCA=as.numeric(dt$ESCA)
dt$GBM=as.numeric(dt$GBM)
dt$HNSC=as.numeric(dt$HNSC)
dt$KIRC=as.numeric(dt$KIRC)
dt$LIHC=as.numeric(dt$LIHC)
dt$LUAD=as.numeric(dt$LUAD)
dt$LUSC=as.numeric(dt$LUSC)
dt$PAAD=as.numeric(dt$PAAD)
dt$PRAD=as.numeric(dt$PRAD)
dt$READ=as.numeric(dt$READ)
dt$SARC=as.numeric(dt$SARC)
dt$SKCM=as.numeric(dt$SKCM)
dt$STAD=as.numeric(dt$STAD)
dt$THCA=as.numeric(dt$THCA)
dt$UCEC=as.numeric(dt$UCEC)




write.table(dt,"图3雷达图.txt",row.names = T,quote = F,sep = "\t")

rownames(dt)=1:3
'#88bd6d', '#354d7d', '#14927c'
library(ggradar)
mycol <- c("#EBD7A5","#91CBCF","#F8CCCD")

ggradar(dt,
        grid.min = 5,
        grid.mid = 45,
        grid.max = 85)

#自定义配色
dt=import("./图3雷达图.xlsx")
dt=dt[,-1]
suncun=import("E:\\iraes-pathway工作总结\\图\\3\\前部\\生存雷达.csv")
dt$group <- factor(dt$group, levels = dt$group) #转化为因子,固定绘图顺序

ggradar(dt,
        grid.min = 5,
        grid.mid = 45,
        grid.max = 85,
        values.radar = NA, #不显示网格线椭圆对应数值，可不显示
        gridline.mid.colour = 'grey', #中间网格线颜色
        axis.label.size = 5.5, #文本标签大小
        group.colours = mycol,
        group.point.size = 4,#点大小
        group.line.width = 1, #线宽
        background.circle.colour = '#D7D6D1', #背景填充色
        background.circle.transparency = 0.2, #背景填充不透明度
        legend.position = 'bottom', #图例位置
        legend.text.size = 14 #图例标签大小
)

# 
# 
# 
# color=c("#91CBCF","#EBD7A5","#F8CCCD")
# ggplot(a, aes(CancerType, weight = IntersectionCount, fill = color[3])) +
#   geom_bar(fill = color[3], color = "black", width = .7, position = 'stack', linewidth = .1) +
#   labs( y = NULL) +
#   scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 90, color = "black"),
#         axis.title.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.line.y = element_blank(),
#         axis.ticks = element_line(color = "black", linewidth = .1),
#         axis.line = element_line(color = "black", linewidth = .1),
#         legend.text = element_text(color = "black"),
#         legend.title = element_text(color = "black"),
#         plot.margin = unit(c(1,1,1,1),'cm'),
#         legend.position = 'none')
# 













