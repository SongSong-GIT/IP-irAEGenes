c()
setwd("E:\\iraes-pathway工作总结\\图\\生存\\图")

dir = "E:\\iraes-pathway工作总结\\图\\生存\\总的"
file_list = list.files(path = dir, pattern = "*",recursive = TRUE,full.names = TRUE)  
result_all <- data.frame()
f=1
for (f in 1:20) {
  exp_data <- read.csv(file_list[f], row.names = 1)
  
  if (f == 1) {
    k <- as.data.frame(rownames(exp_data))
    colnames(k) <- "gene"
    # file_path <- file_list[f]
    # file_name <- gsub(".*/([^/]+)_.*", "\\1", file_path)
    # k$type <- file_name
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
library(rio)
library(progress)
library(RISmed)
write.csv(result_all,"result_gene.csv",quote = F,row.names = F)


setwd("E:\\iraes-pathway工作总结\\图\\生存\\免疫")
# 加载包
library(dplyr)
library(tidyverse)
# 读取基因列表和8个细胞系表达谱数据，创建数据框
gene_list = import("./result_gene.csv")
cell_lines=list.files("./免疫细胞系数据")
cell_lines=sub(".txt", "", cell_lines)
exp_data <- lapply(cell_lines, function(x) read.table(paste0("./免疫细胞系数据/", x,".txt"), header = TRUE)) # 读取每个细胞系的表达谱数据

group=lapply(exp_data,ncol)
group=as.data.frame(group)
group=group-1
group=as.data.frame(t(group))
rownames(group)=1:19
group=rep(cell_lines,group$V1)
group=as.factor(group)
# group <- as.factor(group)
# write.table(levels(group),"anno.txt",row.names = F,quote = F)

exp_data <- lapply(exp_data, function(x) {
  rownames(x)=x$ID_REF
  x=x[,-1]
})

exp_data <- Reduce(cbind, exp_data) # 合并数据
exp <- subset(exp_data, rownames(exp_data) %in% gene_list$gene) # 筛选出基因列表中的基因
intersect(gene_list$gene,rownames(exp_data))

library(ggplot2) 
library(tidyr)
exp=as.data.frame(t(exp))
exp=cbind(exp,group)

# a=exp_data[,c(1,221)]
# ggplot(exp_data,aes(x=group,y=exp_data$`A1BG-AS1`))+geom_boxplot(outlier.shape = NA)

test=function(x){
  a=with(exp,kruskal.test(x , g=group, data = exp,p.adj="bon",group=T))
  return(a$p.value)
}
result=apply(exp,2,test)
result=as.data.frame(result)

boxplot(exp[,"SLFN12L"]~group,data = exp)

exp$group=as.character(exp$group)
gene=aggregate(.~group,data = exp,mean)

max_cell=function(x){
  max_c=max(x)
  index=which(x%in%max_c)
  return(gene$group[index])
}
cell_high=apply(gene,2,max_cell)
cell_high=cbind(colnames(gene),cell_high)
cell_high=cell_high[-1,]
select_cell=c("CD4 T cell activated","CD4 T cell resting","CD8 T cell activated",
              "CD8 T cell resting","NKT activated","T gamma delta","T helper 17")
cell_high=as.data.frame(cell_high)
T_cell_lnc=cell_high[which(cell_high$cell_high%in%select_cell),1]
cell_high[which(cell_high$V1%in%select_gene),]

plot.data=exp[,"SLFN12L"]
plot.data=cbind(plot.data,exp$group)
colnames(plot.data)=c("value","group")
plot.data=as.data.frame(plot.data)
plot.data$value=as.numeric(plot.data$value)
tmp=aggregate(value~group,plot.data,mean)
tmp=tmp[order(tmp$value,decreasing = T),]
plot.data$group=factor(plot.data$group,levels = tmp$group)
library(vioplot)
# vioplot(value~group, data = plot.data,
#         col=c("#d67d79","#9fb9cd","#c5d832","#98afd9","#4d79b6","#ba86b4",
#               "#c2dabd","#d57732","#b4c950","#97c7cb","#8cb0b0","#e2b719",
#               "#dc736b","#7a9ed5","#c5c37b","#94a8c5","#db7aa1","#dfcc7d",
#               "#e4b3b3"),las=2)
vioplot(value~group, data = plot.data,
        col=c("#E99D96","#9DC1E2","#cACEDA","#DAD7C4","#035697","#549BD0",
              "#CFE3F3","#d57732","#b4c950","#D6675B","#FBCCB4","#7EB8D7",
              "#F6B529","#FADC7E","#95CB62","#77B3DF","#FEECAF","#dfcc7d",
              "#DDC8FC"),las=2)