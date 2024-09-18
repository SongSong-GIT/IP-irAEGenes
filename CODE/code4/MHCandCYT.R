dir <- "E:\\iraes-pathway工作总结\\imm\\test\\imm\\"
file_list <- list.files(path = dir, pattern = "*.RData", recursive = TRUE, full.names = TRUE)

TIMER <- data.frame()
SORT <- data.frame()
XCELL <- data.frame()
EPIC <- data.frame()
for (i in 1:length(file_list)) {
  load(file = file_list[i]) 

  # 遍历所有的数据列表，提取名字和数量，并添加到数据框中
  for (j in names(imme_result_list)) {
    if (grepl("^Mac.*(TIMER)$", j)) {
      data_timer <- imme_result_list[[j]]
      gene_timer <- data_timer$gene
      gene_timer <- as.data.frame(gene_timer)
      colnames(gene_timer) <- "gene"
      TIMER <- rbind(TIMER, gene_timer)
    } else if (grepl("^Mac.*(CIBERSORT)$", j)) {
      data_sort <- imme_result_list[[j]]
      gene_sort <- data_sort$gene
      gene_sort <- as.data.frame(gene_sort)
      colnames(gene_sort) <- "gene"
      SORT <- rbind(SORT, gene_sort)
    }else if (grepl("^Mac.*(XCELL)$", j)) {
      data_xcell <- imme_result_list[[j]]
      gene_xcell <- data_xcell$gene
      gene_xcell <- as.data.frame(gene_xcell)
      colnames(gene_xcell) <- "gene"
      XCELL <- rbind(XCELL, gene_xcell)
    }else if (grepl("^Mac.*(EPIC)$", j)) {
      data_epic <- imme_result_list[[j]]
      gene_epic <- data_epic$gene
      gene_epic <- as.data.frame(gene_epic)
      colnames(gene_epic) <- "gene"
      EPIC <- rbind(EPIC, gene_epic)
    }
    
    
    
  }
  TIMER <- unique(TIMER$gene)
  TIMER <- as.data.frame(TIMER)
  colnames(TIMER) <- "gene"
  
  SORT <- unique(SORT$gene)
  SORT <- as.data.frame(SORT)
  colnames(SORT) <- "gene"
  
  XCELL <- unique(XCELL$gene)
  XCELL <- as.data.frame(XCELL)
  colnames(XCELL) <- "gene"
  
  EPIC <- unique(EPIC$gene)
  EPIC <- as.data.frame(EPIC)
  colnames(EPIC) <- "gene"
  
  
  # 这里可以对 TIMER 和 SORT 进行其他操作，例如保存到文件中或进行进一步分析
}
############################    MHC    ##########################################
library(readxl)
library(readr)
library(rio)
gc()
setwd("E:\\iraes-pathway工作总结\\图\\生存\\免疫")
load("E:/iraes-pathway工作总结/imm/result/exp_446.RData")
load("E:/iraes-pathway工作总结/imm/result/sample.RData")
anno=read_csv("33个癌症的MHC和CYT.csv")
anno=as.data.frame(anno)
rownames(anno)=substr(anno$Sample,1,15)
anno$Sample=substr(anno$Sample,1,15)
hub_lnc=read.table("E:\\iraes-pathway工作总结\\图\\3\\前部\\hub_gene.txt",header = T)
data=import("E:\\iraes-pathway工作总结\\图\\生存\\免疫\\DATA.xlsx")
hub_lnc=hub_lnc[which(hub_lnc$gene%in%data$hub_gene),]
colnames(hub_lnc)=c("lnc","count","pubmed")
#hub_lnc=hub_lnc[which(hub_lnc$pubmed>18),]
#hub_lnc=hub_lnc[which(hub_lnc$pubmed!=21),]

other_lnc=import("E:\\iraes-pathway工作总结\\图\\生存\\免疫\\other_lnc.txt")
intersect(hub_lnc$lnc,other_lnc$other_gene)
intersect(hub_lnc$SYMBOL,rownames(dat1))
# hub_lnc=other_lnc

dir <- "E:\\iraes-pathway工作总结\\irAEs1\\exp1_disease/"
file_list <- list.files(path = dir, pattern = "*RData", recursive = TRUE, full.names = TRUE)  
cancer <- c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KIRC","LIHC",
            "LUAD","LUSC","PAAD","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC")

dat_cor <- data.frame()
dat_p <- data.frame()

s_cor <- function(x) {
  test <- cor.test(as.numeric(x), as.numeric(dat2[1,]), method = "spearman")
  estimate <- test$estimate
  pvalue <- test$p.value
  r <- c(R = estimate, p = pvalue)
  return(r)
}

for (i in 1:length(cancer)) {
  #i=14
  cancer_type <- cancer[i]
  dat1 <- import(file_list[i])
  colnames(dat1) <- substr(colnames(dat1), 1, 15)
  coo_sample <- intersect(rownames(anno), colnames(dat1))
  
  # 提取数据
  

  dat1 <- dat1[, coo_sample]
  dat1 <- dat1[other_lnc$other_gene, ]
  na_rows <- grep("^NA", rownames(dat1))
  dat1[na_rows, ] <- 0
  
  rownames(dat1) <- other_lnc$other_gene
  
  # 处理注释数据
  dat2 <- anno[coo_sample, ]
  dat2 <- as.data.frame(t(dat2))
  dat2 <- dat2[-1,]
  
  # 计算相关性
  result <- as.data.frame(t(apply(dat1, 1, s_cor)))
  
  # 设置行名和列名
  tmp <- result$R
  tmp <- as.data.frame(tmp)
  tmp <- as.data.frame(t(tmp))
  rownames(tmp) <- cancer_type
  colnames(tmp) <- rownames(result)
  dat_cor <- rbind(dat_cor, tmp)
  
  tmp <- result$p
  tmp <- as.data.frame(tmp)
  tmp <- as.data.frame(t(tmp))
  rownames(tmp) <- cancer_type
  colnames(tmp) <- rownames(result)
  dat_p <- rbind(dat_p, tmp)
}

plot.data.cor <- list("MHC" = dat_cor)
plot.data.p <- list("MHC" = dat_p)


##############  CYT     #######################
dat_cor <- data.frame()
dat_p <- data.frame()
#cancer_type="PRAD"
s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[2,]), method = "spearman")
  estimate=test$estimate
  pvalue=test$p.value
  r<-c(R=estimate,p=pvalue)
  return(r)
}
for (j in 1:length(cancer)) {
  #i=14
  cancer_type <- cancer[j]
  dat1 <- import(file_list[j])
  colnames(dat1) <- substr(colnames(dat1), 1, 15)
  coo_sample <- intersect(rownames(anno), colnames(dat1))
  
  # 提取数据
  
  
  dat1 <- dat1[, coo_sample]
  dat1 <- dat1[other_lnc$other_gene, ]
  na_rows <- grep("^NA", rownames(dat1))
  dat1[na_rows, ] <- 0
  
  rownames(dat1) <- other_lnc$other_gene
  
  # 处理注释数据
  dat2 <- anno[coo_sample, ]
  dat2 <- as.data.frame(t(dat2))
  dat2 <- dat2[-1,]
  
  # 计算相关性
  result <- as.data.frame(t(apply(dat1, 1, s_cor)))
  
  # 设置行名和列名
  tmp <- result$R
  tmp <- as.data.frame(tmp)
  tmp <- as.data.frame(t(tmp))
  rownames(tmp) <- cancer_type
  colnames(tmp) <- rownames(result)
  dat_cor <- rbind(dat_cor, tmp)
  
  tmp <- result$p
  tmp <- as.data.frame(tmp)
  tmp <- as.data.frame(t(tmp))
  rownames(tmp) <- cancer_type
  colnames(tmp) <- rownames(result)
  dat_p <- rbind(dat_p, tmp)
}
tmp.cor=list("CYT"=dat_cor)
tmp.p=list("CYT"=dat_p)
plot.data.cor=append(plot.data.cor,tmp.cor)
plot.data.p=append(plot.data.p,tmp.p)

##########################################################
heatdata=as.matrix(plot.data.cor$MHC)
anno_melt=as.matrix(plot.data.p$MHC)

library(ggplot2)
library(reshape2)
heatdata_melt=melt(heatdata)
anno_melt=melt(anno_melt)
heatdata_melt$star <- ifelse(anno_melt$value< 0.05, "*", "")
ggplot(data = heatdata_melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colors = c("#5390b5","#FFFFFF","#d56e5e"),
                       limits=c(-1,1)) +  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 6, hjust = 1))+
  coord_fixed()+
  geom_text(aes(label = star), size = 4)

heatdata=as.matrix(plot.data.cor$CYT)
anno_melt=as.matrix(plot.data.p$CYT)

library(ggplot2)
library(reshape2)
heatdata_melt=melt(heatdata)
anno_melt=melt(anno_melt)
heatdata_melt$star <- ifelse(anno_melt$value< 0.05, "*", "")
ggplot(data = heatdata_melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colors = c("#5390b5","#FFFFFF","#d56e5e"),
                       limits=c(-1,1)) +  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 6, hjust = 1))+
  coord_fixed()+
  geom_text(aes(label = star), size = 4)
