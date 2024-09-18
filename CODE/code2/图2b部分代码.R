rm(list = ls())
load(file = "E:\\iraes-pathway工作总结\\irAEs4\\TCGA_hash1.RData")
load(file = "E:\\iraes-pathway工作总结\\irAEs4\\pathwaName.RData")


library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(hash)
library(ConsensusClusterPlus)
library(tinyarray)
library(ggplot2)

kru=c()
pairwil1=c()
pairwil2=c()
esetData1=c()

dir = "E:\\iraes-pathway工作总结\\irAEs1\\exp1_disease/"
file_list = list.files(path =dir, pattern = "*.RData",recursive = TRUE,full.names = TRUE)      
for(h in 1:length(TCGA_hash1)){
  h=1
  TCGA_hash = TCGA_hash1[[h]]
  statistic2 <- matrix(0, nrow = length(hash::keys(TCGA_hash)), ncol = 16)
  rownames(statistic2) <- as.vector(hash::keys(TCGA_hash))
  colnames(statistic2) <- pathwayName
  for(key in 1:length(hash::keys(TCGA_hash))){
    value <- table(as.vector(hash::values(TCGA_hash, hash::keys(TCGA_hash)[key])))
    for(each in 1:length(value)){
      statistic2[key, names(value)[each]] <- as.vector(value[each])
    }
  }
  
  
  ###Consensus Clustering
  name1 = "E:\\iraes-pathway工作总结\\irAEs\\图2"
  title1 = paste0("./", "Consensus_hc")
  title2 = paste0("./", "Consensus_pam")
  setwd("E:\\iraes-pathway工作总结\\irAEs\\图2")
  statistic3 <- apply(statistic2, 1, sum)
  c2 <-  ConsensusClusterPlus(d = dist(statistic3),
                              maxK = 8,
                              reps = 100,
                              pItem = 0.8, pFeature = 1,
                              clusterAlg = "pam",
                              seed = 102,
                              title = title2,
                              innerLinkage = "complete",
                              plot = "pdf")
  
}

