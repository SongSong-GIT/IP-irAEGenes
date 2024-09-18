###-----------------------------------------------------------------------------
types3 = c('G1', 'G2', 'G3')
save(types3, file = "D:/irAEs/types3.RData")
term=c()
dir = paste0("D:/Keer Zhang/irAEs/result2/cluster/")
file_list = list.files(path = dir, pattern = "*.RData",recursive = TRUE,full.names = TRUE)
for(f in 1:length(file_list)){
  load(file = file_list[f])
  if(f==1){
    term1 <- as.vector(statistic3_cluster$gene[statistic3_cluster$class==1])
    term2 <- as.vector(statistic3_cluster$gene[statistic3_cluster$class==2])
    term3 <- as.vector(statistic3_cluster$gene[statistic3_cluster$class==3])
  }else{
    te <- as.vector(statistic3_cluster$gene[statistic3_cluster$class==1])
    term1 <- union(term1, te)
    te <- as.vector(statistic3_cluster$gene[statistic3_cluster$class==2])
    term2 <- union(term2, te)
    te <- as.vector(statistic3_cluster$gene[statistic3_cluster$class==3])
    term3 <- union(term3, te)
  }
}
term <- list(term1, term2, term3)
names(term) = types3



###1.加载包
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(edgeR)
library(limma)
dir = "D:/irAEs/irAEsTCGA/"  

termMat_G1_FC <- matrix(0, nrow = length(term[[1]]), ncol = 20)
colnames(termMat_G1_FC) <- types1
rownames(termMat_G1_FC) <- term[[1]]
termMat_G2_FC <- matrix(0, nrow = length(term[[2]]), ncol = 20)
colnames(termMat_G2_FC) <- types1
rownames(termMat_G2_FC) <- term[[2]]
termMat_G3_FC <- matrix(0, nrow = length(term[[3]]), ncol = 20)
colnames(termMat_G3_FC) <- types1
rownames(termMat_G3_FC) <- term[[3]]

dir = "D:/irAEs/irAEsTCGA/"
file_list = list.files(path = dir, pattern = "*.tsv",recursive = TRUE,full.names = TRUE)  
file_list = file_list[-c(1,11,15,16,25,26)]
for(f in 1:length(file_list)){
  exp1 <- read.table(file_list[f], check.names = F, row.names = 1)
  sample <- exp1[1,]
  exp1$ensembl_ID <- unlist(str_split(rownames(exp1), "[.]", simplify=T))[,1]
  rownames(exp1) <- exp1$ensembl_ID
  ensembl_ID <- rownames(exp1)
  ensembl_ID <- ensembl_ID[-1]
  symbol_ID <- bitr(ensembl_ID, 
                    fromType = "ENSEMBL", 
                    toType = "SYMBOL", 
                    OrgDb="org.Hs.eg.db")
  exp1 <- merge(exp1, symbol_ID, by.x = 'ensembl_ID', by.y = 'ENSEMBL')
  exp1 <- exp1[!duplicated(exp1$SYMBOL),]
  rownames(exp1) <- exp1$SYMBOL
  exp1 <- exp1[,-c(1,ncol(exp1))]
  colnames(exp1) <- sample
  for(i in 1:ncol(exp1)){
    exp1[,i] <- as.numeric(exp1[,i])
  }
  exp1_disease <- exp1[,grep("0[0-9]A", colnames(exp1))]
  exp1_normal <- exp1[,grep("1[0-9]A", colnames(exp1))]
  expall <- cbind(exp1_disease, exp1_normal)
  group <- c(rep("disease",ncol(exp1_disease)), rep("normal",ncol(exp1_normal)))
  eset = expall 
  targets <- cbind(colnames(expall),group)
  targets <- as.data.frame(targets)
  colnames(targets) = c("FileName","Target")
  level <- unique(targets$Target)
  design <- model.matrix(~0+factor(targets$Target, levels=level))
  colnames(design) <- level
  cont.wt <- makeContrasts(disease-normal, levels=design) 
  fit <- lmFit(eset, design)
  fit2 <- contrasts.fit(fit, cont.wt) 
  fit2 <- eBayes(fit2) 
  tT = topTable(fit2, adjust = "fdr", n = Inf)
  tT = subset(tT, select = c("adj.P.Val", "P.Value", "logFC"))
  colnames(tT) = c("FDR", "P.Value", "logFC")
  
  #GROUP1--------------------------------------------------------------------------------------
  for(i in 1:nrow(termMat_G1_FC)){
    gr <- match(rownames(termMat_G1_FC)[i], rownames(tT))
    if(length(gr)!=0){
      termMat_G1_FC[i,f] <- tT$logFC[gr]
    }else
      termMat_G1_FC[i,f] <- 0
  }
  
  
  #GROUP2--------------------------------------------------------------------------------------
  for(i in 1:nrow(termMat_G2_FC)){
    gr <- match(rownames(termMat_G2_FC)[i], rownames(tT))
    if(length(gr)!=0){
      termMat_G2_FC[i,f] <- tT$logFC[gr]
    }else
      termMat_G2_FC[i,f] <- 0
  }
  
  
  #GROUP3--------------------------------------------------------------------------------------
  for(i in 1:nrow(termMat_G3_FC)){
    gr <- match(rownames(termMat_G3_FC)[i], rownames(tT))
    if(length(gr)!=0){
      termMat_G3_FC[i,f] <- tT$logFC[gr]
    }else
      termMat_G3_FC[i,f] <- 0
  }
  
  
}



