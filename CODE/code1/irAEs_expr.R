load(file = "D:/irAEs/UCEC_esetData1.RData")

library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(edgeR)
library(limma)
dir = "D:/Keer Zhang/irAEs/irAEsTCGA/"
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
  eset = eset[rownames(esetData1[[f]]),]
  write.csv(eset, file = paste0("D:/Keer Zhang/irAEs/result2/esetData/", types1[f], "_esetData.csv"))
}
