library(splines)
library(survival)
library(foreign)
library(MASS)
library(nnet)
setwd("E:\\iraes-pathway工作总结\\图\\生存")
types1=c("BLCA" ,"BRCA", "CESC" ,"CHOL" ,"COAD",
         "ESCA", "GBM" , "HNSC", "KIRC", "LIHC",
         "LUAD" ,"LUSC", "PAAD" ,"PRAD", "READ",
         "SARC", "SKCM" ,"STAD", "THCA", "UCEC")
dir = 'E:\\iraes-pathway工作总结\\irAEs6\\result2\\esetData'
file_list = list.files(path = dir, pattern = "*.csv",recursive = TRUE,full.names = TRUE)  
#f=12
for(f in 1:20){
  name = paste0("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival/survival_",types1[f],"_survival.txt")
  survival_data <- read.table(name, header = T,sep = '\t')
  exp_data <- read.csv(file_list[f], row.names = 1)
  
  if(length(grep('.',survival_data$sample))!=0)
    survival_data$sample <- gsub("-",'.',survival_data$sample)
  survival_data1 <- matrix(0, nrow = ncol(exp_data), ncol = 3)
  for(i in 1:nrow(survival_data)){
    #i=1
    n = grep(survival_data[i,1], colnames(exp_data))
    if(length(n)!=0){
      survival_data1[n,1] <- survival_data[i,1]
      survival_data1[n,2] <- survival_data[i,3]
      survival_data1[n,3] <- survival_data[i,4]
    }
  }
  colnames(survival_data1) <- c('sample', 'OS', 'OS.time')
  
  Survival <- as.numeric(survival_data1[,3])
  Events <- as.numeric(survival_data1[,2])
  PVal <- NULL
  HR <- NULL
  BB3 <- NULL
  BB4 <- NULL
  for(i in 1:nrow(exp_data)){
    #i=1
    Class <- as.numeric(exp_data[i,])
    coxph.fit <- coxph(Surv(Survival,Events) ~ Class)
    result <- summary(coxph.fit);
    aa <- result$coefficients
    PVal <- c(PVal,aa[5])
    HR <- c(HR,aa[2])
    bb <- result$conf.int
    BB3 <- c(BB3,bb[3])
    BB4 <- c(BB4,bb[4])
  }
  result <- data.frame(symbol = rownames(exp_data),PVal, HR, BB3, BB4)
  result$P_adj <- p.adjust(result$PVal, method = 'fdr')
  name1 = paste0("E:\\iraes-pathway工作总结\\图\\生存/", types1[f], "_sur_all.csv")
  write.csv(result, file = name1,row.names = F,quote=F)
  sign_result <- result[which(result[,2]<0.05),]
  name2 = paste0("E:\\iraes-pathway工作总结\\图\\生存/", types1[f], "_sur_sign.csv")
  write.csv(sign_result, file = name2,row.names = F,quote=F)
  
}
#-------------------------------------------------------------------------------------------
#统计
dir = 'E:\\iraes-pathway工作总结\\图\\生存'
file_list = list.files(path = dir, pattern = "*.csv",recursive = TRUE,full.names = TRUE) 
#f=1
for(f in 1:20){
  if(f==1){
    sign_result=import(file = file_list[f])
    tmp <- data.frame(column_name = c(types1[f], nrow(sign_result)))
    survival_result = tmp
  }else{
    #f=2
    sign_result=import( file_list[f])
    tmp <- data.frame(column_name = c(types1[f], nrow(sign_result)))
    survival_result = cbind(survival_result, tmp)
  }
  
}
colnames(survival_result)=survival_result[1,]
survival_result=survival_result[-1,]
write.csv(survival_result,file = "生存雷达.csv",row.names = F,quote = F)