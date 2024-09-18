

dir = "D:/irAEs/result2/esetData/"
file_list = list.files(path = dir, pattern = "*.csv",recursive = TRUE,full.names = TRUE)  

for(f in 1:length(file_list)){
  if(f==1){
    data = read.csv(file = file_list[f])
    marker_G1 <- rownames(termMat_G1_FC)
    term1 <- data.frame(X=rownames(termMat_G1_FC),
                        matrix(0,ncol = ncol(data),nrow = length(marker_G1)))
    colnames(term1) <- colnames(data)
    for(i in 1:nrow(term1)){
      if(length(grep(term1[i,1],data$X))!=0){
        term1[i,] <- data[grep(term1[i,1],data$X),]
      }
    }
    name=paste0("D:/irAEs/result2/uMAP/",types1[f],"_term1_uMAP.RData")
    save(term1,file = name)
    target_G1_uMAP = rep(types1[f],ncol(data)-1)
    term_G1_uMAP = term1[,-1]
  }else{
    data = read.csv(file = file_list[f])
    marker_G1 <- rownames(termMat_G1_FC)
    term1 <- data.frame(X=rownames(termMat_G1_FC),
                        matrix(0,ncol = ncol(data),nrow = length(marker_G1)))
    colnames(term1) <- colnames(data)
    for(i in 1:nrow(term1)){
      if(length(grep(term1[i,1],data$X))!=0){
        term1[i,] <- data[grep(term1[i,1],data$X),]
      }
    }
    name=paste0("D:/irAEs/result2/uMAP/",types1[f],"_term1_uMAP.RData")
    save(term1,file = name)
    target_G1_uMAP = c(target_G1_uMAP, rep(types1[f],ncol(data)-1))
    term_G1_uMAP = cbind(term_G1_uMAP, term1[,-1])

  }
}

term_G1_uMAP <- term_G1_uMAP[,!is.na(colnames(term_G1_uMAP))]

name1 = paste0("D:/irAEs/result2/",types1[f],"_target_G1_uMAP.RData")
name2 = paste0("D:/irAEs/result2/",types1[f],"_term_G1_uMAP.RData")
save(target_G1_uMAP, file = name1)
save(term_G1_uMAP, file = name2)



##G2
for(f in 1:length(file_list)){
  if(f==1){
    data = read.csv(file = file_list[f])
    marker_G2 <- rownames(termMat_G2_FC)
    term2 <- data.frame(X=rownames(termMat_G2_FC),
                        matrix(0,ncol = ncol(data),nrow = length(marker_G2)))
    colnames(term2) <- colnames(data)
    for(i in 1:nrow(term2)){
      if(length(grep(term2[i,1],data$X))!=0){
        term2[i,] <- data[grep(term2[i,1],data$X),]
      }
    }
    name=paste0("D:/irAEs/result2/uMAP/",types1[f],"_term2_uMAP.RData")
    save(term2,file = name)
    target_G2_uMAP = rep(types1[f],ncol(data)-1)
    term_G2_uMAP = term2[,-1]
  }else{
    data = read.csv(file = file_list[f])
    marker_G2 <- rownames(termMat_G2_FC)
    term2 <- data.frame(X=rownames(termMat_G2_FC),
                        matrix(0,ncol = ncol(data),nrow = length(marker_G2)))
    colnames(term2) <- colnames(data)
    for(i in 1:nrow(term2)){
      if(length(grep(term2[i,1],data$X))!=0){
        term2[i,] <- data[grep(term2[i,1],data$X),]
      }
    }
    name=paste0("D:/irAEs/result2/uMAP/",types1[f],"_term2_uMAP.RData")
    save(term2,file = name)
    target_G2_uMAP = c(target_G2_uMAP, rep(types1[f],ncol(data)-1))
    term_G2_uMAP = cbind(term_G2_uMAP, term2[,-1])
    
  }
}

term_G2_uMAP <- term_G2_uMAP[,!is.na(colnames(term_G2_uMAP))]

name1 = paste0("D:/irAEs/result2/",types1[f],"_target_G2_uMAP.RData")
name2 = paste0("D:/irAEs/result2/",types1[f],"_term_G2_uMAP.RData")
save(target_G2_uMAP, file = name1)
save(term_G2_uMAP, file = name2)

#rownames(term_G2_uMAP) <- term2[,1] 

##G3
for(f in 1:length(file_list)){
  if(f==1){
    data = read.csv(file = file_list[f])
    marker_G3 <- rownames(termMat_G3_FC)
    term3 <- data.frame(X=rownames(termMat_G3_FC),
                        matrix(0,ncol = ncol(data),nrow = length(marker_G3)))
    colnames(term3) <- colnames(data)
    for(i in 1:nrow(term3)){
      if(length(grep(term3[i,1],data$X))!=0){
        term3[i,] <- data[grep(term3[i,1],data$X),]
      }
    }
    name=paste0("D: /irAEs/result2/uMAP/",types1[f],"_term3_uMAP.RData")
    save(term3,file = name)
    target_G3_uMAP = rep(types1[f],ncol(data)-1)
    term_G3_uMAP = term3[,-1]
  }else{
    data = read.csv(file = file_list[f])
    marker_G3 <- rownames(termMat_G3_FC)
    term3 <- data.frame(X=rownames(termMat_G3_FC),
                        matrix(0,ncol = ncol(data),nrow = length(marker_G3)))
    colnames(term3) <- colnames(data)
    for(i in 1:nrow(term3)){
      if(length(grep(term3[i,1],data$X))!=0){
        term3[i,] <- data[grep(term3[i,1],data$X),]
      }
    }
    name=paste0("D: /irAEs/result2/uMAP/",types1[f],"_term3_uMAP.RData")
    save(term3,file = name)
    target_G3_uMAP = c(target_G3_uMAP, rep(types1[f],ncol(data)-1))
    term_G3_uMAP = cbind(term_G3_uMAP, term3[,-1])
    
  }
}

term_G3_uMAP <- term_G3_uMAP[,!is.na(colnames(term_G3_uMAP))]

name1 = paste0("D: /irAEs/result2/",types1[f],"_target_G3_uMAP.RData")
name2 = paste0("D: /irAEs/result2/",types1[f],"_term_G3_uMAP.RData")
save(target_G3_uMAP, file = name1)
save(term_G3_uMAP, file = name2)


target_uMAP = target_G3_uMAP
term_uMAP = rbind(term_G1_uMAP,term_G2_uMAP,term_G3_uMAP)
target_labels = factor(target_uMAP)



library(pheatmap)
library(Rtsne)
library(ggfortify)
library(mvtnorm)
set.seed(123123)
term_uMAP_T <- t(term_uMAP)
colnames(term_uMAP_T) <- 1:ncol(term_uMAP_T)
tsne_out <- Rtsne(term_uMAP_T, 
                  pca = FALSE, 
                  perplexity=10, 
                  theta=0.0,
                  check_duplicates = F)
tsnes <- tsne_out$Y
colnames(tsnes) <- c("tSNE1", "tSNE2")
tsnes <- as.data.frame(tsnes)

tsnes$group <- target_labels
color = c('#E18A7A', '#9CBCD4', '#C7DD2E', '#8BABD9', '#407DAD', 
          '#C587BB', '#C6DCC3', '#CC7741', '#B9CF47', '#8ECCD2', 
          '#89B2AC', '#DFB315', '#C67A61', '#779DE9', '#C5CA7A', 
          '#96A8C6', '#D981A4', '#DDC87B', '#69ECEB', '#826CF3')
pdf(file = "D: /irAEs/result2/tSNE_plot.pdf", width = 8, height = 7)
ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ 
  geom_point(aes(col=group))+
  scale_color_manual(values = color)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1, 
                             unit = "cm"))
dev.off()


#hub tSNE plot------------------------------------------------------------------------
dir = "D: /irAEs/result2/esetData/"
file_list = list.files(path = dir, pattern = "*.csv",recursive = TRUE,full.names = TRUE)  

for(f in 1:length(file_list)){
  if(f==1){
    data = read.csv(file = file_list[f])
    marker_hub <- hub[,1]
    term_hub <- data.frame(X = hub[,1],
                           matrix(0,ncol = ncol(data),nrow = length(marker_hub)))
    colnames(term_hub) <- colnames(data)
    for(i in 1:nrow(term_hub)){
      if(length(grep(term_hub[i,1],data$X))!=0){
        term_hub[i,] <- data[grep(term_hub[i,1],data$X),]
      }
    }
    name=paste0("D: /irAEs/result2/uMAP/",types1[f],"_term_hub_uMAP.RData")
    save(term_hub,file = name)
    target_hub_uMAP = rep(types1[f],ncol(data)-1)
    term_hub_uMAP = term_hub[,-1]
  }else{
    data = read.csv(file = file_list[f])
    marker_hub <- hub[,1]
    term_hub <- data.frame(X = hub[,1],
                           matrix(0,ncol = ncol(data),nrow = length(marker_hub)))
    colnames(term_hub) <- colnames(data)
    for(i in 1:nrow(term_hub)){
      if(length(grep(term_hub[i,1],data$X))!=0){
        term_hub[i,] <- data[grep(term_hub[i,1],data$X),]
      }
    }
    name=paste0("D: /irAEs/result2/uMAP/",types1[f],"_term_hub_uMAP.RData")
    save(term_hub,file = name)
    target_hub_uMAP = c(target_hub_uMAP, rep(types1[f],ncol(data)-1))
    term_hub_uMAP = cbind(term_hub_uMAP, term_hub[,-1])
    
  }
}

term_hub_uMAP <- term_hub_uMAP[,!is.na(colnames(term_hub_uMAP))]
target_labels = factor(target_hub_uMAP)
name1 = paste0("D: /irAEs/result2/",types1[f],"_target_hub_uMAP.RData")
name2 = paste0("D: /irAEs/result2/",types1[f],"_term_hub_uMAP.RData")
save(target_hub_uMAP, file = name1)
save(term_hub_uMAP, file = name2)



library(pheatmap)
library(Rtsne)
library(ggfortify)
library(mvtnorm)
set.seed(123123)
term_hub_uMAP_T <- t(term_hub_uMAP)
colnames(term_hub_uMAP_T) <- 1:ncol(term_hub_uMAP_T)
tsne_out <- Rtsne(term_hub_uMAP_T, 
                  pca = F, 
                  perplexity = 9, 
                  theta = 0.0,
                  check_duplicates = F)
tsnes <- tsne_out$Y
colnames(tsnes) <- c("tSNE1", "tSNE2")
tsnes <- as.data.frame(tsnes)

tsnes$group <- target_labels
color = c('#E18A7A', '#9CBCD4', '#C7DD2E', '#8BABD9', '#407DAD', 
          '#C587BB', '#C6DCC3', '#CC7741', '#B9CF47', '#8ECCD2', 
          '#89B2AC', '#DFB315', '#C67A61', '#779DE9', '#C5CA7A', 
          '#96A8C6', '#D981A4', '#DDC87B', '#69ECEB', '#826CF3')
pdf(file = "D: /irAEs/result2/tSNE_hub_plot.pdf", width = 8, height = 7)
ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ 
  geom_point(aes(col=group))+
  scale_color_manual(values = color)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1, 
                             unit = "cm"))
dev.off()




