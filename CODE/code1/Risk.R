#计算Risk Score以及p值##################################################################################################################################
###1.加载包
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(edgeR)
library(limma)
dir = "D:/irAEs/irAEsTCGA/"  # 搜索指定文件夹下文件，……填写为你的文件夹路径，注意使用/做目录分隔符
#获得tsv文件列表
file_list = list.files(path =dir, pattern = "*.tsv",recursive = TRUE,full.names = TRUE)  #获得文件列表     
gmt <- read.table("C:/Users/23972/Desktop/irAEs/irAEs2/GeneList.txt",sep = "\t",header = T)
gmt <- gmt[,c(6,1)]
colnames(gmt) <- c("ont","gene")
exp1 <- read.table(file_list[5], check.names = F, row.names = 1)
sample <- exp1[1,]
#将ensemble转换为gene symbol
exp1$ensembl_ID <- unlist(str_split(rownames(exp1), "[.]", simplify=T))[,1]
rownames(exp1) <- exp1$ensembl_ID
keytypes(org.Hs.eg.db) 
ensembl_ID <- rownames(exp1)
ensembl_ID <- ensembl_ID[-1]
# 采用bitr 命令进行ID的转换
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
gmtlast <- c("ont","gene")
for(kk in 1:nrow(gmt)){
  if(gmt$gene[kk] %in% rownames(exp1)){
    gmtlast <- rbind(gmtlast,gmt[kk,])
  } 
}
colnames(gmtlast) <- gmtlast[1,]
gmtlast <- gmtlast[-1,]
gmt <- gmtlast
gmt <- as.data.frame(gmt)
exp1_disease <- exp1[,grep("0[0-9]A", colnames(exp1))]
exp1_normal <- exp1[,grep("1[0-9]A", colnames(exp1))]
  if(ncol(as.matrix(exp1_normal))<2){
    next
  }
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
GSEAdata <- data.frame(gene = rownames(tT), logFC = tT$logFC)
geneList <- GSEAdata$logFC
names(geneList) = GSEAdata$gene 
geneList = sort(geneList, decreasing = T) 
RNA_score <- c("pathway", colnames(exp1_disease))
risk_p <- c("lncRNA", "pathway", "pvalue")
name1 <- unlist(strsplit(file_list[5], split = "/"))[5]
name2 <- unlist(strsplit(name1, split = "\\."))[1]
name3 <- paste("C:/Users/23972/Desktop/irAEs/irAEsResult/Riskscore/", name2, "_riskscore.csv", sep="")
name4 <- paste("C:/Users/23972/Desktop/irAEs/irAEsResult/RiskPvalue/", name2, "_riskpvalue.csv", sep="")

#for(n in 1:17){
  gmt1 <- gmt[which(gmt$ont==unique(gmt$ont)[1]),]
  gmts <- gmt1
#去掉通路i的每一个基因，gmt2是46行（去掉一个gene，第一列为去掉的是第几个gene），gmt1是原始的47行，
#gmts是47+46*47行
for(m in 1:nrow(gmt1)){
  gmt2 <- data.frame(ont=m, gene=gmt1[-m,2])
  gmts <- rbind(gmts,gmt2)
}

set.seed(1234)
GSEAresult2 <- GSEA(geneList, TERM2GENE = gmts, pvalueCutoff=1)
GSEAresult2 <- as.data.frame(GSEAresult2)
ES <- GSEAresult2$enrichmentScore[which(GSEAresult2$ID==unique(gmt$ont)[1])] #gmts中完整的47行
GSEAresult2 <- GSEAresult2[-which(GSEAresult2$ID==unique(gmt$ont)[1]),]
GSEAresult2$ID <- as.numeric(GSEAresult2$ID)
GSEAresult2 <- GSEAresult2[order(GSEAresult2$ID),]
allxishu <- ES-GSEAresult2$enrichmentScore


h<-rownames(exp1)


mat<-t(exp1_disease[gmt1$gene,])

all<-0
for(k in 1:length(allxishu)){
  all <- all+allxishu[k]*mat[,k]
}
if(length(GSEAresult2$enrichmentScore)!=0){
  RNA_score <- rbind(RNA_score,c(unique(gmt$ont)[1],all))
}else{
  RNA_score <- rbind(RNA_score,c(unique(gmt$ont)[1],rep(NA,ncol(exp1_disease)))) 
  next
}


for(j in 1:length(h)){
  group<-ifelse(as.numeric(exp1_disease[h[j],])>median(as.numeric(exp1_disease[h[j],])),"high","low")
  if(sum(group=="low")==ncol(exp1_disease)){
    next
  }
  high<-mean(all[which(group=="high")])#高表达基因们的均值
  low<-mean(all[which(group=="low")])#低表达基因们的均值
  Diff<-high-low
  d=0
  for(kk in 1:1000){
    high1<-mean(all[sample(1:ncol(exp1_disease),0.5*ncol(exp1_disease))])
    low1<-mean(all[-sample(1:ncol(exp1_disease),0.5*ncol(exp1_disease))])
    Diff1<-high1-low1
    if(abs(Diff1)>abs(Diff)){
      d<-d+1
    }
  }
  p<-d/1000
  if(length(GSEAresult2$enrichmentScore)!=0){
    risk_p<-rbind(risk_p,c(h[j],unique(gmt$ont)[1],p))
  }else{
    risk_p<-rbind(risk_p,c(h[j],unique(gmt$ont)[1],NA))
  }
}





for(n in 1:17){
  gmt1 <- gmt[which(gmt$ont==unique(gmt$ont)[n]),]
  gmts <- gmt1
  for(m in 1:nrow(gmt1)){
    gmt2 <- data.frame(ont=m, gene=gmt1[-m,2])
    gmts <- rbind(gmts,gmt2)
  }
  
  set.seed(1234)
  GSEAresult2<-GSEA(geneList, TERM2GENE = gmts,pvalueCutoff=1)
  GSEAresult2<-as.data.frame(GSEAresult2)
  ES<-GSEAresult2$enrichmentScore[which(GSEAresult2$ID==unique(gmt$ont)[n])]
  GSEAresult2<-GSEAresult2[-which(GSEAresult2$ID==unique(gmt$ont)[n]),]
  GSEAresult2$ID<-as.numeric(GSEAresult2$ID)
  GSEAresult2<-GSEAresult2[order(GSEAresult2$ID),]
  allxishu<-ES-GSEAresult2$enrichmentScore
  
  
  h<-rownames(exp1)
  
  
  mat<-t(exp1_disease[gmt1$gene,])
  
  all<-0
  for(k in 1:length(allxishu)){
    all<-all+allxishu[k]*mat[,k]
  }
  if(length(GSEAresult2$enrichmentScore)!=0){
    RNA_score<-rbind(RNA_score,c(unique(gmt$ont)[n],all))
  }else{
    RNA_score<-rbind(RNA_score,c(unique(gmt$ont)[n],rep(NA,ncol(exp1_disease)))) 
  }
  
  if(length(GSEAresult2$enrichmentScore)==0){
    next
  }
  for(j in 1:length(h)){
    group<-ifelse(as.numeric(exp1_disease[h[j],])>median(as.numeric(exp1_disease[h[j],])),"high","low")
    if(sum(group=="low")==ncol(exp1_disease)){
      next
    }
    high<-mean(all[which(group=="high")])
    low<-mean(all[which(group=="low")])
    Diff<-high-low
    d=0
    for(kk in 1:1000){
      high1<-mean(all[sample(1:ncol(exp1_disease),0.5*ncol(exp1_disease))])
      low1<-mean(all[-sample(1:ncol(exp1_disease),0.5*ncol(exp1_disease))])
      Diff1<-high1-low1
      if(abs(Diff1)>abs(Diff)){
        d<-d+1
      }
    }
    p<-d/1000
    if(length(GSEAresult2$enrichmentScore)!=0){
      risk_p<-rbind(risk_p,c(h[j],unique(gmt$ont)[n],p))
    }else{
      risk_p<-rbind(risk_p,c(h[j],unique(gmt$ont)[n],NA))
    }
  }
  
}
write.csv(lnc_score,name3,row.names = F)
write.csv(risk_p,name4,row.names = F)
print(c(i,j))


###计算每个基因的差异得分Diff
diffScoreFunction <- function(j){
  
  geneValue <- as.numeric(exp1_disease[j,])
  group <- ifelse(geneValue > median(geneValue), "high", "low")
  if(sum(group=="low")==ncol(exp1_disease)){next}
  high <- mean(all[which(group=="high")])  #all是∑[δES*exp],一个样本对应一个all值，即所有基因在一个样本的得分
  low <- mean(all[which(group=="low")])  #基因j的高(低)表达样本的表达均值，
  Diff <- high-low  #基因j的表达差异得分Diff
  
  ######随机部分
  ###随机部分函数###
  permutationFunction <- function(){
    d = 0
    high1 <- mean(all[sample(1:ncol(exp1_disease),0.5*ncol(exp1_disease))]) #从所有样本对应的ES*exp中抽取一半
    low1 <- mean(all[-sample(1:ncol(exp1_disease),0.5*ncol(exp1_disease))]) #剩下的一半作为低表达
    Diff1 <- high1-low1 #基因j的随机表达差异得分Diff1
    if(abs(Diff1) > abs(Diff)){d <- d+1}
    return(d)
  }
  ###随机部分实际###
  library(foreach)
  library(doParallel)
  no_cores1 <- detectCores() - 1
  cl1 <- makeCluster(no_cores1)
  registerDoParallel(cl1)
  times = 1000 #随机1000次
  permuteResult <- foreach(i = seq(times),
                           .combine = rbind)  %dopar%  permutationFunction()
  d <- length(which(permuteResult==1)) 
  p <- d/1000
  if(length(GSEAresult2$enrichmentScore) != 0){
    risk_p <- c(risk_p, c(j, unique(gmt$ont)[1], p)) #基因名+路径名+p值
  }else{
    risk_p <- c(risk_p, c(j, unique(gmt$ont)[1], NA)) #基因名+路径名+p值(NA)
  }
  stopCluster(cl1)
  ######
  
  return(risk_p)
}


library(foreach)
library(doParallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)
geneName <- rownames(exp1)
diffScoreResult <- foreach(j = geneName[1:10],
                         .combine = rbind)  %dopar%  diffScoreFunction(j)

write.csv(RNA_score, name3, row.names = F)
write.csv(risk_p, name4, row.names = F)
print(c(i,j))
stopCluster(cl)
####










###
permutationFunction <- function(){
  d = 0
  high1 <- mean(all[sample(1:ncol(exp1_disease),0.5*ncol(exp1_disease))]) #从所有样本对应的ES*exp中抽取一半
  low1 <- mean(all[-sample(1:ncol(exp1_disease),0.5*ncol(exp1_disease))]) #剩下的一半作为低表达
  Diff1 <- high1-low1 #随机差异Diff1
  if(abs(Diff1) > abs(Diff)){d <- d+1}
  return(d)
}

library(foreach)
library(doParallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)
times = 1000 #随机1000次
permuteResult <- foreach(i = seq(times),
                     .combine = rbind)  %dopar%  permutationFunction()
d <- length(which(permuteResult==1)) 
p <- d/1000
if(length(GSEAresult2$enrichmentScore) != 0){
  risk_p <- rbind(risk_p, c(h[j], unique(gmt$ont)[n], p)) #基因名+路径名+p值
}else{
  risk_p <- rbind(risk_p, c(h[j], unique(gmt$ont)[n], NA)) #基因名+路径名+p值(NA)
}
stopCluster(cl)
####










for(i in 1:length(file_list)){
  exp1<-read.csv(file_list[i],check.names = F,row.names = 1)
  exp1<-exp1[,-ncol(exp1)]
  exp1_disease<-exp1[,grep("01A",colnames(exp1))]
  exp1_normal<-exp1[,grep("11A",colnames(exp1))]
  if(ncol(as.matrix(exp1_normal))<2){
    next
  }
  expall<-cbind(exp1_disease,exp1_normal)
  group<-c(rep("disease",ncol(exp1_disease)),rep("normal",ncol(exp1_normal)))
  eset=expall ##将基因表达矩阵赋值给eset
  targets<-cbind(colnames(expall),group)##读入样本数据，包括两列，
  #第一列GSM号，第二列为样本分组
  targets<-as.data.frame(targets)
  colnames(targets)=c("FileName","Target")#更改列名，为了和limma包中的一致
  lev<-unique(targets$Target)
  design <- model.matrix(~0+factor(targets$Target, levels=lev)) #样本矩阵
  colnames(design) <- lev #更改列名为levels名
  ###两两之间的比较，求差异基因使用topTable函数
  cont.wt <- makeContrasts(disease-normal,
                           levels=design) 
  
  fit <- lmFit(eset, design)
  fit2 <- contrasts.fit(fit, cont.wt) 
  fit2 <- eBayes(fit2) 
  tT=topTable(fit2,adjust="fdr",n=Inf)
  tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
  colnames(tT)=c("FDR","P.Value","logFC")
  GSEAdata<-data.frame(gene=rownames(tT),logFC=tT$logFC)
  geneList<-GSEAdata$logFC#第二列可以是folodchange，也可以是logFC
  names(geneList)=GSEAdata$gene #使用转换好的ID
  geneList=sort(geneList,decreasing = T) #从高到低排序
  
  lnc_score<-c("pathway",colnames(exp1_disease))
  risk_p<-c("lncRNA","pathway","pvalue")
  name1<-unlist(strsplit(file_list[i],split = "/"))[5]
  name2<-unlist(strsplit(name1,split = "\\."))[1]
  name3<-paste("D:/TCGA数据/Riskscore/",name2,"_riskscore.csv",sep="")
  name4<-paste("D:/TCGA数据/RiskPvalue/",name2,"_riskpvalue.csv",sep="")
  
  for(n in 1:17){
    gmt1<-gmt[which(gmt$ont==unique(gmt$ont)[n]),]
    gmts<-gmt1
    for(m in 1:nrow(gmt1)){
      gmt2<-data.frame(ont=m,gene=gmt1[-m,2])
      gmts<-rbind(gmts,gmt2)
    }
    
    set.seed(1234)
    GSEAresult2<-GSEA(geneList,TERM2GENE = gmts,pvalueCutoff=1)
    GSEAresult2<-as.data.frame(GSEAresult2)
    ES<-GSEAresult2$enrichmentScore[which(GSEAresult2$ID==unique(gmt$ont)[n])]
    GSEAresult2<-GSEAresult2[-which(GSEAresult2$ID==unique(gmt$ont)[n]),]
    GSEAresult2$ID<-as.numeric(GSEAresult2$ID)
    GSEAresult2<-GSEAresult2[order(GSEAresult2$ID),]
    allxishu<-ES-GSEAresult2$enrichmentScore  #allxishu是每个gene的δES 306*1,allxishu[k]是1个值
    
    
    h<-intersect(rownames(exp1),zongTILncRNA[[i]])
    
    
    mat<-t(exp1_disease[gmt1$gene,]) #mat[,k]是每个gene在所有样本的表达值 36*1
                                     #allxishu[k]*mat[,k]是每个gene在所有样本的表达值*δES
    all<-0
    for(k in 1:length(allxishu)){
      all <- all+allxishu[k]*mat[,k]  #全部gene在所有样本的表达值*δES 36*1
    }
    if(length(GSEAresult2$enrichmentScore)!=0){
      lnc_score<-rbind(lnc_score,c(unique(gmt$ont)[n],all))
    }else{
      lnc_score<-rbind(lnc_score,c(unique(gmt$ont)[n],rep(NA,ncol(exp1_disease)))) 
    }
    
    if(length(GSEAresult2$enrichmentScore)==0){
      next
    }
    for(j in 1:length(h)){
      group<-ifelse(as.numeric(exp1_disease[h[j],])>median(as.numeric(exp1_disease[h[j],])),"high","low")
      if(sum(group=="low")==ncol(exp1_disease)){
        next
      }
      high<-mean(all[which(group=="high")])
      low<-mean(all[which(group=="low")])
      Diff<-high-low
      d=0
      for(kk in 1:1000){
        high1<-mean(all[sample(1:ncol(exp1_disease),0.5*ncol(exp1_disease))])
        low1<-mean(all[-sample(1:ncol(exp1_disease),0.5*ncol(exp1_disease))])
        Diff1<-high1-low1
        if(abs(Diff1)>abs(Diff)){
          d<-d+1
        }
      }
      p<-d/1000
      if(length(GSEAresult2$enrichmentScore)!=0){
        risk_p<-rbind(risk_p,c(h[j],unique(gmt$ont)[n],p))
      }else{
        risk_p<-rbind(risk_p,c(h[j],unique(gmt$ont)[n],NA))
      }
    }
    
  }
  write.csv(lnc_score,name3,row.names = F)
  write.csv(risk_p,name4,row.names = F)
  print(c(i,j))
}

