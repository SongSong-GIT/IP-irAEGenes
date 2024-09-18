#计算Risk Score以及p值##################################################################################################################################
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(edgeR)
library(limma)##加载包
dir = "D:/卵巢癌/TCGA数据/注释为gene name的表达谱"  # 搜索指定文件夹下文件，……填写为你的文件夹路径，注意使用/做目录分隔符
#获得csv文件列表
file_list = list.files(path =dir, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)  #获得文件列表     
gmt<-read.table("D:/卵巢癌/GeneList.txt",sep = "\t",header = T)
gmt<-gmt[,c(6,1)]
colnames(gmt)<-c("ont","gene")
exp1<-read.csv(file_list[1],check.names = F,row.names = 1)
gmtlast<-c("ont","gene")
for(kk in 1:nrow(gmt)){
  if(gmt$gene[kk] %in% rownames(exp1)){
    gmtlast<-rbind(gmtlast,gmt[kk,])
  } 
}
colnames(gmtlast)<-gmtlast[1,]
gmtlast<-gmtlast[-1,]
gmt<-gmtlast
gmt<-as.data.frame(gmt)
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
  name3<-paste("D:/卵巢癌/TCGA数据/Riskscore/",name2,"_riskscore.csv",sep="")
  name4<-paste("D:/卵巢癌/TCGA数据/RiskPvalue/",name2,"_riskpvalue.csv",sep="")
  
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
    allxishu<-ES-GSEAresult2$enrichmentScore
    

    h<-intersect(rownames(exp1),zongTILncRNA[[i]])
   
    
    mat<-t(exp1_disease[gmt1$gene,])
    
    all<-0
    for(k in 1:length(allxishu)){
      all<-all+allxishu[k]*mat[,k]
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

