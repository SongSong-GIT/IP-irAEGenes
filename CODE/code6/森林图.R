library("rio")
exp=import("E:\\iraes-pathway工作总结\\irAEs1\\exp1_disease\\TCGA-SKCM-exp1_disease.RData") 
exp=t(exp)
exp=as.data.frame(exp)
clinical_SKCM=read.csv("E:\\曲课题\\画图相关\\图六\\6A\\临床数据/clinical.tsv",header = T,sep = "\t")
clinical_SKCM=clinical_SKCM[,c(2,4,12)]
clinical_SKCM=unique(clinical_SKCM)
rownames(clinical_SKCM)=clinical_SKCM$case_submitter_id
name=substr(rownames(exp),start = 1,stop = 12)
clinical_SKCM=clinical_SKCM[name,]

##################################################################
install.packages("GSVA")
library(GSVA)
BiocManager::install("gsva")
install.packages("survminer")
library(survminer)

setwd("E:\\iraes-pathway工作总结\\图\\生存")
COAD=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_COAD_survival.txt",sep="\t",header = T)
ESCA=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_ESCA_survival.txt",sep="\t",header = T)
HNSC=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_HNSC_survival.txt",sep="\t",header = T)
KIRC=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_KIRC_survival.txt",sep="\t",header = T)
LIHC=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_LIHC_survival.txt",sep="\t",header = T)
LUAD=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_LUAD_survival.txt",sep="\t",header = T)
LUSC=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_LUSC_survival.txt",sep="\t",header = T)
PRAD=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_PRAD_survival.txt",sep="\t",header = T)
READ=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_READ_survival.txt",sep="\t",header = T)
STAD=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_STAD_survival.txt",sep="\t",header = T)
THCA=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_THCA_survival.txt",sep="\t",header = T)
UCEC=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_UCEC_survival.txt",sep="\t",header = T)
BRCA=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_BRCA_survival.txt",sep="\t",header = T)
BLCA=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_BLCA_survival.txt",sep="\t",header = T)
SKCM=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_SKCM_survival.txt",sep="\t",header = T)
CHOL=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_CHOL_survival.txt",sep="\t",header = T)
GBM=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_GBM_survival.txt",sep="\t",header = T)
CESC=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_CESC_survival.txt",sep="\t",header = T)
SARC=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_SARC_survival.txt",sep="\t",header = T)
PAAD=read.table("E:\\iraes-pathway工作总结\\irAEs3\\TCGA_survival\\survival_PAAD_survival.txt",sep="\t",header = T)

clinical_all=list("BLCA"=BLCA,"BRCA"=BRCA,"CESC"=CESC,"CHOL"=CHOL,
                  "COAD"=COAD,"ESCA"=ESCA,"GBM"=GBM,"HNSC"=HNSC,
                  "KIRC"=KIRC,"LIHC"=LIHC,"LUAD"=LUAD, "LUSC"=LUSC,
                  "PAAD"=PAAD,"PRAD"=PRAD,"READ"=READ,"SKCM"=SKCM,
                  "SARC"=SARC,"STAD"=STAD,"THCA"=THCA,"UCEC"=UCEC)
remove(COAD,BLCA,BRCA,CHOL,ESCA,GBM,HNSC,KIRC,LIHC,
       LUAD,LUSC,PAAD,CESC,SARC,PRAD,READ,SKCM,STAD,THCA,UCEC)
SKCM=clinical_all$SKCM

exp=cbind(clinical_all$SKCM$OS.time,clinical_all$SKCM$OS,clinical_SKCM$age_at_index,clinical_SKCM$gender,exp)
colnames(exp)[1:4]=c("OS.time","OS","age","gender")
exp=exp_pre
x=c("MME","TRPV4","CCL17")#ENSG00000105246(EBI3)
x=c(  "BAG3" ,   "FAM111A", "TRPV4" ,  "CCL17"  , "CD109"  , "KCTD9",   "MTCP1" ,  "PTPN14" , "TPSD1" ,  "ATL3"  ,
      "CLEC2D" , "ESYT2"  , "MYSM1"  , "P4HA3"  , "CYP4B1" , "DNASE2" , "TENT5C",  "PTPN12" , "CLEC10A" ,"CPLX1" ,
      "MME" ,   "RARRES1", "CD300E" , "SPATA25",  "BCAR3")
formula <- as.formula(paste("Surv(OS.time, OS) ~", 
                            paste(x, collapse = "+")))
fit <- coxph(formula, data = exp)

coef <- as.numeric(fit$coef) # 回归系数
score1=coef[1]*(as.numeric(exp[,x[1]]))
score2=coef[2]*(as.numeric(exp[,x[2]]))
score3=coef[3]*(as.numeric(exp[,x[3]]))
score4=coef[4]*(as.numeric(exp[,x[4]]))
score5=coef[5]*(as.numeric(exp[,x[5]]))
score6=coef[6]*(as.numeric(exp[,x[6]]))
score7=coef[7]*(as.numeric(exp[,x[7]]))
score8=coef[8]*(as.numeric(exp[,x[8]]))
score9=coef[9]*(as.numeric(exp[,x[9]]))
score10=coef[10]*(as.numeric(exp[,x[10]]))
score11=coef[11]*(as.numeric(exp[,x[11]]))
score12=coef[12]*(as.numeric(exp[,x[12]]))
score13=coef[13]*(as.numeric(exp[,x[13]]))
score14=coef[14]*(as.numeric(exp[,x[14]]))
score15=coef[15]*(as.numeric(exp[,x[15]]))
score16=coef[16]*(as.numeric(exp[,x[16]]))
score17=coef[17]*(as.numeric(exp[,x[17]]))
score18=coef[18]*(as.numeric(exp[,x[18]]))
score19=coef[19]*(as.numeric(exp[,x[19]]))
score20=coef[20]*(as.numeric(exp[,x[20]]))
score21=coef[21]*(as.numeric(exp[,x[21]]))
score22=coef[22]*(as.numeric(exp[,x[22]]))
score23=coef[23]*(as.numeric(exp[,x[23]]))
score24=coef[24]*(as.numeric(exp[,x[24]]))
score25=coef[25]*(as.numeric(exp[,x[25]]))
score=score1+score2+score3+score4+score5+score6+score7+score8+score9+score10+
  score11+score12+score13+score14+score15+score16+score17+score18+score19+
  score20+score21+score22+score23+score24+score25

train_time=clinical_all$SKCM$OS.time
train_status=clinical_all$SKCM$OS
train_group=matrix(0,2,nrow(exp))                                                                                                                  
train_group[1,]=score
colnames(train_group)=rownames(exp)
cutoff=median(train_group[1,])
index_low=train_group[1,]<=cutoff
index_high=!(index_low)
train_group[2,index_low]=1
train_group[2,index_high]=2
train=list(train_time,train_status,train_group[2,])
y=Surv(train_time,train_status)
dif=survfit(y~train_group[2,])
r=survdiff(y~train_group[2,],train)
1-pchisq(r$chisq,length(r$n)-1)
remove(fit,formula,index_high,index_low,r,score1,score2,score3,train,train_time
       ,train_status,x,y,cutoff,coef,dif)

risk=exp[,1:4]
train_group=as.data.frame(t(train_group))
risk=cbind(risk,train_group)
colnames(risk)[5:6]=c("risk","group")
risk$group=factor(risk$group,levels = c(1,2),labels = c("low","high"))
library(tidyverse)

risk=risk %>% drop_na(OS.time)
risk$age=as.numeric(risk$age)
risk <- risk[!is.na(risk$age), ]

fit <- coxph(Surv(OS.time, OS) ~age+gender+risk, data = risk)



#str(risk)

#risk$OS.time=as.character(risk$OS.time)
#risk$OS=as.character(risk$OS)


#risk <- risk[(risk$age!="'--"), ]

ggforest(fit, data = risk)



train_group=cbind(clinical_all$SKCM$OS.time,clinical_all$SKCM$OS,train_group)
colnames(train_group)=c("OS.time","OS","Risk","group")
fit <- survfit(Surv(OS.time, OS)~group,data = train_group)
## 绘制logrank的KM曲线
library(survminer)
library(survival)
ggsurvplot(fit, data = train_group,
           pval = T,conf.int = T,
           surv.median.line = "hv",
           palette = c("#E7B800","#2E9FDF"))

ggsurvplot(fit, data = train_group,
           pval = T,conf.int = T,
           surv.median.line = "hv",
           risk.table = TRUE, 
           palette = c("#E7B800","#2E9FDF"))

vioplot(Risk~group, data = train_group,names=c("low","high"),col=c("#E7B800","#2E9FDF"))
write.table(train_group,"h图.txt",sep = "\t",row.names = F,quote = F)







