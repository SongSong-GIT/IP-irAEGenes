setwd(path)
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
#########################################################################


remove_last_char <- function(x) {
  # 使用substr函数截取列名的第一个字符到倒数第二个字符
  substr(x, 1, nchar(x) - 1)
}
rownames(exp) <- remove_last_char(rownames(exp))

gene=import("E:\\iraes-pathway工作总结\\图\\生存\\总的\\SKCM_sur_all.csv")
comm_gene=import("E:\\iraes-pathway工作总结\\图\\生存\\图\\comm_gene.csv")
x=intersect(gene$symbol,comm_gene$gene)
gene_symbols <- gene$symbol

# 从 data 数据框中选取列名在 gene_symbols 中的列
exp<- exp[, colnames(exp) %in% gene_symbols]
 name=rownames(exp)
exp_pre
tmp=intersect(rownames(exp),clinical_all$SKCM$sample)
clinical_all$SKCM=clinical_all$SKCM[which(clinical_all$SKCM$sample%in%tmp),]
exp=exp[clinical_all$SKCM$sample,]

exp=cbind(clinical_all$SKCM$OS.time,clinical_all$SKCM$OS,clinical_SKCM$age_at_index,clinical_SKCM$gender,exp)
colnames(exp)[1:4]=c("OS.time","OS","age","gender")
#################################################################
################################# 绘制校准曲线 ################################
library(rms)
##########################
# 打包数据
dd<-datadist(risk)
options(datadist='dd')
Cox_nomo1<-cph(Surv(OS.time, OS) ~risk,data = risk, x=T, y=T, surv=T)
# 3.生成1/3/5年的死亡概率
surv<-Survival(Cox_nomo1)
surv1 <-function(x)surv(1*365,lp=1-x) #“lp=1-x”计算死亡率，“lp=x”计算生存率
surv3 <-function(x)surv(1*1095,lp=1-x)
surv5 <-function(x)surv(1*1825,lp=1-x)

# 4.绘制nomogram
Cox_nomo1<-cph(Surv(OS.time, OS) ~age+group,data = risk, x=T, y=T, surv=T)
nomo_2a<-nomogram(Cox_nomo1, fun=list(surv1,surv3,surv5), lp=F,
                  funlabel =c("1-year Death", "3-year Death", "5-year Death"),
                  fun.at =c(0.05, seq(0.1,0.9, by=0.1), 0.95)
)
plot(nomo_2a,
     #col.grid=c("pink","cyan"),
     xfrac =0.3, #设置变量名与线段的横向占比
     
     cex.var =1, # 加粗变量字体
     
     cex.axis =1, #设置数据轴字体的代销
     
     lmgp =0.3) # 设置文字与数据轴刻度的距离
## 校准曲线
# 一年校准曲线
coxm_triplet_1 <- cph(Surv(OS.time, OS) ~risk,data=risk,
                      surv=T,x=T,y=T,time.inc = 365)
mmm<-floor(nrow(risk)/3)
cal_triplet_1<-calibrate(coxm_triplet_1,u=365,cmethod='KM',m=mmm,B=1000)

# 三年校准曲线
coxm_triplet_2 <- cph(Surv(OS.time, OS) ~risk,data=risk,
                      surv=T,x=T,y=T,time.inc = 3*365)
cal_triplet_2<-calibrate(coxm_triplet_2,u=3*365,cmethod='KM',m=mmm,B=1000)

# 五年校准曲线
coxm_triplet_3 <- cph(Surv(OS.time, OS) ~risk+group,data=risk,
                      surv=T,x=T,y=T,time.inc = 5*365)
cal_triplet_3<-calibrate(coxm_triplet_3,u=5*365,cmethod='KM',m=mmm,B=1000)

############# lnc #############
exp=exp_pre[,c("OS.time", "OS","ENSG00000234883")]
dd<-datadist(exp)
options(datadist='dd')
Cox_nomo1<-cph(Surv(OS.time, OS) ~ENSG00000234883,data = exp, x=T, y=T, surv=T)
surv<-Survival(Cox_nomo1)
surv1 <-function(x)surv(1*365,lp=1-x) #“lp=1-x”计算死亡率，“lp=x”计算生存率
surv3 <-function(x)surv(1*1095,lp=1-x)
surv5 <-function(x)surv(1*1825,lp=1-x)
coxm_lnc_1 <- cph(Surv(OS.time, OS) ~ENSG00000234883,data = exp,
                  surv=T,x=T,y=T,time.inc = 365)
mmm<-floor(nrow(risk)/3)
cal_lnc_1<-calibrate(coxm_lnc_1,u=365,cmethod='KM',m=mmm,B=1000)
coxm_lnc_2 <- cph(Surv(OS.time, OS) ~ENSG00000234883,data = exp,
                  surv=T,x=T,y=T,time.inc = 3*365)
cal_lnc_2<-calibrate(coxm_lnc_2,u=3*365,cmethod='KM',m=mmm,B=1000)
coxm_lnc_3 <- cph(Surv(OS.time, OS) ~ENSG00000234883,data = exp,
                  surv=T,x=T,y=T,time.inc = 5*365)
cal_lnc_3<-calibrate(coxm_lnc_3,u=5*365,cmethod='KM',m=mmm,B=1000)

plot(cal_triplet_1,lwd=2,lty=0, ##设置线条形状和尺寸
     errbar.col=c("#eb4d4b"), ##设置一个颜色
     xlab='Predicted Probability of 1-year OS',#便签
     ylab='Actual 1-year OS(proportion)',#标签
     col=c("#eb4d4b"),#设置一个颜色
     xlim = c(0.89,0.98),ylim = c(0.87,0.98),
     riskdist = F,subtitles = F) ##x轴和y轴范围

lines(cal_triplet_1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#eb4d4b"), pch = 16)
	  
legend("bottomright", #图例的位置
       legend = "Genes involved in survival", #图例文字
       col ="#eb4d4b", #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 0.7,#图例字体大小
       bty = "n")#不显示图例边框


# 3年
plot(cal_triplet_2,lwd=2,lty=0, ##设置线条形状和尺寸
     errbar.col=c("#eb4d4b"), ##设置一个颜色
     xlab='Predicted Probability of 3-year OS',#便签
     ylab='Actual 3-year OS(proportion)',#标签
     col=c("#eb4d4b"),#设置一个颜色
     xlim = c(0.55,0.85),ylim = c(0.55,0.85),
     riskdist = F,subtitles = F) ##x轴和y轴范围

lines(cal_triplet_2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#eb4d4b"), pch = 16)
	  
	  legend("bottomright", #图例的位置
       legend ="Genes involved in survival" , #图例文字
       col ="#eb4d4b", #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 0.7,#图例字体大小
       bty = "n")#不显示图例边框


# 5年
plot(cal_triplet_3,lwd=2,lty=0, ##设置线条形状和尺寸
     errbar.col=c("#eb4d4b"), ##设置一个颜色
     xlab='Predicted Probability of 5-year OS',#便签
     ylab='Actual 5-year OS(proportion)',#标签
     col=c("#eb4d4b"),#设置一个颜色
     xlim = c(0.4,0.8),ylim = c(0.4,0.8),
     riskdist = F,subtitles = F) ##x轴和y轴范围

lines(cal_triplet_3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#eb4d4b"), pch = 16)
	  
	  legend("bottomright", #图例的位置
       legend ="Genes involved in survival", #图例文字
       col ="#eb4d4b", #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 0.7,#图例字体大小
       bty = "n")#不显示图例边框
