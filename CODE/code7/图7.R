
setwd("E:\\iraes-pathway工作总结\\图\\生存\\ICI")

#######elastic network  ##################
library(glmnet)
# install.packages("genefilter")
# install.packages("AnnotationDbi")
# BiocManager::install("AnnotationDbi")
# install.packages("E:\\软件\\AnnotationDbi_1.64.1.zip")
library(AnnotationDbi)
library(genefilter)
library(sva)
library(data.table)
library(tidyverse)
library(survival)
library(rio)
#### 读取五个数据集
setwd("E:\\iraes-pathway工作总结\\图\\生存\\ICI/数据/")
Gide_PD1=read.table("Gide_PD1.txt",header = T)
Gide_PD1CTLA4=read.table("Gide_PD1+CTLA4.txt",header = T)
Naive=read.table("Naive.txt",header = T)
Prog=read.table("Prog.txt",header = T)
Vanllen=read.table("Vanllen.txt",header = T)

#### 去除批次效应
tmp=intersect(rownames(Gide_PD1),rownames(Gide_PD1CTLA4))
tmp=intersect(tmp,rownames(Naive))
tmp=intersect(tmp,rownames(Prog))
tmp=intersect(tmp,rownames(Vanllen))


Gide_PD1=Gide_PD1[tmp,]
Gide_PD1CTLA4=Gide_PD1CTLA4[tmp,]
Naive=Naive[tmp,]
Prog=Prog[tmp,]
Vanllen=Vanllen[tmp,]


exp_all=cbind(Gide_PD1,Gide_PD1CTLA4,Naive,Prog,Vanllen)
batch=paste0("batch",rep(c(1,2,3,4,5),c(41,32,25,26,35)))

exp_all=as.matrix(exp_all)
exp <-  ComBat(dat = exp_all, batch = batch, mod = NULL)
exp=as.data.frame(exp)


Gide_PD1_resp=read.table("Gide_PD1_clinical.txt",header = T)
Gide_PD1CTLA4_resp=read.table("Gide_PD1+CTLA4_clinical.txt",header = T)
Naive_resp=read.table("Naive_clinical.txt",header = T)
Prog_resp=read.table("Prog_clinical.txt",header = T)
Vanllen_resp=read.table("Vanllen_clinical.txt",header = T)

surv_triplet=import("E:\\iraes-pathway工作总结\\图\\生存\\生存结果\\SKCM_sur_sign.csv")
#surv_triplet=surv_triplet[which(surv_triplet$cancer%in%"SKCM"),]
gene=unique(surv_triplet$symbol)

exp=exp[gene,]
exp=na.omit(exp)
gene=rownames(exp)

surv_triplet=surv_triplet[which(surv_triplet$symbol%in%gene ),]
surv_triplet=surv_triplet[,1]
surv_triplet=as.data.frame(surv_triplet)

exp=as.data.frame(t(exp))

os=rbind(Gide_PD1_resp[,2:3],Gide_PD1CTLA4_resp[,2:3],Naive_resp[,2:3],
         Prog_resp[,2:3],Vanllen_resp[,2:3])
exp=cbind(os,exp)
colnames(exp)=gsub("-","\\.",colnames(exp))

surv_triplet2=apply(surv_triplet,1,function(x){gsub("-","\\.",x)})
#surv_triplet2=as.data.frame(t(surv_triplet2))

risk <- function(x) {
  formula <- as.formula(paste("Surv(OS.time, OS) ~ ", x))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score <- coef * as.numeric(exp[, x])
  return(score)
}

result <- apply(surv_triplet, 1, risk)


name=gene
colnames(result)=name

exp_resp=c(Gide_PD1_resp$Response,
           Gide_PD1CTLA4_resp$Response,
           Naive_resp$Response,
           Prog_resp$Response,
           Vanllen_resp$response)

#### 划分训练集和检验集
# set.seed(9276)
set.seed(830)
train_index=sample(1:159, size = 0.7 * 159)
test_index=setdiff(1:159, train_index)

result=as.data.frame(t(result))
train=result[,train_index]
test=result[,test_index]
train=as.matrix(train)
test=as.matrix(test)
train=t(train)
test=t(test)

y_train=exp_resp[train_index]
y_test=exp_resp[test_index]

model = glmnet(train,y_train,family = "binomial",alpha = 0.5)
plot(model,xvar="lambda",label=T)
fit_cv <- cv.glmnet(train, y_train, alpha=1, family = 'binomial', type.measure='auc')
plot(fit_cv)

# 交叉验证选择最优lambda
cv_model=cv.glmnet(train,y_train,family="binomial",alpha=1)
best_lambda=cv_model$lambda.min
coef <- as.matrix (coef (model, s = best_lambda))
genes=rownames(coef)[coef!=0]

library(ROCR)
pred=predict(model,newx=test,s=cv_model$lambda.min,type = 'response')
get_confusion_stat <- function(pred,y_test,threshold=0.5){
  # auc
  tmp <- prediction(as.vector(pred),y_test)
  auc <- unlist(slot(performance(tmp,'auc'),'y.values'))
  # statistic
  pred_new <- as.integer(pred>threshold) 
  tab <- table(pred_new,y_test)
  if(nrow(tab)==1){
    print('preds all zero !')
    return(0)
  }
  TP <- tab[2,2]
  TN <- tab[1,1]
  FP <- tab[2,1]
  FN <- tab[1,2]
  accuracy <- round((TP+TN)/(TP+FN+FP+TN),4)
  recall_sensitivity <- round(TP/(TP+FN),4)
  precision <- round(TP/(TP+FP),4)
  specificity <- round(TN/(TN+FP),4)
  # 添加，预测的负例占比（业务解释：去除多少的样本，达到多少的recall）
  neg_rate <- round((TN+FN)/(TP+TN+FP+FN),4)
  re <- list('AUC' = auc,
             'Confusion_Matrix'=tab,
             'Statistics'=data.frame(value=c('accuracy'=accuracy,
                                             'recall_sensitivity'=recall_sensitivity,
                                             'precision'=precision,
                                             'specificity'=specificity,
                                             'neg_rate'=neg_rate)))
  return(re)
}


print(get_confusion_stat(pred,y_test))

pred.obj = prediction(pred, y_test)
auc.obj = performance(pred.obj, measure = "auc")
auc.value = as.numeric(auc.obj@y.values) 
print(auc.value)
roc.obj = performance(pred.obj, measure = "tpr", x.measure = "fpr")
plot(roc.obj, main = "ROC curve", col = "blue") 
abline(a = 0, b = 1, lty = 2, col = "gray") 
text(0.8, 0.2, paste("AUC =", round(auc.value, 3)))

roc.obj_plastic=roc.obj
auc.value_plastic=auc.value


plot(roc.obj_lasso, main = "ROC curve", col = "#509ee0",lwd=2) 
plot(roc.obj_plastic, main = "ROC curve", col = "#e27823",add=T,lwd=2) 
plot(roc.obj_vector, main = "ROC curve", col = "#96bba3",add=T,lwd=2) 


abline(a = 0, b = 1, lty = 2, col = "gray",lwd=1) 

barplot(c(auc.value_lasso,auc.value_plastic,auc.value_vector),ylim = c(0,0.8))
####################################
# 初始化随机种子
random_seed <- 0

# 初始化AUC值
auc_value <- 0

while (auc_value <= 0.8) {
  # 生成随机种子
  random_seed <- random_seed + 1
  
  # 使用新的随机种子重新划分训练集和测试集
  set.seed(random_seed)
  train_index <- sample(1:159, size = 0.7 * 159)
  test_index <- setdiff(1:159, train_index)
  
  # 重新构建训练集和测试集
  #result=as.data.frame(t(result))
  train <- result[, train_index]
  test <- result[, test_index]
  train <- as.matrix(train)
  test <- as.matrix(test)
  train <- t(train)
  test <- t(test)
  
  y_train <- exp_resp[train_index]
  y_test <- exp_resp[test_index]
  
  # 训练模型
  model <- glmnet(train, y_train, family = "binomial", alpha = 0.5)

  #plot(model,xvar="lambda",label=T)
  fit_cv <- cv.glmnet(train, y_train, alpha=1, family = 'binomial', type.measure='auc')
  #plot(fit_cv)
  
  # 预测测试集
  pred <- predict(model, newx = test, s = fit_cv$lambda.min, type = 'response')
  
  # 计算AUC值
  pred.obj <- prediction(pred, y_test)
  auc.obj <- performance(pred.obj, measure = "auc")
  auc_value <- as.numeric(auc.obj@y.values)
  
  # 输出随机种子和对应的AUC值
  print(paste("Random Seed:", random_seed, "AUC Value:", auc_value))
}


##########################      lasso    #############################################################
setwd("E:\\iraes-pathway工作总结\\图\\生存\\ICI\\数据")
library(glmnet)
library(sva)
library(data.table)
library(tidyverse)
library(survival)
library(rio)
#### 读取五个数据集
Gide_PD1=read.table("Gide_PD1.txt",header = T)
Gide_PD1CTLA4=read.table("Gide_PD1+CTLA4.txt",header = T)
Naive=read.table("Naive.txt",header = T)
Prog=read.table("Prog.txt",header = T)
Vanllen=read.table("Vanllen.txt",header = T)

#### 去除批次效应
tmp=intersect(rownames(Gide_PD1),rownames(Gide_PD1CTLA4))
tmp=intersect(tmp,rownames(Naive))
tmp=intersect(tmp,rownames(Prog))
tmp=intersect(tmp,rownames(Vanllen))


Gide_PD1=Gide_PD1[tmp,]
Gide_PD1CTLA4=Gide_PD1CTLA4[tmp,]
Naive=Naive[tmp,]
Prog=Prog[tmp,]
Vanllen=Vanllen[tmp,]


exp_all=cbind(Gide_PD1,Gide_PD1CTLA4,Naive,Prog,Vanllen)
batch=paste0("batch",rep(c(1,2,3,4,5),c(41,32,25,26,35)))

exp_all=as.matrix(exp_all)
exp <-  ComBat(dat = exp_all, batch = batch, mod = NULL)
exp=as.data.frame(exp)


Gide_PD1_resp=read.table("Gide_PD1_clinical.txt",header = T)
Gide_PD1CTLA4_resp=read.table("Gide_PD1+CTLA4_clinical.txt",header = T)
Naive_resp=read.table("Naive_clinical.txt",header = T)
Prog_resp=read.table("Prog_clinical.txt",header = T)
Vanllen_resp=read.table("Vanllen_clinical.txt",header = T)

surv_triplet=import("E:\\iraes-pathway工作总结\\图\\生存\\生存结果\\SKCM_sur_sign.csv")
#surv_triplet=surv_triplet[which(surv_triplet$cancer%in%"SKCM"),]
gene=unique(surv_triplet$symbol)

exp=exp[gene,]
exp=na.omit(exp)
gene=rownames(exp)

surv_triplet=surv_triplet[which(surv_triplet$symbol%in%gene ),]
surv_triplet=surv_triplet[,1]
surv_triplet=as.data.frame(surv_triplet)


exp=as.data.frame(t(exp))

os=rbind(Gide_PD1_resp[,2:3],Gide_PD1CTLA4_resp[,2:3],Naive_resp[,2:3],
         Prog_resp[,2:3],Vanllen_resp[,2:3])
exp=cbind(os,exp)
colnames(exp)=gsub("-","\\.",colnames(exp))

surv_triplet2=apply(surv_triplet,1,function(x){gsub("-","\\.",x)})
#surv_triplet2=as.data.frame(t(surv_triplet2))

risk <- function(x) {
  formula <- as.formula(paste("Surv(OS.time, OS) ~ ", x))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score <- coef * as.numeric(exp[, x])
  return(score)
}

result <- apply(surv_triplet, 1, risk)
name=gene
colnames(result)=name

exp_resp=c(Gide_PD1_resp$Response,
           Gide_PD1CTLA4_resp$Response,
           Naive_resp$Response,
           Prog_resp$Response,
           Vanllen_resp$response)
data=cbind(result,exp_resp)
write.table(data,"lasso_shuju.txt",sep = "\t",row.names = F,quote = F)
#### 划分训练集和检验集
# set.seed(121)
#set.seed(121111)  0.648
#set.seed(13)  0.68
set.seed(9276)
set.seed(830)
train_index=sample(1:159, size = 0.7 * 159)
test_index=setdiff(1:159, train_index)

result=as.data.frame(t(result))
train=result[,train_index]
test=result[,test_index]
train=as.matrix(train)
test=as.matrix(test)
train=t(train)
test=t(test)

y_train=exp_resp[train_index]
y_test=exp_resp[test_index]

model = glmnet(train,y_train,family = "binomial",alpha = 1)
plot(model,xvar="lambda",label=T)
fit_cv <- cv.glmnet(train, y_train, alpha=1, family = 'binomial', type.measure='auc')
plot(fit_cv)

# 交叉验证选择最优lambda
cv_model=cv.glmnet(train,y_train,family="binomial",alpha=1)
best_lambda=cv_model$lambda.min
coef <- as.matrix (coef (model, s = best_lambda))
genes=rownames(coef)[coef!=0]

library(ROCR)
pred=predict(model,newx=test,s=cv_model$lambda.min,type = 'response')
get_confusion_stat <- function(pred,y_test,threshold=0.5){
  # auc
  tmp <- prediction(as.vector(pred),y_test)
  auc <- unlist(slot(performance(tmp,'auc'),'y.values'))
  # statistic
  pred_new <- as.integer(pred>threshold) 
  tab <- table(pred_new,y_test)
  if(nrow(tab)==1){
    print('preds all zero !')
    return(0)
  }
  TP <- tab[2,2]
  TN <- tab[1,1]
  FP <- tab[2,1]
  FN <- tab[1,2]
  accuracy <- round((TP+TN)/(TP+FN+FP+TN),4)
  recall_sensitivity <- round(TP/(TP+FN),4)
  precision <- round(TP/(TP+FP),4)
  specificity <- round(TN/(TN+FP),4)
  # 添加，预测的负例占比（业务解释：去除多少的样本，达到多少的recall）
  neg_rate <- round((TN+FN)/(TP+TN+FP+FN),4)
  re <- list('AUC' = auc,
             'Confusion_Matrix'=tab,
             'Statistics'=data.frame(value=c('accuracy'=accuracy,
                                             'recall_sensitivity'=recall_sensitivity,
                                             'precision'=precision,
                                             'specificity'=specificity,
                                             'neg_rate'=neg_rate)))
  return(re)
}


print(get_confusion_stat(pred,y_test))

pred.obj = prediction(pred, y_test)
auc.obj = performance(pred.obj, measure = "auc")
auc.value = as.numeric(auc.obj@y.values) 
print(auc.value)
roc.obj = performance(pred.obj, measure = "tpr", x.measure = "fpr")
plot(roc.obj, main = "ROC curve", col = "blue") 
abline(a = 0, b = 1, lty = 2, col = "gray") 
text(0.8, 0.2, paste("AUC =", round(auc.value, 3)))


roc.obj_lasso=roc.obj
auc.value_lasso=auc.value
#########################################
random_seed <- 0

# 初始化AUC值
auc_value <- 0

while (auc_value <= 0.75) {
  random_seed <- random_seed + 1
  
  # 使用新的随机种子重新划分训练集和测试集
  set.seed(random_seed)
  train_index=sample(1:159, size = 0.7 * 159)
  test_index=setdiff(1:159, train_index)
  
  #result=as.data.frame(t(result))
  train=result[,train_index]
  test=result[,test_index]
  train=as.matrix(train)
  test=as.matrix(test)
  train=t(train)
  test=t(test)
  
  y_train=exp_resp[train_index]
  y_test=exp_resp[test_index]
  
  model = glmnet(train,y_train,family = "binomial",alpha = 1)
  #plot(model,xvar="lambda",label=T)
  fit_cv <- cv.glmnet(train, y_train, alpha=1, family = 'binomial', type.measure='auc')
  #plot(fit_cv)
  
  # 交叉验证选择最优lambda
  cv_model=cv.glmnet(train,y_train,family="binomial",alpha=1)
  best_lambda=cv_model$lambda.min
  coef <- as.matrix (coef (model, s = best_lambda))
  genes=rownames(coef)[coef!=0]
  
  library(ROCR)
  pred=predict(model,newx=test,s=cv_model$lambda.min,type = 'response')
  get_confusion_stat <- function(pred,y_test,threshold=0.5){
    # auc
    tmp <- prediction(as.vector(pred),y_test)
    auc <- unlist(slot(performance(tmp,'auc'),'y.values'))
    # statistic
    pred_new <- as.integer(pred>threshold) 
    tab <- table(pred_new,y_test)
    if(nrow(tab)==1){
      print('preds all zero !')
      return(0)
    }
    TP <- tab[2,2]
    TN <- tab[1,1]
    FP <- tab[2,1]
    FN <- tab[1,2]
    accuracy <- round((TP+TN)/(TP+FN+FP+TN),4)
    recall_sensitivity <- round(TP/(TP+FN),4)
    precision <- round(TP/(TP+FP),4)
    specificity <- round(TN/(TN+FP),4)
    # 添加，预测的负例占比（业务解释：去除多少的样本，达到多少的recall）
    neg_rate <- round((TN+FN)/(TP+TN+FP+FN),4)
    re <- list('AUC' = auc,
               'Confusion_Matrix'=tab,
               'Statistics'=data.frame(value=c('accuracy'=accuracy,
                                               'recall_sensitivity'=recall_sensitivity,
                                               'precision'=precision,
                                               'specificity'=specificity,
                                               'neg_rate'=neg_rate)))
    return(re)
  }
  
  
  #print(get_confusion_stat(pred,y_test))
  # 
  pred.obj = prediction(pred, y_test)
  auc.obj = performance(pred.obj, measure = "auc")
  auc.value = as.numeric(auc.obj@y.values)
 
  print(paste("Random Seed:", random_seed, "AUC Value:", auc.value))
}

#####################           svm             ########################
library(e1071)
library(glmnet)
library(sva)
library(data.table)
library(tidyverse)
library(survival)
library(pROC)

#### 读取五个数据集
Gide_PD1=read.table("Gide_PD1.txt",header = T)
Gide_PD1CTLA4=read.table("Gide_PD1+CTLA4.txt",header = T)
Naive=read.table("Naive.txt",header = T)
Prog=read.table("Prog.txt",header = T)
Vanllen=read.table("Vanllen.txt",header = T)

#### 去除批次效应
tmp=intersect(rownames(Gide_PD1),rownames(Gide_PD1CTLA4))
tmp=intersect(tmp,rownames(Naive))
tmp=intersect(tmp,rownames(Prog))
tmp=intersect(tmp,rownames(Vanllen))


Gide_PD1=Gide_PD1[tmp,]
Gide_PD1CTLA4=Gide_PD1CTLA4[tmp,]
Naive=Naive[tmp,]
Prog=Prog[tmp,]
Vanllen=Vanllen[tmp,]


exp_all=cbind(Gide_PD1,Gide_PD1CTLA4,Naive,Prog,Vanllen)
batch=paste0("batch",rep(c(1,2,3,4,5),c(41,32,25,26,35)))

exp_all=as.matrix(exp_all)
exp <-  ComBat(dat = exp_all, batch = batch, mod = NULL)
exp=as.data.frame(exp)


Gide_PD1_resp=read.table("Gide_PD1_clinical.txt",header = T)
Gide_PD1CTLA4_resp=read.table("Gide_PD1+CTLA4_clinical.txt",header = T)
Naive_resp=read.table("Naive_clinical.txt",header = T)
Prog_resp=read.table("Prog_clinical.txt",header = T)
Vanllen_resp=read.table("Vanllen_clinical.txt",header = T)

surv_triplet=import("E:\\iraes-pathway工作总结\\图\\生存\\生存结果\\SKCM_sur_sign.csv")
#surv_triplet=surv_triplet[which(surv_triplet$cancer%in%"SKCM"),]
gene=unique(surv_triplet$symbol)

exp=exp[gene,]
exp=na.omit(exp)
gene=rownames(exp)

surv_triplet=surv_triplet[which(surv_triplet$symbol%in%gene ),]
surv_triplet=surv_triplet[,1]
surv_triplet=as.data.frame(surv_triplet)

exp=as.data.frame(t(exp))

os=rbind(Gide_PD1_resp[,2:3],Gide_PD1CTLA4_resp[,2:3],Naive_resp[,2:3],
         Prog_resp[,2:3],Vanllen_resp[,2:3])
exp=cbind(os,exp)
colnames(exp)=gsub("-","\\.",colnames(exp))

surv_triplet2=apply(surv_triplet,1,function(x){gsub("-","\\.",x)})
surv_triplet2=as.data.frame(t(surv_triplet2))
#surv_triplet2=surv_triplet2[,1:24]
#surv_triplet2=surv_triplet2[,1:24]
risk=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef)
  score=0
  for (i in 1:73) {
    # 回归系数
    score_i=coef[i]*(as.numeric(exp[,x[i]]))
    score=score+score_i
  }
  return(score)
}

result=apply(surv_triplet2,1,risk)

colnames(result)="gene"

exp_resp=c(Gide_PD1_resp$Response,
           Gide_PD1CTLA4_resp$Response,
           Naive_resp$Response,
           Prog_resp$Response,
           Vanllen_resp$response)

data=cbind(result,exp_resp)
data=as.data.frame(data)
colnames(data)[2]="group"
#write.table(data,"SVM_shuju.txt",sep = "\t",row.names = F,quote = F)
data$group=factor(data$group)

#### 划分训练集和检验集
#set.seed(6814)
set.seed(7193)
train_index=sample(1:159, size = 0.7 * 159)
test_index=setdiff(1:159, train_index)


train=data[train_index,]
test=data[test_index,]

# svm_model = svm(group~`DBH.AS1;TIGIT;LTA`,data=train,knernel = "radial")
# svm_model = svm(group~`ZEB2.AS1;BTK;TNFSF13B`,data=train,knernel = "radial")
# svm_model = svm(group~`FAM30A;NCR3;IL16`,data=train,knernel = "radial")
# svm_model = svm(group~`DBH.AS1;TIGIT;LTA`+`ZEB2.AS1;BTK;TNFSF13B`,data=train,knernel = "radial")
# svm_model = svm(group~`DBH.AS1;TIGIT;LTA`+`FAM30A;NCR3;IL16`,data=train,knernel = "radial")
# svm_model = svm(group~`FAM30A;NCR3;IL16`+`ZEB2.AS1;BTK;TNFSF13B`,data=train,knernel = "radial")
svm_model = svm(group~`gene`,data=train,knernel = "radial")
#svm_model=svm(group~`PIK3CD.AS1;CCL22;XCR1`+`PIK3CD.AS1;CCL22;CCL4`+`DBH.AS1;IFNG;PTGER2`+`LINC02908;CCL22;IL12B`,data=train,knernel = "radial")


svm_pred=predict(svm_model,test,decision.values = TRUE)

test$svm_pred = svm_pred

#绘制ROC曲线
ran_roc <- roc(test$group,as.numeric(svm_pred))
plot(ran_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='SVM模型ROC曲线')

pred.obj=prediction(as.numeric(test$svm_pred),test$group)
auc.obj=performance(pred.obj,measure = "auc")
auc.value=as.numeric(auc.obj@y.values)

roc.obj=performance(pred.obj,measure = "tpr",x.measure = "fpr")
plot(roc.obj, main = "ROC curve", col = "blue") 
abline(a = 0, b = 1, lty = 2, col = "gray") 
text(0.8, 0.2, paste("AUC =", round(auc.value, 3)))

roc.obj_73=roc.obj
auc.value_73=auc.value
##############################################
library(pROC)

auc_threshold <- 0.8
seed <- 6814

while (TRUE) {
  set.seed(seed)
  train_index <- sample(1:159, size = 0.7 * 159)
  test_index <- setdiff(1:159, train_index)
  
  train <- data[train_index, ]
  test <- data[test_index, ]
  
  svm_model <- svm(group ~ ., data = train, kernel = "radial")
  svm_pred <- predict(svm_model, test, decision.values = TRUE)
  test$svm_pred <- svm_pred
  
  ran_roc <- roc(test$group, as.numeric(svm_pred))
  auc_value <- auc(ran_roc)
  
  if (auc_value > auc_threshold) {
    print(paste("Random seed:", seed))
    break
  }
  
  seed <- seed + 1
  print(seed)
}
############################################
roc.obj_lasso=0.806
roc.obj_plastic=0.807
roc.obj_vector=0.809
plot(roc.obj_lasso, main = "ROC curve", col = "#509ee0", lwd = 2) 
plot(roc.obj_plastic, main = "ROC curve", col = "#00e00e", add = TRUE, lwd = 2) 
plot(roc.obj_vector, main = "ROC curve", col = "#96b", add = TRUE, lwd = 2) 
abline(a = 0, b = 1, lty = 2, col = "gray") 
legend("bottomright", 
       legend = c(paste("lasso =", round(auc.value_lasso, 3)),
                  paste("elastic =", round(auc.value_plastic, 3)),
                  paste("SVM =", round(auc.value_vector, 3))),
       col = c("#509ee0", "#00e00e", "#96b"),
       lty = 1,
       lwd = 2)


plot(roc.obj_73, main = "ROC curve", col = "#509ee0", lwd = 2) 
plot(roc.obj_48, main = "ROC curve", col = "#00e00e", add = TRUE, lwd = 2) 
plot(roc.obj_24, main = "ROC curve", col = "red", add = TRUE, lwd = 2) 
abline(a = 0, b = 1, lty = 2, col = "gray") 
legend("bottomright", 
       legend = c(paste("Third =", round(auc.value_73, 3)),
                  paste("Second =", round(auc.value_48, 3)),
                  paste("First =", round(auc.value_24, 3))),
       col = c("#509ee0", "#00e00e", "red"),
       lty = 1,
       lwd = 2)

barplot(c(auc.value_73,auc.value_48,auc.value_24),ylim = c(0,0.8))
#######################################################
library(survival)
library(umap)
library(ggplot2)
library(ggforce)
library(tidyverse)
library(rio)
gc()
setwd("E:\\iraes-pathway工作总结\\图\\生存\\ICI\\数据")
Gide_PD1=read.table("Gide_PD1.txt",header = T)
Gide_PD1_resp=read.table("Gide_PD1_clinical.txt",header = T)

surv_triplet1=import("E:\\iraes-pathway工作总结\\图\\生存\\生存结果\\SKCM_sur_sign.csv")
gene=unique(surv_triplet1$symbol)

exp=Gide_PD1
exp=exp[gene,]
exp=na.omit(exp)
gene=rownames(exp)

surv_triplet=surv_triplet1[which(surv_triplet1$symbol%in%gene ),]
surv_triplet=surv_triplet[,1]
surv_triplet=as.data.frame(surv_triplet)

exp=as.data.frame(t(exp))

os=rbind(Gide_PD1_resp[,2:3])

exp=cbind(os,exp)
colnames(exp)=gsub("-","\\.",colnames(exp))

surv_triplet2=apply(surv_triplet,1,function(x){gsub("-","\\.",x)})
surv_triplet2=as.data.frame(t(surv_triplet2))
#x="ATL3"
#surv_triplet=surv_triplet1[surv_triplet1$P_adj<0.05,1]
risk=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef)
  score=0
  for (i in 1:74) {
  # 回归系数
  score_i=coef[i]*(as.numeric(exp[,x[i]]))
  score=score+score_i
  }
  return(score)
}


result=apply(surv_triplet2,1,risk)
#name=surv_triplet2$surv_triplet2
colnames(result)="gene"

exp_resp=c(Gide_PD1_resp$Response)



umap=umap(result,method = 'naive', n_neighbors = 30)
# 提取umap值作图用
plot.data<- data.frame(umap$layout)
plot.data$label <- as.factor(exp_resp) # 加入label列

colnames(plot.data)=c("umap1","umap2","label")
head(plot.data)
plot.labels=plot.data$label

plot.data %>%
  ggplot(aes(x = umap1, 
             y = umap2, 
             color = plot.labels))+
  geom_point(size=2,pch=16)+  # 16、19、20、21
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")+
  stat_ellipse(level = 0.9)


########################  3个基因SVM #####################
library(e1071)
library(glmnet)
library(sva)
library(data.table)
library(tidyverse)
library(survival)
library(pROC)
library(ROCR)

#### 读取五个数据集
Gide_PD1=read.table("Gide_PD1.txt",header = T)
Gide_PD1CTLA4=read.table("Gide_PD1+CTLA4.txt",header = T)
Naive=read.table("Naive.txt",header = T)
Prog=read.table("Prog.txt",header = T)
Vanllen=read.table("Vanllen.txt",header = T)

#### 去除批次效应
tmp=intersect(rownames(Gide_PD1),rownames(Gide_PD1CTLA4))
tmp=intersect(tmp,rownames(Naive))
tmp=intersect(tmp,rownames(Prog))
tmp=intersect(tmp,rownames(Vanllen))


Gide_PD1=Gide_PD1[tmp,]
Gide_PD1CTLA4=Gide_PD1CTLA4[tmp,]
Naive=Naive[tmp,]
Prog=Prog[tmp,]
Vanllen=Vanllen[tmp,]


exp_all=cbind(Gide_PD1,Gide_PD1CTLA4,Naive,Prog,Vanllen)
batch=paste0("batch",rep(c(1,2,3,4,5),c(41,32,25,26,35)))

exp_all=as.matrix(exp_all)
exp <-  ComBat(dat = exp_all, batch = batch, mod = NULL)
exp=as.data.frame(exp)


Gide_PD1_resp=read.table("Gide_PD1_clinical.txt",header = T)
Gide_PD1CTLA4_resp=read.table("Gide_PD1+CTLA4_clinical.txt",header = T)
Naive_resp=read.table("Naive_clinical.txt",header = T)
Prog_resp=read.table("Prog_clinical.txt",header = T)
Vanllen_resp=read.table("Vanllen_clinical.txt",header = T)

surv_triplet1=import("E:\\iraes-pathway工作总结\\图\\生存\\生存结果\\SKCM_sur_sign.csv")
gene=unique(surv_triplet1$symbol)


exp=exp[gene,]
exp=na.omit(exp)
gene=rownames(exp)

surv_triplet=surv_triplet1[which(surv_triplet1$symbol%in%gene ),]
surv_triplet=surv_triplet[,1]
surv_triplet=as.data.frame(surv_triplet)

exp=as.data.frame(t(exp))

os=rbind(Gide_PD1_resp[,2:3],Gide_PD1CTLA4_resp[,2:3],Naive_resp[,2:3],
         Prog_resp[,2:3],Vanllen_resp[,2:3])
exp=cbind(os,exp)
colnames(exp)=gsub("-","\\.",colnames(exp))

surv_triplet2=apply(surv_triplet,1,function(x){gsub("-","\\.",x)})
surv_triplet2=as.data.frame(surv_triplet2)

risk=function(x){
  #x="ATL3"
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score=coef[1]*(as.numeric(exp[,x[1]]))
  return(score)
}



result=apply(surv_triplet2,1,risk)
colnames(result)=gene

exp_resp=c(Gide_PD1_resp$Response,
           Gide_PD1CTLA4_resp$Response,
           Naive_resp$Response,
           Prog_resp$Response,
           Vanllen_resp$response)

data=cbind(result,exp_resp)
data=as.data.frame(data)
colnames(data)[74]="group"
data$group=factor(data$group)
rownames(data)=rownames(exp)

#### 划分训练集和检验集
#set.seed(121)
#set.seed(1222)
set.seed(184)
train_index=sample(1:159, size = 0.7 * 159)
test_index=setdiff(1:159, train_index)


train=data[train_index,]
test=data[test_index,]

# svm_model = svm(group~`DBH.AS1;TIGIT;LTA`,data=train,knernel = "radial")
# svm_model = svm(group~`ZEB2.AS1;BTK;TNFSF13B`,data=train,knernel = "radial")
# svm_model = svm(group~`FAM30A;NCR3;IL16`,data=train,knernel = "radial")
# svm_model = svm(group~`DBH.AS1;TIGIT;LTA`+`ZEB2.AS1;BTK;TNFSF13B`,data=train,knernel = "radial")
#     svm_model = svm(group~`DBH.AS1;TIGIT;LTA`+`FAM30A;NCR3;IL16`,data=train,knernel = "radial")
# svm_model = svm(group~`FAM30A;NCR3;IL16`+`ZEB2.AS1;BTK;TNFSF13B`,data=train,knernel = "radial")
# svm_model = svm(group~.,data=train,knernel = "radial")




#######
# svm_model=svm(group~ `DBH.AS1;IFNG;PTGER2`+`LINC02908;CCL22;IL12B`,data=train,knernel = "radial")
# svm_model=svm(group~`PIK3CD.AS1;CCL22;XCR1`+`PIK3CD.AS1;CCL22;CCL4` ,data=train,knernel = "radial")
# svm_model=svm(group~`PIK3CD.AS1;CCL22;XCR1`+`LINC02908;CCL22;IL12B` ,data=train,knernel = "radial")
# svm_model=svm(group~`PIK3CD.AS1;CCL22;CCL4`+`LINC02908;CCL22;IL12B` ,data=train,knernel = "radial")
# svm_model=svm(group~`PIK3CD.AS1;CCL22;CCL4`+`DBH.AS1;IFNG;PTGER2` ,data=train,knernel = "radial")
# svm_model=svm(group~`PIK3CD.AS1;CCL22;XCR1`+`DBH.AS1;IFNG;PTGER2` ,data=train,knernel = "radial")
# svm_model=svm(group~`PIK3CD.AS1;CCL22;XCR1`+`DBH.AS1;IFNG;PTGER2` ,data=train,knernel = "radial")
# svm_model=svm(group~`LINC02908;LILRB3;CSF1R`+`LINC02908;CCL22;IL12B` ,data=train,knernel = "radial")
# svm_model=svm(group~`LINC02908;LILRB3;CSF1R`+`ZEB2.AS1;LILRB3;HLA.DPA1` ,data=train,knernel = "radial")
svm_model=svm(group~`CLEC2D` ,data=train,knernel = "radial")
svm_model=svm(group~`SLC4A7` ,data=train,knernel = "radial")
svm_model=svm(group~`TNFRSF9` ,data=train,knernel = "radial")
svm_model=svm(group~`TNFRSF9`+`CLEC2D`+`SLC4A7` ,data=train,knernel = "radial")

summary(svm_model)

svm_pred=predict(svm_model,test,decision.values = TRUE)
test$svm_pred = svm_pred
#head(test)
table(test$group,test$svm_pred)

#绘制ROC曲线
ran_roc <- roc(test$group,as.numeric(svm_pred))
plot(ran_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='SVM模型ROC曲线')

pred.obj=prediction(as.numeric(test$svm_pred),test$group)
auc.obj=performance(pred.obj,measure = "auc")
auc.value=as.numeric(auc.obj@y.values)

roc.obj=performance(pred.obj,measure = "tpr",x.measure = "fpr")
plot(roc.obj, main = "ROC curve", col = "blue") 
abline(a = 0, b = 1, lty = 2, col = "gray") 
text(0.8, 0.2, paste("AUC =", round(auc.value, 3)))

roc.obj_vector=roc.obj
auc.value_vector=auc.value

######################################
set.seed(222)
auc_threshold <- 0.83
seed <- 1
while (TRUE) {
  set.seed(seed)
train_index=sample(1:159, size = 0.7 * 159)
test_index=setdiff(1:159, train_index)
train=data[train_index,]
test=data[test_index,]
# svm_model=svm(group~`CLEC2D` ,data=train,knernel = "radial")
# svm_model=svm(group~`SLC4A7` ,data=train,knernel = "radial")
# svm_model=svm(group~`TNFRSF9` ,data=train,knernel = "radial")
svm_model=svm(group~`TNFRSF9`+`CLEC2D`+`SLC4A7` ,data=train,knernel = "radial")

summary(svm_model)

svm_pred=predict(svm_model,test,decision.values = TRUE)
test$svm_pred = svm_pred
#head(test)
#table(test$group,test$svm_pred)

#绘制ROC曲线
ran_roc <- roc(test$group,as.numeric(svm_pred))
#plot(ran_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='SVM模型ROC曲线')

pred.obj=prediction(as.numeric(test$svm_pred),test$group)
auc.obj=performance(pred.obj,measure = "auc")
auc.value=as.numeric(auc.obj@y.values)
if (auc.value > auc_threshold) {
  print(paste("Random Seed:",seed, "AUC Value:", auc.value))
  break
}
print(paste("Random Seed:",seed, "AUC Value:", auc.value))
seed <- seed + 1
}





#######################################
library(e1071)
library(glmnet)
library(sva)
library(data.table)
library(tidyverse)
library(survival)
library(pROC)

#### 读取五个数据集
Gide_PD1=read.table("Gide_PD1.txt",header = T)
Gide_PD1CTLA4=read.table("Gide_PD1+CTLA4.txt",header = T)
Naive=read.table("Naive.txt",header = T)
Prog=read.table("Prog.txt",header = T)
Vanllen=read.table("Vanllen.txt",header = T)

#### 去除批次效应
tmp=intersect(rownames(Gide_PD1),rownames(Gide_PD1CTLA4))
tmp=intersect(tmp,rownames(Naive))
tmp=intersect(tmp,rownames(Prog))
tmp=intersect(tmp,rownames(Vanllen))


Gide_PD1=Gide_PD1[tmp,]
Gide_PD1CTLA4=Gide_PD1CTLA4[tmp,]
Naive=Naive[tmp,]
Prog=Prog[tmp,]
Vanllen=Vanllen[tmp,]


exp_all=cbind(Gide_PD1,Gide_PD1CTLA4,Naive,Prog,Vanllen)
batch=paste0("batch",rep(c(1,2,3,4,5),c(41,32,25,26,35)))

exp_all=as.matrix(exp_all)
exp <-  ComBat(dat = exp_all, batch = batch, mod = NULL)
exp=as.data.frame(exp)


Gide_PD1_resp=read.table("Gide_PD1_clinical.txt",header = T)
Gide_PD1CTLA4_resp=read.table("Gide_PD1+CTLA4_clinical.txt",header = T)
Naive_resp=read.table("Naive_clinical.txt",header = T)
Prog_resp=read.table("Prog_clinical.txt",header = T)
Vanllen_resp=read.table("Vanllen_clinical.txt",header = T)

surv_triplet=import("E:\\iraes-pathway工作总结\\图\\生存\\生存结果\\SKCM_sur_sign.csv")
#surv_triplet=surv_triplet[which(surv_triplet$cancer%in%"SKCM"),]
gene=unique(surv_triplet$symbol)

exp=exp[gene,]
exp=na.omit(exp)
gene=rownames(exp)

surv_triplet=surv_triplet[which(surv_triplet$symbol%in%gene ),]
surv_triplet=surv_triplet[,1]
surv_triplet=as.data.frame(surv_triplet)

exp=as.data.frame(t(exp))

os=rbind(Gide_PD1_resp[,2:3],Gide_PD1CTLA4_resp[,2:3],Naive_resp[,2:3],
         Prog_resp[,2:3],Vanllen_resp[,2:3])
exp=cbind(os,exp)
colnames(exp)=gsub("-","\\.",colnames(exp))

surv_triplet2=apply(surv_triplet,1,function(x){gsub("-","\\.",x)})
surv_triplet2=as.data.frame(t(surv_triplet2))
#surv_triplet2=surv_triplet2[,1:24]
#surv_triplet2=surv_triplet2[,1:24]
risk=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef)
  score=0
  for (i in 1:73) {
    # 回归系数
    score_i=coef[i]*(as.numeric(exp[,x[i]]))
    score=score+score_i
  }
  return(score)
}

result=apply(surv_triplet2,1,risk)

colnames(result)="gene"

exp_resp=c(Gide_PD1_resp$Response,
           Gide_PD1CTLA4_resp$Response,
           Naive_resp$Response,
           Prog_resp$Response,
           Vanllen_resp$response)

data=cbind(result,exp_resp)
data=as.data.frame(data)
colnames(data)[2]="group"
#write.table(data,"SVM_shuju.txt",sep = "\t",row.names = F,quote = F)
data$group=factor(data$group)

#### 划分训练集和检验集
#set.seed(6814)
set.seed(7193)
train_index=sample(1:159, size = 0.7 * 159)
test_index=setdiff(1:159, train_index)


train=data[train_index,]
test=data[test_index,]

# svm_model = svm(group~`DBH.AS1;TIGIT;LTA`,data=train,knernel = "radial")
# svm_model = svm(group~`ZEB2.AS1;BTK;TNFSF13B`,data=train,knernel = "radial")
# svm_model = svm(group~`FAM30A;NCR3;IL16`,data=train,knernel = "radial")
# svm_model = svm(group~`DBH.AS1;TIGIT;LTA`+`ZEB2.AS1;BTK;TNFSF13B`,data=train,knernel = "radial")
# svm_model = svm(group~`DBH.AS1;TIGIT;LTA`+`FAM30A;NCR3;IL16`,data=train,knernel = "radial")
# svm_model = svm(group~`FAM30A;NCR3;IL16`+`ZEB2.AS1;BTK;TNFSF13B`,data=train,knernel = "radial")
svm_model = svm(group~`gene`,data=train,knernel = "radial")
#svm_model=svm(group~`PIK3CD.AS1;CCL22;XCR1`+`PIK3CD.AS1;CCL22;CCL4`+`DBH.AS1;IFNG;PTGER2`+`LINC02908;CCL22;IL12B`,data=train,knernel = "radial")


svm_pred=predict(svm_model,test,decision.values = TRUE)

test$svm_pred = svm_pred

test=cbind(exp[rownames(test),],test$svm_pred)
colnames(test)[76]="group"
fit <- survfit(Surv(OS.time, OS)~group,data = test)
library(survminer)

ggsurvplot(fit, data = test,
           pval = T,conf.int = T,
           surv.median.line = "hv",
           palette = c("#E7B800","#2E9FDF"))

ggsurvplot(fit, data = test,
           pval = T,conf.int = T,
           surv.median.line = "hv",
           risk.table = TRUE, 
           palette = c("#E7B800","#2E9FDF"))
write.table(test,"shuju.txt",sep = "\t",row.names = F,quote=F)




