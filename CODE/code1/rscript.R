install.packages("foreach")
install.packages("doParallel")
library(foreach)
library(doParallel)

pathwayIterationFunction2 <- function(y, exp1_disease){
  
  all <- pathwayIterationResult[y,]
  all <- all[-1]
  all <- as.numeric(all[-length(all)])
  enrichmentScore <- pathwayIterationResult[y,ncol(pathwayIterationResult)]
  ###计算每个基因的差异得分Diff
  diffScoreFunction <- function(j, exp1_disease){
    geneName <- rownames(exp1_disease)[j]
    geneValue <- as.numeric(exp1_disease[j,])
    group <- ifelse(geneValue > median(geneValue), "high", "low")
    #if(sum(group=="low")==ncol(exp1_disease)){next}
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
    times = 1000 #随机1000次
    permuteResult <- foreach(i = seq(times),
                             .combine = rbind)  %dopar%  permutationFunction()
    d <- length(which(permuteResult==1)) 
    p <- d/1000
    if(enrichmentScore){
      risk_p <- t(c(geneName, unique(gmt$ont)[y], p)) #基因名+路径名+p值
    }else{
      risk_p <- t(c(geneName, unique(gmt$ont)[y], NA)) #基因名+路径名+p值(NA)
    }
    ######
    
    return(risk_p)
  }
  
  diffScoreResult <- foreach(j = 1:length(rownames(exp1_disease)),
                             .combine = rbind,
                             .packages = c('foreach', 'doParallel'))  %dopar%  diffScoreFunction(j, exp1_disease)
  return(diffScoreResult)
}

library(foreach)
library(doParallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

pathwayIterationResult2 <- foreach(y = 1:17,
                                   .combine = rbind,
                                   .packages = c('foreach', 'doParallel'))  %dopar%  pathwayIterationFunction2(y, exp1_disease)
stopCluster(cl)
save(pathwayIterationResult2, file = "./new.RData")










start <- proc.time()
pathwayIterationFunction2 <- function(y, exp1_disease){
  
  all <- pathwayIterationResult[y,]
  all <- all[-1]
  all <- as.numeric(all[-length(all)])
  enrichmentScore <- pathwayIterationResult[y,ncol(pathwayIterationResult)]
  if(length(which(is.na(all)==TRUE))!=length(all)){
    ###计算每个基因的差异得分Diff
    diffScoreFunction <- function(j, exp1_disease){
      geneName <- rownames(exp1_disease)[j]
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
      times = 1000 #随机1000次
      permuteResult <- foreach(i = seq(times),
                               .combine = rbind)  %dopar%  permutationFunction()
      d <- length(which(permuteResult==1)) 
      p <- d/1000
      if(enrichmentScore){
        risk_p <- t(c(geneName, unique(gmt$ont)[y], p)) #基因名+路径名+p值
      }else{
        risk_p <- t(c(geneName, unique(gmt$ont)[y], NA)) #基因名+路径名+p值(NA)
      }
      ######
      
      return(risk_p)
    }
    
    diffScoreResult <- foreach(j = 1:5,
                               .combine = rbind,
                               .packages = c('foreach', 'doParallel'))  %dopar%  diffScoreFunction(j, exp1_disease)
    return(diffScoreResult)
  }
 
}

library(foreach)
library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

pathwayIterationResult2 <- foreach(y = 1:17,
                                   .combine = rbind,
                                   .packages = c('foreach', 'doParallel'))  %dopar%  pathwayIterationFunction2(y, exp1_disease)
stopCluster(cl)
end <- proc.time()
running1 <- end-start 
save(pathwayIterationResult2, file = "./new.RData")








start <- Sys.time()
pathwayIterationFunction2 <- function(y, exp1_disease){
  
  all <- pathwayIterationResult[y,]
  all <- all[-1]
  all <- as.numeric(all[-length(all)])
  enrichmentScore <- pathwayIterationResult[y,ncol(pathwayIterationResult)]
  risk_p <- NULL
  if(length(which(is.na(all)==TRUE))!=length(all)){
    ###计算每个基因的差异得分Diff
    for(j in 1:5){
      geneName <- rownames(exp1_disease)[j]
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
      times = 1000 #随机1000次
      permuteResult <- foreach(i = seq(times),
                               .combine = rbind)  %dopar%  permutationFunction()
      d <- length(which(permuteResult==1)) 
      p <- d/1000
      if(enrichmentScore){
        risk_p <- rbind(risk_p, t(c(geneName, unique(gmt$ont)[y], p))) #基因名+路径名+p值
      }else{
        risk_p <- rbind(risk_p, t(c(geneName, unique(gmt$ont)[y], NA))) #基因名+路径名+p值(NA)
      }
      ######
      
    }
    
    diffScoreResult <- risk_p
    return(diffScoreResult)
  }
  
}

library(foreach)
library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

pathwayIterationResult2 <- foreach(y = 1:17,
                                   .combine = rbind,
                                   .packages = c('foreach', 'doParallel'))  %dopar%  pathwayIterationFunction2(y, exp1_disease)
stopCluster(cl)
end <- Sys.time()
running2 <- end-start 
save(pathwayIterationResult2, file = "./new.RData")