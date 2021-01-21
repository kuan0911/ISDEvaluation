#------------------------------------------------------------
#Cumulative model with seperate network at each timepoint
#include all data in structure learning. Na if dead or censored
#
#
#------------------------------------------------------------

library(bnlearn)
library(plyr)
library(gtools)
#library(e1071)
library(prodlim)
library(MTLR)
library(Rcpp)
library(glmnet)
source("Models/cleanDataForBayesNet.R")
source("Models/bayesianNetHelper.R")
source("Models/fractionallyInclude.R")
source("Models/changesource.R")
source("ValidateCleanCV/createFoldsAndNormalize.R")
library(keras)
library(tensorflow)
library(magrittr)
library(dplyr)


DHBNglm = function(training, testing, timePoints = 0,debug = FALSE){
  originalTraining = training
  originalTesting = testing
  
  # C1 = NULL
  # if(is.null(C1)){
  #   C1 = mtlr_cv(Surv(time,delta)~.,data=originalTraining, loss= "conc", C1_vec = c(0.001,0.01,0.1,1,10,100,1000)
  #                , train_biases = F,train_uncensored = F)$best_C1
  # }
  #mod = mtlr(Surv(time,delta)~., data = originalTraining, C1=C1, train_biases = F, train_uncensored = F)
  mod = NULL
  kmMod = prodlim(Surv(time,delta)~1, data = originalTraining)
  
  queryMethod = 'exact'
  prior = F
  weighted = F

  variableList = variableTypes(originalTraining,10)
  
  print(timePoints)
  
  numTimepoint = length(timePoints)
  
  dataListNA = createTimeSplit(training,timePoints,includeNa = T)

  print(length(dataListNA))
  
  fitList <- vector("list", numTimepoint)
  print('learn starting graph')
  
  for(iter in 1:1) {
    curveCache = matrix(NA,nrow(dataListNA[[1]]),length(timePoints))
    if(iter>1) {
      for(k in 1:nrow(dataListNA[[1]])) {
        if(is.na(dataListNA[[1]][k,'TIMEPOINT'])) {
          evidence = dataListNA[[1]][k,]
          evidence[c('PREVTIMEPOINT','TIMEPOINT','time','delta','id')] = NULL
          curveCache[k,] = predictLRCurve(fitList,evidence,length(timePoints))
        }
      }
    }
    cbind(dataListNA[[1]]$id,curveCache)
    
    for(i in 1:numTimepoint) {
      if(isTRUE(debug)){cat(i)
        cat(' ')}
      
      data = dataListNA[[i]]
      
      dataHazardNA = data
      dataHazardNA = dataHazardNA[is.na(dataHazardNA$PREVTIMEPOINT)|dataHazardNA$PREVTIMEPOINT==0,]
      resKM = weightedImputeKM(kmMod,dataHazardNA,timePoints,i)
      weightsNAKM = resKM$weight
      dataHazardNAKM = resKM$data
      resLR = weightedImputeLR(fitList,dataHazardNA,timePoints,i,curveCache)
      weightsNALR = resLR$weight
      dataHazardNALR = resLR$data
      
      dataHazard = data
      dataHazard = dataHazard[!is.na(dataHazard$PREVTIMEPOINT),]
      dataHazard = dataHazard[dataHazard$PREVTIMEPOINT == 0,]
      weights = rep(1,nrow(dataHazard))
      weights[is.na(dataHazard$TIMEPOINT)] = 0.5
      dataHazard[is.na(dataHazard$TIMEPOINT),'TIMEPOINT'] = 0
      
      if(iter==1) {
        trainingData = dataHazard
        weights = weights
      }else{
        trainingData = dataHazardNALR
        weights = weightsNALR
      }
      
      
      cvFoldIndex = createFoldsOfData(trainingData, numberOfFolds=5)[[1]]
      transCvFoldIndex = rep(1,nrow(trainingData))
      for(cvFoldIter in 1:length(cvFoldIndex)){transCvFoldIndex[cvFoldIndex[[cvFoldIter]]]=cvFoldIter}
      
      trainingData[,c('time','delta','id','PREVTIMEPOINT')] = NULL
      y = trainingData$TIMEPOINT
      trainingData$TIMEPOINT = NULL
      x = as.matrix(trainingData)
      
      if(sum(y == 0)<4 | sum(y == 1)<4) {
        fit = NULL
        cat('skip ')
      }else {
        fit = cv.glmnet(x,y,family='binomial', weights=weights,foldid=transCvFoldIndex)
        #fit <- glmnet(x,y,family='binomial',weight=weight)
      }
      fitList[[i]] <- fit
    }
  }
  
  # aveCoef = rep(0,length(fit$glmnet.fit$beta[,1]))
  # fitcount = 0
  # for(i in 1:length(fitList)) {
  #   fit = fitList[[i]]
  #   if(!is.null(fit)) {
  #     wh = which(fit$lambda==fit$lambda.min)
  #     aveCoef = aveCoef + fit$glmnet.fit$beta[,wh]
  #     fitcount = fitcount + 1
  #   }
  # }
  # aveCoef = aveCoef/fitcount
  # 
  # for(i in 1:length(fitList)) {
  #   fit = fitList[[i]]
  #   if(!is.null(fit)) {
  #     wh = which(fit$lambda==fit$lambda.min)
  #     delta = fit$glmnet.fit$beta[,wh] - aveCoef
  #     fit$glmnet.fit$beta[,wh] = fit$glmnet.fit$beta[,wh] - delta/5
  #   }
  #   #fitList[[i]] = fit
  # }
  # tail2 = length(fitList)-1
  # for(i in 2:tail2) {
  #   localCoef = rep(0,length(fit$glmnet.fit$beta[,1]))
  #   fitcount = 0
  #   if(!is.null(fitList[[i-1]])) {
  #     fit = fitList[[i-1]]
  #     wh = which(fit$lambda==fit$lambda.min)
  #     localCoef = localCoef+fit$glmnet.fit$beta[,wh]
  #     fitcount = fitcount + 1
  #   }
  #   if(!is.null(fitList[[i+1]])) {
  #     fit = fitList[[i+1]]
  #     wh = which(fit$lambda==fit$lambda.min)
  #     localCoef = localCoef+fit$glmnet.fit$beta[,wh]
  #     fitcount = fitcount + 1
  #   }
  #   if(!is.null(fitList[[i]])) {
  #     fit = fitList[[i]]
  #     wh = which(fit$lambda==fit$lambda.min)
  #     localCoef = localCoef+fit$glmnet.fit$beta[,wh]
  #     fitcount = fitcount + 1
  #     delta = fit$glmnet.fit$beta[,wh] - localCoef
  #     fit$glmnet.fit$beta[,wh] = fit$glmnet.fit$beta[,wh] - delta/5
  #   }
  #   fitList[[i]] = fit
  # }

  print('start predict')
  #prediction
  survivalFunctionTesting = predictFunctionLR(fitList,originalTesting,timePoints)
  survivalFunctionTraining = predictFunctionLR(fitList,training,timePoints)
  
  testCurvesToReturn = survivalFunctionTesting
  timesAndCensTest = cbind.data.frame(time = originalTesting$time, delta = originalTesting$delta)
  timesAndCensTrain = cbind.data.frame(time = originalTraining$time, delta = originalTraining$delta)
  trainingCurvesToReturn = survivalFunctionTraining
  
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}

predictFunctionLR <- function(fitList,testing,timePoints) {
  numTimepoint = length(timePoints)
  numReturnNA = 0
  numNotDecreasing = 0
  testing[,c('time','delta')] = NULL
  
  previousTimepointProb = rep(1,nrow(testing))
  
  survivalFunction <- data.frame(matrix(ncol = nrow(testing), nrow = numTimepoint))
  for(i in 1:numTimepoint) {
    if(i>length(fitList)){
      survivalFunction[i,] = previousTimepointProb
      break
    }
    fitted = fitList[[i]]
    if(is.null(fitted)) {
      prob = 1
    }else {
      testing = as.matrix(testing)
      prob = 1 - predict(fitted,newx=testing,type='response',s="lambda.min")
    }

    survivalFunction[i,] = prob*previousTimepointProb
    
    previousTimepointProb = survivalFunction[i,]
  }
  if(numReturnNA>0) {cat('return NA: ',numReturnNA)}
  if(numNotDecreasing>0) {cat('Not decreasing: ',numNotDecreasing)}
  
  colnames(survivalFunction) = 1:nrow(testing)
  return(survivalFunction)
}
