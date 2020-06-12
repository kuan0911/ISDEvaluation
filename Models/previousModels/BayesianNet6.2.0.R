#------------------------------------------------------------
#Cumulative model with seperate network at each timepoint
#learning with all time
#
#
#------------------------------------------------------------

library(bnlearn)
library(plyr)
library(gtools)
library(e1071)
source("Models/cleanDataForBayesNet.R")

BayesianNet = function(training, testing){
  originalTraining = training
  originalTesting = testing
  
  timesplitNumber = 40
  queryMethod = 'predict'
  prior = F
  
  timesplit = timeSplitFunction(numTimepoint = timesplitNumber, method = 'quantile', data = training, debug=T)
  numTimepoint = length(timesplit)
  print(numTimepoint)
  
  allData = rbind(testing,training)
  
  cleaned = cleanDataType(nrow(testing),allData,discretize_level=5)
  testing = cleaned$test
  training = cleaned$train
  noCensorTraining = training[training$delta == 1,]
  noCensoredData = createTimeSplitIntegrated(noCensorTraining,timesplit)
  data = createTimeSplitIntegrated(training,timesplit)
  
  noCensoredData$time = NULL
  noCensoredData$delta = NULL
  data$time = NULL
  data$delta = NULL
  dag = hc(data,score='bic')
  dag = cextend(dag, strict = TRUE, debug = FALSE)
  fit = bn.fit(dag, data, method='mle')
  #plotDag(dag)
  
  
  print('start predict')
  #prediction
  survivalFunctionTesting = predictFunction(fit,testing,timesplit,method = queryMethod)
  survivalFunctionTraining = predictFunction(fit,training,timesplit,method = queryMethod)
  
  testCurvesToReturn = survivalFunctionTesting
  timesAndCensTest = cbind.data.frame(time = originalTesting$time, delta = originalTesting$delta)
  timesAndCensTrain = cbind.data.frame(time = originalTraining$time, delta = originalTraining$delta)
  trainingCurvesToReturn = survivalFunctionTraining
  
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}

predictFunction <- function(fit,testing,timesplit,method) {
  testing$time = NULL
  testing$delta = NULL
  numTimepoint = length(timesplit)
  numReturnNA = 0
  numTooLow = 0
  numNotDecreasing = 0
  previousTimepointProb = rep(1,nrow(testing))
  
  #first value in cpt
  survivalFunction <- data.frame(matrix(ncol = nrow(testing), nrow = numTimepoint))
  for(i in 1:numTimepoint) {
    print(i)
    for(j in 1:nrow(testing)) {
      
      if(method == 'lw') {
        evidence = as.list(testing[j,])
        prob = cpquery(fit, event = (TIMEPOINT == 0), evidence = evidence,method = 'lw')
      }else if(method == 'predict') {
        evidence = testing[j,]
        predicted = predict(fit, node=paste0("t_",i), evidence,method = "bayes-lw", prob = TRUE, n=2000)
        prob = attr(predicted, "prob")[1]
      }else{
        print('No such query method')
        break
      }
      
      survivalFunction[i,j] = prob
      
      if(is.na(survivalFunction[i,j])) {
        survivalFunction[i,j] = previousTimepointProb[j]-0.01
        if(survivalFunction[i,j] <0){survivalFunction[i,j] = 0}
        numReturnNA = numReturnNA+1
      }
      if(i==1 & survivalFunction[i,j]<0.6){
        survivalFunction[i,j] = 0.6
        numTooLow = numTooLow+1
      }
      if(previousTimepointProb[j]<=survivalFunction[i,j]){
        survivalFunction[i,j] = previousTimepointProb[j]-0.00001
        if(survivalFunction[i,j] <0){survivalFunction[i,j] = 0}
        numNotDecreasing = numNotDecreasing+1
        #print('probability not decreasing')
      }
      previousTimepointProb[j] = survivalFunction[i,j]
      
    }
  }
  if(numReturnNA>0) {cat('return NA: numReturnNA',numReturnNA)}
  if(numTooLow>0) {cat('too Low: ',numTooLow)}
  if(numNotDecreasing>0) {cat('Not decreasing: ',numNotDecreasing)}
  colnames(survivalFunction) = 1:nrow(testing)
  survivalFunction = cbind(time = timesplit, survivalFunction) 
  return(survivalFunction)
}

plotDag <-function(dag){
  #plot with bigger fontsize
  g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(dag))
  graph::nodeRenderInfo(g) <- list(fontsize=80)
  Rgraphviz::renderGraph(g)
  
  #graphviz.compare(dagList[[1]],start, shape = "rectangle")
  #graphviz.compare(startDagList[[10]],startDagList[[11]], shape = "rectangle",diff.args = list(tp.lwd = 2, fp.col = "orange"))
}
blacklistFunction <- function(nodeNames){
  blackList = matrix(ncol = 2, nrow = length(nodeNames))
  blackList2 = matrix(ncol = 2, nrow = length(nodeNames))
  blackList3 = matrix(ncol = 2, nrow = length(nodeNames))
  blackList4 = matrix(ncol = 2, nrow = length(nodeNames))
  for(i in 1:length(nodeNames)) {
    blackList[i,] = c(nodeNames[i],"PREVTIMEPOINT")
    blackList2[i,] = c("TIMEPOINT", nodeNames[i])
    blackList3[i,] = c("PREVTIMEPOINT", nodeNames[i])
    blackList4[i,] = c(nodeNames[i],"TIMEPOINT")
  }
  blackList = blackList[blackList[,1] != 'PREVTIMEPOINT',]
  blackList2 = blackList2[blackList2[,2] != 'TIMEPOINT',]
  blackList2 = blackList2[blackList2[,2] != 'PREVTIMEPOINT',]
  blackList3 = blackList3[blackList3[,2] != 'TIMEPOINT',]
  blackList3 = blackList3[blackList3[,2] != 'PREVTIMEPOINT',]
  blackList4 = blackList4[blackList4[,1] != 'TIMEPOINT',]
  blackList4 = blackList4[blackList4[,1] != 'PREVTIMEPOINT',]
  
  #blackList = rbind(blackList, blackList2)
  #blackList = rbind(blackList, blackList3)
  
  return(list(blackList1=blackList2,blackList2=blackList4))
}

blacklistFunctionStart <- function(nodeNames){
  blackList = matrix(ncol = 2, nrow = length(nodeNames))
  blackList2 = matrix(ncol = 2, nrow = length(nodeNames))
  for(i in 1:length(nodeNames)) {
    blackList[i,] = c(nodeNames[i],"TIMEPOINT")
    blackList2[i,] = c("TIMEPOINT", nodeNames[i])
  }
  blackList = blackList[blackList[,1] != 'TIMEPOINT',]
  blackList2 = blackList2[blackList2[,2] != 'TIMEPOINT',]
  
  blackList = rbind(blackList, blackList2)
  
  return(blackList)
}

commaStringFunction <- function(p){
  commaString1 = replace(p,p!="TIMEPOINT"&p!="PREVTIMEPOINT",'')
  commaString1 = replace(commaString1,commaString1=="TIMEPOINT",'1')
  #commaString1 = paste0(replace(commaString1,commaString1=="PREVTIMEPOINT",'1'),collapse = ',')
  commaString1 = paste0(commaString1,collapse = ',')
  commaString2 = replace(p,p!="TIMEPOINT"&p!="PREVTIMEPOINT",'')
  commaString2 = replace(commaString2,commaString2=="TIMEPOINT",'2')
  #commaString2 = paste0(replace(commaString2,commaString2=="PREVTIMEPOINT",'1'),collapse = ',')
  commaString2 = paste0(commaString2,collapse = ',')
  commaString3 = replace(p,p!="TIMEPOINT"&p!="PREVTIMEPOINT",'')
  commaString3 = replace(commaString3,commaString3=="TIMEPOINT",'1')
  commaString3 = paste0(replace(commaString3,commaString3=="PREVTIMEPOINT",'2'),collapse = ',')
  commaString4 = replace(p,p!="TIMEPOINT"&p!="PREVTIMEPOINT",'')
  commaString4 = replace(commaString4,commaString4=="TIMEPOINT",'2')
  commaString4 = paste0(replace(commaString4,commaString4=="PREVTIMEPOINT",'2'),collapse = ',')
  return(list(commaString1,commaString2,commaString3,commaString4))
}

independent_test <-function(data){
  
  independence_test = data.frame(matrix(ncol = ncol(data), nrow = 3))
  names(independence_test) = names(data)
  independence_test$TIMEPOINT = NULL
  row.names(independence_test) = c('mi','df','p.value')
  for(c in names(data)) {
    if(c != 'TIMEPOINT') {
      res = ci.test(x="TIMEPOINT", y=c, data = data, test = "mi")
      independence_test[[c]] = c(res$statistic, res$parameter, res$p.value)
    }
  }
  
  topname = names(independence_test)[independence_test['mi',]>100]
  if(length(topname)==0) {
    topname = names(independence_test)[head(order(independence_test['mi',],decreasing=TRUE),3)]
  }
  return(topname)
}