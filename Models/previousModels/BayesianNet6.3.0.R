#------------------------------------------------------------
#Cumulative model with seperate network at each timepoint
#learning with all time
#EM structure learning
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
  
  timesplitNumber = 10
  queryMethod = 'predict'
  
  timesplit = timeSplitFunction(numTimepoint = timesplitNumber, method = 'fixrate', data = training, debug=T)
  numTimepoint = length(timesplit)
  print(numTimepoint)
  
  allData = rbind(testing,training)
  
  cleaned = cleanDataType(nrow(testing),allData,discretize_level=3)
  testing = cleaned$test
  training = cleaned$train
  noCensorTraining = training[training$delta == 1,]
  noCensoredData = createTimeSplitIntegrated(noCensorTraining,timesplit)
  learningData = createTimeSplitIntegrated(training,timesplit)
  
  noCensoredData$time = NULL
  noCensoredData$delta = NULL
  learningData$time = NULL
  learningData$delta = NULL
  
  #blacklist2 = blacklistFunction(numTimepoint,colnames(data))
  blacklist = blacklistFunctionStart(numTimepoint)
  #blacklist = rbind(blacklist,blacklist2)
  
  #eval <- cnSearchSA(data = learningData)
  #bnet <- cnFindBIC(object = eval)
  
  #dag = bnlearn::empty.graph(names(learningData))
  #arcs(dag) = cnMatEdges(bnet)
  #dag = hc(learningData,blacklist=blacklist,score='bic')
  dag = structural.em(learningData, maximize = "hc",maximize.args = list(), fit = "mle",return.all = T,start = NULL, max.iter = 50, debug = FALSE)
  #dag = cextend(dag, strict = TRUE, debug = FALSE)
  #fit = bn.fit(dag, learningData, method='mle')
  fit = dag$fitted
  plotDag(dag$dag)
  
  
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
        str2 = paste0(paste0("t_",i),'==0')
        prob = cpquery(fit, event = eval(parse(text = str2)), evidence = evidence,method = 'lw')
      }else if(method == 'predict') {
        evidence = testing[j,]
        predicted = predict(fit, node=paste0("t_",i), evidence,method = "bayes-lw", prob = TRUE, n=500)
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
  if(numReturnNA>0) {cat('return NA: ',numReturnNA)}
  if(numTooLow>0) {cat('too Low: ',numTooLow)}
  if(numNotDecreasing>0) {cat('Not decreasing: ',numNotDecreasing)}
  colnames(survivalFunction) = 1:nrow(testing)
  survivalFunction = cbind(time = timesplit, survivalFunction) 
  return(survivalFunction)
}

plotDag <-function(dag){
  #plot with bigger fontsize
  g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(dag))
  graph::nodeRenderInfo(g) <- list(fontsize=100)
  Rgraphviz::renderGraph(g)
  
  #graphviz.compare(dagList[[1]],start, shape = "rectangle")
  #graphviz.compare(startDagList[[10]],startDagList[[11]], shape = "rectangle",diff.args = list(tp.lwd = 2, fp.col = "orange"))
}
blacklistFunction <- function(numTimepoint,nodeNames){
  numCol = numTimepoint*length(nodeNames)
  blackList = matrix(ncol = 2, nrow = numCol)
  index = 1
  for(i in 1:numTimepoint) {
    for(j in 1:length(nodeNames)) {
      blackList[index,] = c(paste0("t_",i),nodeNames[j])
      index = index+1
    }
  }
  for(i in 1:numTimepoint) {
    blackList = blackList[blackList[,2] != paste0("t_",i),]
  }
  
  return(blackList)
}

blacklistFunctionStart <- function(numTimepoint){
  numCol = (numTimepoint-2)*numTimepoint+1
  blackList = matrix(ncol = 2, nrow = numCol)
  index = 1
  for(i in 1:numTimepoint) {
    for(j in 1:numTimepoint) {
      if(i!=j & j!=i+1) {
        blackList[index,] = c(paste0("t_",i),paste0("t_",j))
        index = index+1
      }
    }
  }
  
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