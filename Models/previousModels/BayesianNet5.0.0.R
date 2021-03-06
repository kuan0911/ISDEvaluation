#------------------------------------------------------------
#Cumulative model with seperate network at each timepoint
#add number of timepoint
#smooth the line
#hierarchical
#------------------------------------------------------------

library(bnlearn)
library(plyr)
library(gtools)
library(e1071)
source("Models/cleanDataForBayesNet.R")

BayesianNet = function(training, testing){
  originalTraining = training
  originalTesting = testing
  
  numTimepoint = 10
  queryMethod = 'naivebayes'
  prior = F
  
  timesplit = timeSplitFunction(numTimepoint = numTimepoint, method = 'fixrate', data = training, debug=T)
  allData = rbind(testing,training)
  
  numTimepoint2 = 5
  timesplit2 = timeSplitFunction(numTimepoint = numTimepoint2, method = 'fixrate', data = training, debug=T)
  
  numTimepoint3 = 2
  timesplit3 = timeSplitFunction(numTimepoint = numTimepoint3, method = 'fixrate', data = training, debug=T)
  
  numTimepoint4 = 4
  timesplit4 = timeSplitFunction(numTimepoint = numTimepoint4, method = 'quantile', data = training, debug=T)
  
  cleaned = cleanDataType(nrow(testing),allData,discretize_level=7)
  testing = cleaned$test
  training = cleaned$train
  #noCensorTraining = training[training$delta == 1,]
  dataList = createTimeSplit(training,timesplit)
  dataList2 = createTimeSplit(training,timesplit2)
  dataList3 = createTimeSplit(training,timesplit3)
  dataList4 = createTimeSplit(training,timesplit4)
  
  #learning structure and parameter
  fitList <- vector("list", numTimepoint)
  dagList <- vector("list", numTimepoint)
  nbList <- vector("list", numTimepoint)
  childrenList <- vector("list", numTimepoint)
  
  fitList2 <- vector("list", numTimepoint2)
  dagList2 <- vector("list", numTimepoint2)
  nbList2 <- vector("list", numTimepoint2)
  childrenList2 <- vector("list", numTimepoint2)
  
  fitList3 <- vector("list", numTimepoint3)
  dagList3 <- vector("list", numTimepoint3)
  nbList3 <- vector("list", numTimepoint3)
  childrenList3 <- vector("list", numTimepoint3)
  
  fitList4 <- vector("list", numTimepoint4)
  dagList4 <- vector("list", numTimepoint4)
  nbList4 <- vector("list", numTimepoint4)
  childrenList4 <- vector("list", numTimepoint4)
  
  #blackList = blacklistFunction(colnames(dataList[[1]]))
  
  print('learn starting graph')
  data = training[training$delta==1,]
  discretize_level = 100
  breaks = seq(0,100,100/discretize_level)
  breakList = c(-Inf, unique(quantile(data$time, probs = breaks[2:discretize_level]/100)), Inf)
  data$time = cut(data$time,breaks = breakList,labels = 1:(length(breakList)-1))
  data$delta = NULL
  data$TIMEPOINT = data$time
  data$time = NULL
  start = hc(data,score='bic',k=2)
  
  #whiteList = c('PREVTIMEPOINT','TIMEPOINT')
  for(i in 1:length(dataList)) {
    data = dataList[[i]]
    data = data[data[,'PREVTIMEPOINT']==0,]
    data$PREVTIMEPOINT = NULL
    
    dag = hc(data,score='bic',k=1)
    dag = cextend(dag, strict = TRUE, debug = FALSE)
    dagList[[i]] = dag 
    fit = bn.fit(dag, data, method='bayes')
    fitList[[i]] <- fit
    nb_model = naiveBayes(as.factor(TIMEPOINT) ~., data=data)
    nbList[[i]] <- nb_model
    childrenList[[i]] <- children(fitList[[i]],'TIMEPOINT')
    #print(mb(fitList[[i]],'TIMEPOINT'))
  }
  
  for(i in 1:length(dataList2)) {
    data = dataList2[[i]]
    data = data[data[,'PREVTIMEPOINT']==0,]
    data$PREVTIMEPOINT = NULL
    
    dag = hc(data,score='bic',k=1)
    dag = cextend(dag, strict = TRUE, debug = FALSE)
    dagList2[[i]] = dag 
    fit = bn.fit(dag, data, method='bayes')
    fitList2[[i]] <- fit
    nb_model = naiveBayes(as.factor(TIMEPOINT) ~., data=data)
    nbList2[[i]] <- nb_model
    childrenList2[[i]] <- children(fitList2[[i]],'TIMEPOINT')
    #print(mb(fitList[[i]],'TIMEPOINT'))
  }
  
  for(i in 1:length(dataList3)) {
    data = dataList3[[i]]
    data = data[data[,'PREVTIMEPOINT']==0,]
    data$PREVTIMEPOINT = NULL
    
    dag = hc(data,score='bic',k=1)
    dag = cextend(dag, strict = TRUE, debug = FALSE)
    dagList3[[i]] = dag 
    fit = bn.fit(dag, data, method='bayes')
    fitList3[[i]] <- fit
    nb_model = naiveBayes(as.factor(TIMEPOINT) ~., data=data)
    nbList3[[i]] <- nb_model
    childrenList3[[i]] <- children(fitList3[[i]],'TIMEPOINT')
    #print(mb(fitList[[i]],'TIMEPOINT'))
  }
  
  for(i in 1:length(dataList4)) {
    data = dataList4[[i]]
    data = data[data[,'PREVTIMEPOINT']==0,]
    data$PREVTIMEPOINT = NULL
    
    dag = hc(data,score='bic',k=1)
    dag = cextend(dag, strict = TRUE, debug = FALSE)
    dagList4[[i]] = dag 
    fit = bn.fit(dag, data, method='bayes')
    fitList4[[i]] <- fit
    nb_model = naiveBayes(as.factor(TIMEPOINT) ~., data=data)
    nbList4[[i]] <- nb_model
    childrenList4[[i]] <- children(fitList4[[i]],'TIMEPOINT')
    #print(mb(fitList[[i]],'TIMEPOINT'))
  }
  
  if(prior) {
    print('Impose prior to TIMEPOINT')
    for(i in 1:length(dataList)) {
      data = dataList[[i]]
      data = data[data[,'PREVTIMEPOINT']==0,]
      data$PREVTIMEPOINT = NULL
      cpt = coef(fitList[[i]]$TIMEPOINT)
      p = colnames(as.data.frame(cpt))
      p = p[p!="Freq"]
      if(p != "Var1"){
        newcpt = table(data[,p])
        commaString = commaStringFunction(p)
        codeStringTotal1 = paste0('newcpt[',commaString[1],']',collapse = '')
        codeStringTotal2 = paste0('newcpt[',commaString[2],']',collapse = '')
        #total = eval(parse(text=codeStringTotal1))+eval(parse(text=codeStringTotal2))
        timepoint1 = eval(parse(text=codeStringTotal1))
        timepoint2 = eval(parse(text=codeStringTotal2))
        numData = nrow(data)
        weight = 5
        for(j in 1:length(timepoint1)) {
          timepoint1[j] = rbeta(1,(numTimepoint-1)*weight*5+timepoint1[j]*5,1*weight*5+timepoint2[j]*5)
          #while(timepoint1[j]<0.5){
          #print(timepoint1[j])
          #  timepoint1[j] = rbeta(1,9+timepoint1[j]*5,1+timepoint2[j]*5)
          #}
        }
        
        codeString1 = paste0('newcpt[',commaString[1],']=','timepoint1',collapse = '')
        codeString2 = paste0('newcpt[',commaString[2],']=1-','newcpt[',commaString[1],']',collapse = '')
        eval(parse(text=codeString1))
        eval(parse(text=codeString2))
        fitList[[i]]$TIMEPOINT = newcpt
      }
      
    }
  }
  
  
  #plot with bigger fontsize
  #plotDag(dagList[[1]])
  print('start predict')
  #prediction
  survivalFunctionTesting = predictFunction(fitList,nbList,childrenList,testing,timesplit,method = queryMethod)
  survivalFunctionTraining = predictFunction(fitList,nbList,childrenList,training,timesplit,method = queryMethod)
  
  #survivalFunctionTesting = smoothCurve(survivalFunctionTesting)
  #survivalFunctionTraining = smoothCurve(survivalFunctionTraining)
  
  survivalFunctionTesting2 = predictFunction(fitList2,nbList2,childrenList2,testing,timesplit2,method = queryMethod)
  survivalFunctionTraining2 = predictFunction(fitList2,nbList2,childrenList2,training,timesplit2,method = queryMethod)
  
  #survivalFunctionTesting2 = smoothCurve(survivalFunctionTesting2)
  #survivalFunctionTraining2 = smoothCurve(survivalFunctionTraining2)
  
  survivalFunctionTesting3 = predictFunction(fitList3,nbList3,childrenList3,testing,timesplit3,method = queryMethod)
  survivalFunctionTraining3 = predictFunction(fitList3,nbList3,childrenList3,training,timesplit3,method = queryMethod)
  
  #survivalFunctionTesting3 = smoothCurve(survivalFunctionTesting3)
  #survivalFunctionTraining3 = smoothCurve(survivalFunctionTraining3)
  
  survivalFunctionTesting4 = predictFunction(fitList4,nbList4,childrenList4,testing,timesplit4,method = queryMethod)
  survivalFunctionTraining4 = predictFunction(fitList4,nbList4,childrenList4,training,timesplit4,method = queryMethod)
  
  #survivalFunctionTesting4 = smoothCurve(survivalFunctionTesting4)
  #survivalFunctionTraining4 = smoothCurve(survivalFunctionTraining4)
  
  testCurvesToReturn = smoothCurve(overlapCurve(survivalFunctionTesting,survivalFunctionTesting2,survivalFunctionTesting3,survivalFunctionTesting4))
  timesAndCensTest = cbind.data.frame(time = originalTesting$time, delta = originalTesting$delta)
  timesAndCensTrain = cbind.data.frame(time = originalTraining$time, delta = originalTraining$delta)
  trainingCurvesToReturn = smoothCurve(overlapCurve(survivalFunctionTraining,survivalFunctionTraining2,survivalFunctionTraining3,survivalFunctionTraining4))
  
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}

predictFunction <- function(fitList,nbList,childrenList,testing,timesplit,method) {
  testing$time = NULL
  testing$delta = NULL
  numTimepoint = length(timesplit)
  numReturnNA = 0
  numReturnNoisy = 0
  threshold = 0.05
  
  #first value in cpt
  previousTimepointProb = rep(1,nrow(testing))
  previousProb = rep(1-(1/numTimepoint),nrow(testing))
  previousProb2 = rep(1-(1/numTimepoint),nrow(testing))
  previousProb3 = rep(1-(1/numTimepoint),nrow(testing))
  previousProb4 = rep(1-(1/numTimepoint),nrow(testing))
  previousProb5 = rep(1-(1/numTimepoint),nrow(testing))
  
  survivalFunction <- data.frame(matrix(ncol = nrow(testing), nrow = numTimepoint))
  for(i in 1:numTimepoint) {
    tempFit = fitList[[i]]
    tempNb = nbList[[i]]
    for(j in 1:nrow(testing)) {
      
      if(method == 'lw') {
        eviList = as.list(testing[j,])
        prob = cpquery(tempFit, event = (TIMEPOINT == 0), evidence = as.list(testing[j,]),method = 'lw')
      }else if(method == 'predict') {
        evidence = testing[j,]
        predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "parents", prob = TRUE, n=10000)
        attr(predicted, "prob")
        prob = attr(predicted, "prob")[1]
      }else if(method == 'naivebayes') {
        evidence = testing[j,]
        prob = naivebayesPred(evidence,tempFit,tempNb,childrenList[[i]])
      }else{
        print('No such query method')
      }
      previousProb5[j] = previousProb4[j]
      previousProb4[j] = previousProb3[j]
      previousProb3[j] = previousProb2[j]
      previousProb2[j] = previousProb[j]
      previousProb[j] = prob
      #prob = (previousProb4[j]+previousProb3[j]*2+previousProb2[j]*3+previousProb[j]*4)/10
      
      survivalFunction[i,j] = prob*previousTimepointProb[j]
      
      if(is.na(survivalFunction[i,j])) {
        survivalFunction[i,j] = previousTimepointProb[j]-0.01
        if(survivalFunction[i,j] <0){survivalFunction[i,j] = 0}
        numReturnNA = numReturnNA+1
      }
      if(i==1 & survivalFunction[i,j]<0.6){
        survivalFunction[i,j] = 0.6
        numReturnNoisy = numReturnNoisy+1
      }
      if(previousTimepointProb[j]<=survivalFunction[i,j]){
        survivalFunction[i,j] = previousTimepointProb[j]-0.00001
        if(survivalFunction[i,j] <0){survivalFunction[i,j] = 0}
        numReturnNoisy = numReturnNoisy+1
        #print('probability not decreasing')
      }
      previousTimepointProb[j] = survivalFunction[i,j]
      
    }
  }
  if(numReturnNA>0) {cat(numReturnNA, 'cpquery return NA. impute by previous timepoint')}
  if(numReturnNoisy>0) {cat(numReturnNoisy, 'cpquery return incorrect probability. impute by previous timepoint')}
  colnames(survivalFunction) = 1:nrow(testing)
  survivalFunction = cbind(time = timesplit, survivalFunction) 
  return(survivalFunction)
}

plotDag <-function(dag){
  #plot with bigger fontsize
  g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(dag))
  graph::nodeRenderInfo(g) <- list(fontsize=80)
  Rgraphviz::renderGraph(g)
}
blacklistFunction <- function(nodeNames){
  blackList = matrix(ncol = 2, nrow = length(nodeNames))
  blackList2 = matrix(ncol = 2, nrow = length(nodeNames))
  blackList3 = matrix(ncol = 2, nrow = length(nodeNames))
  for(i in 1:length(nodeNames)) {
    blackList[i,] = c(nodeNames[i],"PREVTIMEPOINT")
    blackList2[i,] = c("TIMEPOINT", nodeNames[i])
    blackList3[i,] = c("PREVTIMEPOINT", nodeNames[i])
  }
  blackList = blackList[blackList[,1] != 'PREVTIMEPOINT',]
  blackList2 = blackList2[blackList2[,2] != 'TIMEPOINT',]
  blackList2 = blackList2[blackList2[,2] != 'PREVTIMEPOINT',]
  blackList3 = blackList3[blackList3[,2] != 'TIMEPOINT',]
  blackList3 = blackList3[blackList3[,2] != 'PREVTIMEPOINT',]
  
  #blackList = rbind(blackList, blackList2)
  #blackList = rbind(blackList, blackList3)
  
  return(blackList2)
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

smoothCurve <-function(survivalCurve) {
  newSurvivalCurve = survivalCurve
  for(i in 3:(nrow(survivalCurve)-2)) {
    for(j in 2:ncol(survivalCurve)) {
      newSurvivalCurve[i,j] = (survivalCurve[i-2,j]+survivalCurve[i-1,j]+2*survivalCurve[i,j]+survivalCurve[i+1,j]+survivalCurve[i+2,j])/6
    }
  }
  for(j in 2:ncol(newSurvivalCurve)) {
    for(i in 1:(nrow(newSurvivalCurve)-1)) {
      if(newSurvivalCurve[i+1,j] > newSurvivalCurve[i,j]){
        print('fix curve')
        newSurvivalCurve[i+1,j] = newSurvivalCurve[i,j]
      }
    }
  }
  return(newSurvivalCurve)
}

overlapCurve <- function(curve1, curve2, curve3, curve4) {
  newCurve = curve1
  i2 = 1
  i3 = 1
  i4 = 1
  for(i in 2:(nrow(curve1)-1)) {
    if(curve1[i,1]>curve2[i2,1]) {i2 = i2+1}
    if(curve1[i,1]>curve3[i3,1]) {i3 = i3+1}
    if(curve1[i,1]>curve4[i4,1]) {i4 = i4+1}
    for(j in 2:ncol(curve1)) {
      newCurve[i,j] = (curve4[i4,j]+curve3[i3,j]+curve2[i2,j]+curve1[i,j])/4
    }
  }
  
  return(newCurve)
}

naivebayesPred <- function(evidence,tempFit,tempNb,childrenList) {
  threshold = 0.05
  predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "parents", prob = TRUE, n=10000)
  ap = attr(predicted, "prob")
  
  p0 = 0
  p1 = 0
  variables = childrenList
  for(v in variables) {
    e = as.integer(evidence[[v]])
    t = tempNb$tables[v]
    
    if(t[[v]][1,e]>threshold) {rawp0=t[[v]][1,e]} else {rawp0=threshold}
    if(t[[v]][2,e]>threshold) {rawp1=t[[v]][2,e]} else {rawp1=threshold}
    p0 = p0 + log(rawp0)
    p1 = p1 + log(rawp1)
  }
  
  
  p0 = p0 + log(ap[1])
  p1 = p1 + log(ap[2])
  p0_Prob = exp(p0)/(exp(p0)+exp(p1))
  
  prob = p0_Prob
  return(prob)
}

