#------------------------------------------------------------
#Cumulative model with seperate network at each timepoint
#add number of timepoint
#structure em
#independent t_i
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
  
  numTimepoint = 10
  queryMethod = 'predict'
  prior = F
  
  timesplit = timeSplitFunction(numTimepoint = numTimepoint, method = 'quantile', data = training, debug=T)
  numTimepoint = length(timesplit)
  
  allData = rbind(testing,training)
  
  cleaned = cleanDataType(nrow(testing),allData,discretize_level=5)
  testing = cleaned$test
  training = cleaned$train
  noCensorTraining = training[training$delta == 1,]
  dataList = createTimeSplit(training,timesplit,includeNa = F)
  
  #learning structure and parameter
  fitList <- vector("list", numTimepoint)
  dagList <- vector("list", numTimepoint)
  nbList <- vector("list", numTimepoint)
  childrenList <- vector("list", numTimepoint)
  print('learn starting graph')
  data = training
  discretize_level = 10
  breaks = seq(0,100,100/discretize_level)
  breakList = c(-Inf, unique(quantile(data$time, probs = breaks[2:discretize_level]/100)), Inf)
  data$time = cut(data$time,breaks = breakList,labels = 1:(length(breakList)-1))
  data$delta = NULL
  data$TIMEPOINT = data$time
  #data$PREVTIMEPOINT = data$time
  data$time = NULL
  start = hc(data,score='bic')
  
  #plotDag(start)
  
  blackList = blacklistFunction(colnames(dataList[[1]]))
  
  #whiteList = c('PREVTIMEPOINT','TIMEPOINT')
  for(i in 1:numTimepoint) {
    print(i)
    
    data = dataList[[i]]
    # if(i<numTimepoint) {
    #   data$NEXT = dataList[[i+1]]$TIMEPOINT
    # }
    
    #data = data[data$PREVTIMEPOINT == 0,]
    # if(i>2 & i < numTimepoint-1){
    #   datap = dataList[[i-1]]
    #   datap = datap[datap$PREVTIMEPOINT == 0,]
    #   datan = dataList[[i+1]]
    #   datan = datan[datan$PREVTIMEPOINT == 0,]
    #   data = rbind(datap,dataList[[i]])
    #   data = rbind(data,datan)
    # }else if(i==1){
    #   datan = dataList[[i+1]]
    #   datan = datan[datan$PREVTIMEPOINT == 0,]
    #   data = rbind(dataList[[i]],datan)
    # }else if(i==numTimepoint){
    #   datap = dataList[[i-1]]
    #   datap = datap[datap$PREVTIMEPOINT == 0,]
    #   data = rbind(datap,dataList[[i]])
    # }
    #data = data[data[,'PREVTIMEPOINT']==0,]
    #data$PREVTIMEPOINT = NULL
    #numMissing = nrow(data)/numTimepoint - nrow(data[data$TIMEPOINT==1,])
    #print(numMissing)
    #print(nrow(data[data$TIMEPOINT==0,])/(nrow(data[data$TIMEPOINT==1,])+nrow(data[data$TIMEPOINT==0,])))
    #censoredData = training[originalTraining$delta == 0 && originalTraining$time<timesplit[i]]
    #topname = independent_test(data)
    #topname = c(topname,'PREVTIMEPOINT')
    #whiteList = matrix(ncol = 2, nrow = length(topname))
    #for(j in 1:length(topname)) {
    #  whiteList[j,] = c(topname[j],"TIMEPOINT")
    #}
    dag = hc(data,start=NULL,score='bic')
    #dag = cextend(dag, strict = TRUE, debug = FALSE)
    
    #dag = structural.em(data, maximize = "hc",maximize.args = list(restart = 5,blacklist=NULL), fit = "mle",impute='bayes-lw',return.all = T,start = NULL, max.iter = 5, debug = FALSE)
    #dag = dag$dag
    dag = cextend(dag, strict = TRUE, debug = FALSE)
    dagList[[i]] = dag 
    #start = dag
    fit = bn.fit(dag, data, method='mle')
    fitList[[i]] <- fit
    nb_model = naiveBayes(as.factor(TIMEPOINT) ~., data=data)
    nbList[[i]] <- nb_model
    childrenList[[i]] <- children(fitList[[i]],'TIMEPOINT')
    #print(parents(fitList[[i]],'TIMEPOINT'))
    plotDag(dag)
  }
  
  # for(i in 1:length(dataList)) {
  #   
  #   data = dataList[[i]]
  #   data$PREVTIMEPOINT = NULL
  #   if(i<numTimepoint) {
  #     start2 = dagList[[i+1]]
  #   }else {
  #     start2 = start
  #   }
  #   dag = structural.em(data, maximize = "hc",maximize.args = list(blacklist=NULL), fit = "mle",impute='bayes-lw',return.all = T,start = start2, max.iter = 5, debug = FALSE)
  #   dag = dag$dag
  #   dag = cextend(dag, strict = TRUE, debug = FALSE)
  #   dagList[[i]] = dag
  #   fit = bn.fit(dag, data, method='mle')
  #   fitList[[i]] <- fit
  # }
  
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
  
  testCurvesToReturn = survivalFunctionTesting
  timesAndCensTest = cbind.data.frame(time = originalTesting$time, delta = originalTesting$delta)
  timesAndCensTrain = cbind.data.frame(time = originalTraining$time, delta = originalTraining$delta)
  trainingCurvesToReturn = survivalFunctionTraining
  
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}

predictFunction <- function(fitList,nbList,childrenList,testing,timesplit,method) {
  testing$time = NULL
  testing$delta = NULL
  testing$PREVTIMEPOINT = factor(integer(nrow(testing)),levels = c('0','1'))
  numTimepoint = length(timesplit)
  numReturnNA = 0
  numTooLow = 0
  numNotDecreasing = 0
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
    print(i)
    for(j in 1:nrow(testing)) {
      
      if(method == 'lw') {
        eviList = as.list(testing[j,])
        prob = cpquery(tempFit, event = (TIMEPOINT == 0), evidence = as.list(testing[j,]),method = 'lw')
      }else if(method == 'predict') {
        evidence = testing[j,]
        predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "bayes-lw", prob = TRUE, n=5000)
        #attr(predicted, "prob")
        prob = attr(predicted, "prob")[1]
      }else if(method == 'naivebayes') {
        evidence = testing[j,]
        predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "parents", prob = TRUE, n=1000)
        ap = attr(predicted, "prob")
        
        p0 = 0
        p1 = 0
        variables = childrenList[[i]]
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
      }else{
        print('No such query method')
      }
      previousProb5[j] = previousProb4[j]
      previousProb4[j] = previousProb3[j]
      previousProb3[j] = previousProb2[j]
      previousProb2[j] = previousProb[j]
      previousProb[j] = prob
      #prob = (previousProb3[j]+previousProb2[j]*2+previousProb[j]*3)/6
      #prob = previousProb[j]
      
      survivalFunction[i,j] = prob*previousTimepointProb[j]
      #survivalFunction[i,j] = prob
      
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
  blackList = rbind(blackList, blackList3)
  
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