#------------------------------------------------------------
#Cumulative model with seperate network at each timepoint
#add number of timepoint
#structure em
#connected to BayesianNet Upper
#
#------------------------------------------------------------

library(bnlearn)
library(plyr)
library(gtools)
library(e1071)
source("Models/cleanDataForBayesNet.R")

BayesianNet = function(training, testing, BayesianC1 = 1, timesplit = 0,debug = FALSE){
  originalTraining = training
  originalTesting = testing

  queryMethod = 'predict'
  prior = T
  
  if(timesplit == 0) {
    #numTimepoint = floor(sqrt(nrow(training)))
    numTimepoint = 10
    timesplit = timeSplitFunction(numTimepoint = numTimepoint, method = 'quantile', data = training, debug=T)
  }
  print(timesplit)
  
  numTimepoint = length(timesplit)
  
  allData = rbind(testing,training)
  
  cleaned = cleanDataType(nrow(testing),allData,discretize_level=7)
  testing = cleaned$test
  training = cleaned$train
  noCensorTraining = training[training$delta == 1,]
  dataList = createTimeSplit(training,timesplit,includeNa = F)
  dataListNA = createTimeSplit(training,timesplit,includeNa = T)
  
  numTimepoint = length(dataList)
  
  #learning structure and parameter
  fitList <- vector("list", numTimepoint)
  dagList <- vector("list", numTimepoint)
  nbList <- vector("list", numTimepoint)
  childrenList <- vector("list", numTimepoint)
  print('learn starting graph')
  data = training
  discretize_level = 8
  breaks = seq(0,100,100/discretize_level)
  breakList = c(-Inf, unique(quantile(data$time, probs = breaks[2:discretize_level]/100)), Inf)
  data$time = cut(data$time,breaks = breakList,labels = 1:(length(breakList)-1))
  data$delta = NULL
  data$TIMEPOINT = data$time
  data$PREVTIMEPOINT = data$time
  data$time = NULL
  start = hc(data,score='bic')
  
  
  #plotDag(start)
  
  blackList = blacklistFunction(colnames(dataList[[1]]))
  
  #dataList[[floor(numTimepoint/2)]]$PREVTIMEPOINT = NULL
  start = structural.em(dataList[[floor(numTimepoint/2)]], maximize = "hc",maximize.args = list(blacklist=blackList,whitelist=c('PREVTIMEPOINT','TIMEPOINT')), fit = "mle",impute='bayes-lw',return.all = T,start = start, max.iter = 10, debug = FALSE)
  start = start$dag
  plotDag(start)
  #whiteList = c('PREVTIMEPOINT','TIMEPOINT')
  rcount = nrow(training)
  for(i in 1:numTimepoint) {
    if(isTRUE(debug)){cat(i)
    cat(' ')}
    
    data = dataList[[i]]
    #rcount = nrow(data[data$PREVTIMEPOINT == 0,])
    #rcount = nrow(data)

    if(i>1 & i<numTimepoint) {
      structureData = rbind(dataList[[i-1]],dataList[[i]])
      structureData = rbind(structureData,dataList[[i+1]])
    }else if(i==1) {
      structureData = rbind(dataList[[i]],dataList[[i+1]])
    }else if(i==numTimepoint) {
      structureData = rbind(dataList[[i-1]],dataList[[i]])
    }
    structureData = structureData[structureData$PREVTIMEPOINT == 0,]
    #structureData = rbind(structureData,data)
    structureData = data
    
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
    dag = hc(data,restart=500,blacklist=blackList,whitelist = c('PREVTIMEPOINT','TIMEPOINT'),start=start,score='bic')
    dag = cextend(dag, strict = TRUE, debug = FALSE)
    
    #dag = structural.em(dataListNA[[i]], maximize = "hc",maximize.args = list(restart=100,blacklist=blackList,whitelist=c('PREVTIMEPOINT','TIMEPOINT')), fit = "mle",impute='bayes-lw',return.all = T,start = start, max.iter = 1000, debug = FALSE)
    #dag = dag$dag
    #dag = cextend(dag, strict = TRUE, debug = FALSE)
    #dag = start
    dagList[[i]] = dag 
    #start = dag
    fit = bn.fit(dag, data,method='mle',iss = floor(rcount*BayesianC1))
    #fit = bn.fit(dag, data, method='bayes')
    fitList[[i]] <- fit
    nb_model = naiveBayes(as.factor(TIMEPOINT) ~., data=data)
    nbList[[i]] <- nb_model
    childrenList[[i]] <- children(fitList[[i]],'TIMEPOINT')
    #print(parents(fitList[[i]],'TIMEPOINT'))
    #plotDag(dag)
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
      dataPrev = data[data[,'PREVTIMEPOINT']==0,]
      dataPrev$PREVTIMEPOINT = NULL
      
      #pweight = nrow(dataPrev)
      pweight = 10
      w = 10
      prior0 = nrow(dataPrev[dataPrev$TIMEPOINT==0,])/nrow(dataPrev)
      prior1 = nrow(dataPrev[dataPrev$TIMEPOINT==1,])/nrow(dataPrev)
      #print(prior0)
      
      cpt = coef(fitList[[i]]$TIMEPOINT)
      p = colnames(as.data.frame(cpt))
      p = p[p!="Freq"]

      if(length(p)<=2){
        newcpt = cpt
        newcpt[1,2] = 0
        newcpt[2,2] = 1
        #newcpt[1,1] = rbeta(1,floor((prior0*pweight)*w),floor((prior1*pweight)*w))
        newcpt[1,1] = prior0
        newcpt[2,1] = 1 - newcpt[1,1]
        fitList[[i]]$TIMEPOINT = newcpt
      }else{
        newcpt = table(data[,p])
        commaString = commaStringFunction(p)
        codeStringTotal1 = paste0('newcpt[',commaString[1],']',collapse = '')
        codeStringTotal2 = paste0('newcpt[',commaString[2],']',collapse = '')
        #total = eval(parse(text=codeStringTotal1))+eval(parse(text=codeStringTotal2))
        timepoint1 = eval(parse(text=codeStringTotal1))
        timepoint2 = eval(parse(text=codeStringTotal2))
        
        for(j in 1:length(timepoint1)) {
          timepoint1[j] = rbeta(1,(prior0*pweight+timepoint1[j])*w,(prior1*pweight+timepoint2[j])*w)
          #while(timepoint1[j]<0.5){
          #print(timepoint1[j])
          #  timepoint1[j] = rbeta(1,9+timepoint1[j]*5,1+timepoint2[j]*5)
          #}
        }
        
        codeString1 = paste0('newcpt[',commaString[1],']=','timepoint1',collapse = '')
        codeString2 = paste0('newcpt[',commaString[2],']=1-','newcpt[',commaString[1],']',collapse = '')
        eval(parse(text=codeString1))
        eval(parse(text=codeString2))
        
        codeString3 = paste0('newcpt[',commaString[3],']=','0',collapse = '')
        codeString4 = paste0('newcpt[',commaString[4],']=1-','newcpt[',commaString[3],']',collapse = '')
        eval(parse(text=codeString3))
        eval(parse(text=codeString4))
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
    cat(i)
    cat(' ')
    for(j in 1:nrow(testing)) {
      
      if(method == 'lw') {
        eviList = as.list(testing[j,])
        prob = cpquery(tempFit, event = (TIMEPOINT == 0), evidence = as.list(testing[j,]),method = 'lw')
      }else if(method == 'predict') {
        evidence = testing[j,]
        predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "bayes-lw", prob = TRUE, n=500)
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
      #prob = (previousProb3[j]+previousProb2[j]+previousProb[j])/3
      #prob = previousProb[j]

      survivalFunction[i,j] = prob*previousTimepointProb[j]
      #survivalFunction[i,j] = prob
      
      if(is.na(survivalFunction[i,j])) {
        survivalFunction[i,j] = previousTimepointProb[j]-0.01
        if(survivalFunction[i,j] <0){survivalFunction[i,j] = 0}
        numReturnNA = numReturnNA+1
      }
      
      if(previousTimepointProb[j]<survivalFunction[i,j]){
        survivalFunction[i,j] = previousTimepointProb[j]-0.0001
        if(survivalFunction[i,j] <0){survivalFunction[i,j] = 0}
        numNotDecreasing = numNotDecreasing+1
        #print('probability not decreasing')
      }
      
      previousTimepointProb[j] = survivalFunction[i,j]
      
    }
  }
  if(numReturnNA>0) {cat('return NA: ',numReturnNA)}
  if(numNotDecreasing>0) {cat('Not decreasing: ',numNotDecreasing)}
  #survivalFunction = rbind(rep(1,nrow(testing)),survivalFunction)
  
  colnames(survivalFunction) = 1:nrow(testing)
  #survivalFunction = cbind(time = c(0,timesplit), survivalFunction) 
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
    blackList2[i,] = c("TIMEPOINT",nodeNames[i])
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
  commaString1 = paste0(replace(commaString1,commaString1=="PREVTIMEPOINT",'1'),collapse = ',')
  commaString1 = paste0(commaString1,collapse = ',')
  commaString2 = replace(p,p!="TIMEPOINT"&p!="PREVTIMEPOINT",'')
  commaString2 = replace(commaString2,commaString2=="TIMEPOINT",'2')
  commaString2 = paste0(replace(commaString2,commaString2=="PREVTIMEPOINT",'1'),collapse = ',')
  commaString2 = paste0(commaString2,collapse = ',')
  commaString3 = replace(p,p!="TIMEPOINT"&p!="PREVTIMEPOINT",'')
  commaString3 = replace(commaString3,commaString3=="TIMEPOINT",'1')
  commaString3 = paste0(replace(commaString3,commaString3=="PREVTIMEPOINT",'2'),collapse = ',')
  commaString4 = replace(p,p!="TIMEPOINT"&p!="PREVTIMEPOINT",'')
  commaString4 = replace(commaString4,commaString4=="TIMEPOINT",'2')
  commaString4 = paste0(replace(commaString4,commaString4=="PREVTIMEPOINT",'2'),collapse = ',')
  return(list(commaString1,commaString2,commaString3,commaString4))
}

fixtime <- function(curve) {
  newcurve = curve
  newcurve[1,1] = curve[1,1]/2
  for(i in 2:nrow(curve)) {
    newcurve[i,1] = (curve[i-1,1] + curve[i,1])/2
  }
  return(newcurve)
}

smoothCurve <-function(survivalCurve) {
  newSurvivalCurve = survivalCurve
  for(i in 2:(nrow(survivalCurve)-1)) {
    for(j in 2:ncol(survivalCurve)) {
      newSurvivalCurve[i,j] = (survivalCurve[i-1,j]+survivalCurve[i,j]+survivalCurve[i+1,j])/3
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