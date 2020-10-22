#------------------------------------------------------------
#Cumulative model with seperate network at each timepoint
#add effective smaple size correction
#
#
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
  
  queryMethod = 'naivebayes'
  prior = T
  
  if(timesplit == 0) {
    #numTimepoint = floor(sqrt(nrow(training)))
    numTimepoint = 10
    timesplit = timeSplitFunction(numTimepoint = numTimepoint, method = 'quantile', data = training, debug=T)
  }
  print(timesplit)
  
  numTimepoint = length(timesplit)
  
  allData = rbind(testing,training)
  
  cleaned = cleanDataType(nrow(testing),allData,discretize_level=5)
  testing = cleaned$test
  training = cleaned$train
  noCensorTraining = training[training$delta == 1,]
  imputeTraining = training
  maxTime = max(imputeTraining$time)
  for(i in 1:nrow(training)) {
    if(imputeTraining[i,]$delta == 0) {
      imputeTime = runif(1, imputeTraining[i,]$time, maxTime)
      imputeTraining[i,]$time = imputeTime
      imputeTraining[i,]$delta = 1
    }
  }
  dataList = createTimeSplit(training,timesplit,assumption = F,includeNa = F)
  #dataListNA = createTimeSplit(training,timesplit,includeNa = T)
  print(length(dataList))
  
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
  blackliststart = blacklistFunctionStart(colnames(dataList[[1]]))
  start = hc(data,score='bic')
  
  
  #plotDag(start)
  
  blackList = blacklistFunction(colnames(dataList[[1]]))
  
  #dataList[[floor(numTimepoint/2)]]$PREVTIMEPOINT = NULL
  data = dataList[[floor(numTimepoint/2)]]
  data$PREVTIMEPOINT = NULL
  start = structural.em(data, maximize = "hc",maximize.args = list(blacklist=NULL,whitelist=NULL), fit = "mle",impute='bayes-lw',return.all = T,start = NULL, max.iter = 10, debug = FALSE)
  start = start$dag
  #plotDag(start)

  allmb = c()
  #whiteList = c('PREVTIMEPOINT','TIMEPOINT')
  for(i in 1:numTimepoint) {
    if(isTRUE(debug)){cat(i)
      cat(' ')}
    
    data = dataList[[i]]
    rcount = nrow(data[data$PREVTIMEPOINT == 0,])
    #rcount = nrow(data)
    data = data[data$PREVTIMEPOINT == 0,]
    data$PREVTIMEPOINT = NULL
    
    # if(i>1 & i<numTimepoint) {
    #   structureData = rbind(dataList[[i-1]],dataList[[i]])
    #   structureData = rbind(structureData,dataList[[i+1]])
    # }else if(i==1) {
    #   structureData = rbind(dataList[[i]],dataList[[i+1]])
    # }else if(i==numTimepoint) {
    #   structureData = rbind(dataList[[i-1]],dataList[[i]])
    # }
    structureData = dataList[[i]]
    structureData$PREVTIMEPOINT = NULL
    structureData = structureData[complete.cases(structureData),]
    #structureData = rbind(structureData,data)
    #structureData = data
    

    dag = hc(structureData,restart=50,blacklist=NULL,whitelist = NULL,start=NULL,score='bic')
    #dag = cextend(dag, strict = TRUE, debug = FALSE)
    #whitelist = rbind(c('LDH_SERUM','TIMEPOINT'),c('GRANULOCYTES','TIMEPOINT'))
    #dag = structural.em(structureData, maximize = "hc",maximize.args = list(restart=10,blacklist=NULL,whitelist=whitelist), fit = "bayes",impute='bayes-lw',return.all = T,start = start, max.iter = 10, debug = FALSE)
    #dag = dag$dag
    dag = cextend(dag, strict = TRUE, debug = FALSE)
    #dag = start
    #dagList[[i]] = dag 
    #start = dag
    fit = bn.fit(dag, data,method='bayes',iss=50)
    #fit = bn.fit(dag, data, method='bayes')
    fitList[[i]] <- fit
    #data = dataListNA[[i]]
    nbdata=data
    #nbdata$PREVTIMEPOINT = NULL
    nb_model = naiveBayes(as.factor(TIMEPOINT) ~., nbdata, laplace = 20)
    nbList[[i]] <- nb_model
    childrenList[[i]] <- children(fitList[[i]],'TIMEPOINT')
    #print(parents(fitList[[i]],'TIMEPOINT'))
    allmb = unique(c(allmb,mb(fitList[[i]],'TIMEPOINT')))
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
        
        timepointNA = timepoint1
        for(k in 1:length(timepointNA)){timepointNA[k]=0}
        for(k in 1:nrow(data)) {
          if(is.na(data[k,'TIMEPOINT'])) {
            
            index = data[k,p[-1]]
            value = 
            codeStringTotalNA = paste0('timepointNA[',toString(index),']=timepointNA[',toString(index),']+1',collapse = '')
            eval(parse(text=codeStringTotalNA))
            }
        }
        
        numData = nrow(data)
        weight = 1

        alive = nrow(data[data$TIMEPOINT==0 & complete.cases(data),])
        dead = nrow(data[data$TIMEPOINT==1 & complete.cases(data),])
        cencered = nrow(data[is.na(data$TIMEPOINT),])
        priorAlive = 1 - (dead/(alive+dead+0.5*cencered))
        nPriorAlive = floor(priorAlive*100)
        nPriorDeath = floor((1-priorAlive)*100)
        
        timepointProb = timepoint1
        for(j in 1:length(timepoint1)) {
          #timepoint1[j] = rbeta(1,(numTimepoint-1)*weight*15+timepoint1[j]*5,1*weight*15+timepoint2[j]*5)
          survivalRate = 1 - (timepoint2[j])/(timepoint1[j] + timepoint2[j] + timepointNA[j]*0.5)
          if(is.nan(survivalRate)) {
            survivalRate = priorAlive
            print('divided by zero')
          }
          nalive = floor(survivalRate*100)
          ndeath = floor((1-survivalRate)*100)
          
          #prob = priorAlive
          prob = rbeta(1,nPriorAlive*weight+nalive*10,nPriorDeath*weight+ndeath*10)
          timepointProb[j] = prob
        }
        
        codeString1 = paste0('newcpt[',commaString[1],']=','timepointProb',collapse = '')
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
  survivalFunctionTesting = predictFunction(fitList,nbList,childrenList,testing,timesplit,method = queryMethod,dataList)
  survivalFunctionTraining = predictFunction(fitList,nbList,childrenList,training,timesplit,method = queryMethod,dataList)
  
  testCurvesToReturn = survivalFunctionTesting
  timesAndCensTest = cbind.data.frame(time = originalTesting$time, delta = originalTesting$delta)
  timesAndCensTrain = cbind.data.frame(time = originalTraining$time, delta = originalTraining$delta)
  trainingCurvesToReturn = survivalFunctionTraining
  
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}

predictFunction <- function(fitList,nbList,childrenList,testing,timesplit,method,dataList) {
  testing$time = NULL
  testing$delta = NULL
  #testing$PREVTIMEPOINT = factor(integer(nrow(testing)),levels = c('0','1'))
  numTimepoint = length(timesplit)
  numReturnNA = 0
  numNotDecreasing = 0
  threshold = 0.000001
  
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
    
    data = dataList[[i]]
    data = data[data$PREVTIMEPOINT == 0,]    
    alive = nrow(data[data$TIMEPOINT==0 & complete.cases(data),])
    dead = nrow(data[data$TIMEPOINT==1 & complete.cases(data),])
    cencered = nrow(data[is.na(data$TIMEPOINT),])
    total = alive + dead + 0.5*cencered
    
    for(j in 1:nrow(testing)) {
      
      if(method == 'lw') {
        eviList = as.list(testing[j,])
        prob = cpquery(tempFit, event = (TIMEPOINT == 0), evidence = as.list(testing[j,]),method = 'lw')
      }else if(method == 'ls') {
        evi = generateEvidenceString(colnames(testing),testing[j,])
        evi <<- evi
        prob = cpquery(tempFit, event = (TIMEPOINT == 0), evidence = eval(parse(text=evi)),n = 10000,method = 'ls')
        #print(evi)
        #prob = cpquery(tempFit, event = (TIMEPOINT == 0), evidence = (BOX1_SCORE==-0.847670404481496)&(BOX2_SCORE==0.0849319942719063)&(BOX3_SCORE==4)&(PERFORMANCE_STATUS==-0.221318340332475)&(BMI==7)&(NO_PROBLEM==1.08932317940723)&(SITE_BRUNCHUS_LUNG==1.62055863299408),n = 1000,method = 'ls')
        
      }else if(method == 'predict') {
        evidence = testing[j,]
        predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "parents", prob = TRUE, n=500)
        attr(predicted, "prob")
        prob = attr(predicted, "prob")[1]
      }else if(method == 'naivebayes') {
        evidence = testing[j,]
        predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "parents", prob = TRUE)
        ap = attr(predicted, "prob")

        p0 = 1
        p1 = 1
        variables = childrenList[[i]]
        for(v in variables) {
          e = as.integer(evidence[[v]])
          t = tempNb$tables[v]
          t = tempFit[[v]][['prob']]
          temp_evidence = evidence
          temp_evidence[[v]] = NULL
          temp_evidence$TIMEPOINT = factor(0,levels = c('0','1'))
          t1 = predict(tempFit, node=v, temp_evidence,method = "parents", prob = TRUE)
          t1 = attr(t1, "prob")
          temp_evidence$TIMEPOINT = factor(1,levels = c('0','1'))
          t2 = predict(tempFit, node=v, temp_evidence,method = "parents", prob = TRUE)
          t2 = attr(t2, "prob")
          
          if(t1[e]>threshold) {rawp0=t1[e]} else {rawp0=threshold}
          if(t2[e]>threshold) {rawp1=t2[e]} else {rawp1=threshold}
          p0 = p0 * rawp0
          p1 = p1 * rawp1
        }
        #print(t)
        p0 = p0 * ap[1]
        p1 = p1 * ap[2]
        p0_Prob = p0/(p0+p1)
        
        prob = p0_Prob
      }else{
        print('No such query method')
      }
      previousProb5[j] = previousProb4[j]
      previousProb4[j] = previousProb3[j]
      previousProb3[j] = previousProb2[j]
      previousProb2[j] = previousProb[j]
      previousProb[j] = prob
      #prob = (previousProb3[j]+previousProb2[j]*2+previousProb[j]*3)/3
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
    blackList2[i,] = c(nodeNames[i],"TIMEPOINT")
    blackList3[i,] = c("PREVTIMEPOINT", nodeNames[i])
  }
  blackList = blackList[blackList[,1] != 'PREVTIMEPOINT',]
  blackList2 = blackList2[blackList2[,1] != 'TIMEPOINT',]
  blackList2 = blackList2[blackList2[,1] != 'PREVTIMEPOINT',]
  blackList3 = blackList3[blackList3[,2] != 'TIMEPOINT',]
  blackList3 = blackList3[blackList3[,2] != 'PREVTIMEPOINT',]
  
  #blackList = rbind(blackList, blackList2)
  #blackList = rbind(blackList, blackList3)
  
  return(blackList2)
}

blacklistFunctionStart <- function(nodeNames){
  blackList = matrix(ncol = 2, nrow = length(nodeNames))
  blackList2 = matrix(ncol = 2, nrow = length(nodeNames))
  blackList3 = matrix(ncol = 2, nrow = length(nodeNames))
  blackList4 = matrix(ncol = 2, nrow = length(nodeNames))
  for(i in 1:length(nodeNames)) {
    blackList[i,] = c(nodeNames[i],"PREVTIMEPOINT")
    blackList2[i,] = c("TIMEPOINT",nodeNames[i])
    blackList3[i,] = c("PREVTIMEPOINT", nodeNames[i])
    blackList4[i,] = c(nodeNames[i],"TIMEPOINT")
  }
  blackList = blackList[blackList[,1] != 'PREVTIMEPOINT',]
  blackList2 = blackList2[blackList2[,2] != 'TIMEPOINT',]
  blackList2 = blackList2[blackList2[,2] != 'PREVTIMEPOINT',]
  blackList3 = blackList3[blackList3[,2] != 'TIMEPOINT',]
  blackList3 = blackList3[blackList3[,2] != 'PREVTIMEPOINT',]
  
  blackList = rbind(blackList, blackList2)
  blackList = rbind(blackList, blackList3)
  blackList = rbind(blackList, blackList4)
  
  return(blackList)
}

indexingFunction <- function(index) {
  
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

generateEvidenceString <- function(nodeNames,dataInstance) {
  colnamelist = nodeNames
  for(k in 1:length(colnamelist)) {
    if(class(dataInstance[,k]) == 'factor') {
      colnamelist[k] = paste(colnamelist[k], dataInstance[1,k], sep="==")
    }else {
      colnamelist[k] = paste(colnamelist[k], dataInstance[1,k], sep="=")
    }
  }
  evi = paste(colnamelist,collapse = ')&(')
  evi = paste( c( '(' , evi , ')' ), collapse = '')
  
  return(evi)
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