#------------------------------------------------------------
#Cumulative model with seperate network at each timepoint
#add effective smaple size correction
#add weighted censored data
#
#
#------------------------------------------------------------

library(bnlearn)
library(plyr)
library(gtools)
library(e1071)
library(prodlim)
library(MTLR)
library(Rcpp)
source("Models/cleanDataForBayesNet.R")

BayesianNet = function(training, testing, BayesianC1 = 1, timesplit = 0,debug = FALSE){
  originalTraining = training
  originalTesting = testing

  queryMethod = 'naivebayes'
  prior = T
  weighted = T
  
  print(timesplit)
  
  numTimepoint = length(timesplit)
  
  allData = rbind(testing,training)
  
  cleaned = cleanDataType(nrow(testing),allData,discretize_level=5)
  testing = cleaned$test
  training = cleaned$train
  noCensorTraining = training[training$delta == 1,]

  dataList = createTimeSplit(training,timesplit,assumption = F,includeNa = F)
  dataListNA = createTimeSplit(training,timesplit,includeNa = T)
  print(length(dataList))

  numTimepoint = length(dataList)
  
  fitList <- vector("list", numTimepoint)
  dagList <- vector("list", numTimepoint)
  nbList <- vector("list", numTimepoint)
  childrenList <- vector("list", numTimepoint)
  print('learn starting graph')

  allmb = c()
  #whiteList = c('PREVTIMEPOINT','TIMEPOINT')
  for(iter in 1:1) {
  
  for(i in 1:numTimepoint) {
    if(isTRUE(debug)){cat(i)
      cat(' ')}
    
    data = dataList[[i]]
    rcount = nrow(data[data$PREVTIMEPOINT == 0,])
    #rcount = nrow(data)
    data = data[data$PREVTIMEPOINT == 0,]
    data$PREVTIMEPOINT = NULL
    data$time = NULL
    data$delta = NULL
    data$id = NULL

    structureData = dataList[[i]]
    structureData$PREVTIMEPOINT = NULL
    structureData$time = NULL
    structureData$delta = NULL
    structureData$id = NULL
    
    #structureData = structureData[complete.cases(structureData),]
    #structureData = rbind(structureData,data)
    #structureData = data
    
    #whitelist = rbind(c('LDH_SERUM','TIMEPOINT'),c('GRANULOCYTES','TIMEPOINT'))
    #dag = hc(structureData,restart=50,blacklist=NULL,whitelist = NULL,start=NULL,score='bic')
    dag = structural.em(structureData, maximize = "hc",maximize.args = list(restart=20,blacklist=NULL,whitelist=NULL), fit = "bayes",fit.args=list(iss=10),impute='bayes-lw',return.all = T,start = fitList[[i]], max.iter = 20, debug = FALSE)
    emResult = structural.em(data, maximize = "hc",maximize.args = list(restart=20,blacklist=NULL,whitelist=NULL), fit = "bayes",fit.args=list(iss=10),impute='bayes-lw',return.all = T,start = fitList[[i]], max.iter = 20, debug = FALSE)
    #dag = emResult$dag
    #dag = cextend(dag, strict = TRUE, debug = FALSE)
    #dag = start
    #dagList[[i]] = dag 
    #start = dag
    fit = bn.fit(dag$dag, emResult$imputed,method='bayes',iss=20)
    #fit = emResult$fitted
    #fit = bn.fit(dag, data, method='bayes')
    fitList[[i]] <- fit
    #data = dataListNA[[i]]
    nbdata=data
    #nbdata$PREVTIMEPOINT = NULL
    nb_model = naiveBayes(as.factor(TIMEPOINT) ~., nbdata, laplace = 20)
    nbList[[i]] <- nb_model
    childrenList[[i]] <- children(fitList[[i]],'TIMEPOINT')
    #print(children(fitList[[i]],'TIMEPOINT'))
    allmb = unique(c(allmb,mb(fitList[[i]],'TIMEPOINT')))
    #plotDag(dag)
  }
    
  if(prior) {
    fitList = timepointPrior(dataList,fitList,training,originalTraining,timesplit)
  }
    
  
    if(weighted) {
      #fitList = weightedLearning(dataList,dataListNA,fitList,timesplit,training,childrenList,first=F)
      if(iter==1) {
        fitList = weightedLearning(dataList,dataListNA,fitList,timesplit,training,childrenList,first=T,originalTraining)
      }else {
        fitList = weightedLearning(dataList,dataListNA,fitList,timesplit,training,childrenList,first=F,originalTraining)
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
    censored = nrow(data[is.na(data$TIMEPOINT),])
    total = alive + dead + 0.5*censored
    
    for(j in 1:nrow(testing)) {
      
      if(method == 'lw') {
        eviList = as.list(testing[j,])
        prob = cpquery(tempFit, event = (TIMEPOINT == 0), evidence = as.list(testing[j,]),method = 'lw')
      }else if(method == 'ls') {
        
        
        iss_total = (alive+dead+censored*0.5)*0.1
        iss_dead = dead*0.1
        prob = 1-(iss_dead/iss_total)
      }else if(method == 'predict') {
        evidence = testing[j,]
        predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "bayes-lw", prob = TRUE)
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
      #prob = (previousProb3[j]+previousProb2[j]*2+previousProb[j]*3)/6
      #prob = previousProb[j]
      
      #prob = 1-(dead/total)

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

predictSingle <- function(fitList,childrenList,timesplit,evidence,time,method='naivebayes') {
  evidence$time = NULL
  evidence$TIMEPOINT = NULL
  evidence$PREVTIMEPOINT = NULL
  evidence$delta = NULL
  #testing$PREVTIMEPOINT = factor(integer(nrow(testing)),levels = c('0','1'))
  threshold = 0.000001
  
  #first value in cpt
  previousTimepointProb = 1
  
  survivalFunction = rep(0,length(fitList))
  for(i in 1:length(fitList)) {
    tempFit = fitList[[i]]
    
    if(method == 'lw') {
      eviList = as.list(evidence)
      prob = cpquery(tempFit, event = (TIMEPOINT == 0), evidence = as.list(testing[j,]),method = 'lw')
    }else if(method == 'predict') {
      predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "parents", prob = TRUE, n=500)
      attr(predicted, "prob")
      prob = attr(predicted, "prob")[1]
    }else if(method == 'naivebayes') {
      predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "parents", prob = TRUE)
      ap = attr(predicted, "prob")
      
      p0 = 1
      p1 = 1
      variables = childrenList[[i]]
      for(v in variables) {
        e = as.integer(evidence[[v]])
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
    
    survivalFunction[i] = prob*previousTimepointProb
    if(time<timesplit[i]) {
      toReturn = previousTimepointProb + (survivalFunction[i]-previousTimepointProb) * (time - timesplit[i-1])/(timesplit[i] - timesplit[i-1])
      return(toReturn)
    }
    
    
    if(is.na(survivalFunction[i])) {
      survivalFunction[i] = previousTimepointProb-0.01
      if(survivalFunction[i] <0){survivalFunction[i] = 0}
    }
    
    if(previousTimepointProb<survivalFunction[i]){
      survivalFunction[i] = previousTimepointProb[j]-0.0001
      if(survivalFunction[i] <0){survivalFunction[i] = 0}
      #print('probability not decreasing')
    }
    
    previousTimepointProb = survivalFunction[i]
  }
  return(survivalFunction[length(fitList)])
}

timepointPrior <- function(dataList,fitList,training,originalTraining,timesplit) {
  print('Impose prior to TIMEPOINT')
  kmMod = prodlim(Surv(time,delta)~1, data = training)
  C1 = NULL
  if(is.null(C1)){
    C1 = mtlr_cv(Surv(time,delta)~.,data=originalTraining, loss= "conc", C1_vec = c(0.001,0.01,0.1,1,10,100,1000)
                 , train_biases = F,train_uncensored = F)$best_C1
  }
  mod = mtlr(Surv(time,delta)~., data = originalTraining, C1=C1, train_biases = F, train_uncensored = F)
  testCurvesToReturn = predict(mod)
  
  for(i in 1:length(dataList)) {
    data = dataList[[i]]
    #data = data[data[,'PREVTIMEPOINT']==0,]
    data = data[data[,'PREVTIMEPOINT']==0 | is.na(data[,'PREVTIMEPOINT']),]
    #data$PREVTIMEPOINT = NULL
    data$time = NULL
    data$delta = NULL
    #data$id = NULL
    dataComplete = data[complete.cases(data),]
    cpt = coef(fitList[[i]]$TIMEPOINT)
    p = colnames(as.data.frame(cpt))
    p = p[p!="Freq"]
    if(p != "Var1"){
      newcpt = table(dataComplete[,p])
      
      commaString = commaStringFunction(p)
      codeStringTotal1 = paste0('newcpt[',commaString[1],']',collapse = '')
      codeStringTotal2 = paste0('newcpt[',commaString[2],']',collapse = '')
      #total = eval(parse(text=codeStringTotal1))+eval(parse(text=codeStringTotal2))
      timepoint1 = eval(parse(text=codeStringTotal1))
      timepoint2 = eval(parse(text=codeStringTotal2))
      
      timepointNA = timepoint1
      timepointNATotal = timepoint1
      for(k in 1:length(timepointNA)){timepointNA[k]=0}
      for(k in 1:length(timepointNA)){timepointNATotal[k]=0}
      for(k in 1:nrow(data)) {
        #if(is.na(data[k,'TIMEPOINT'])) {
        if(T) {
          id = data[k,'id']
          curveSpline = splinefun(testCurvesToReturn[,1],testCurvesToReturn[,id+1],method='hyman')
          
          if(i>1) {
            #weight = (data[k,'time']-timesplit[i-1])/(timesplit[i] - timesplit[i-1])
            weight = curveSpline(timesplit[i])/curveSpline(timesplit[i-1])
            #weight = predict(kmMod,timesplit[i])/predict(kmMod,timesplit[i-1])
          }else {
            #weight = data[k,'time']/timesplit[i]
            weight = curveSpline(timesplit[i])
            #weight = predict(kmMod,timesplit[i])
          }
          if(weight<0 | is.na(weight)){weight = 0}
          if(curveSpline(timesplit[i])<0){weight = 0}
          weight = 1 - weight
          # if(is.na(data[k,'PREVTIMEPOINT'])) {
          #   weight = weight*0.1
          # }else {
          #   weight = weight*0.5
          # }
          
          index = data[k,p[-1]]
          flatIndex = convertToFlatIndex(index,dim(timepointNA))
          timepointNA[flatIndex] = timepointNA[flatIndex] + weight
          flatIndex = convertToFlatIndex(index,dim(timepointNATotal))
          timepointNATotal[flatIndex] = timepointNATotal[flatIndex] + 1
          
        }
      }
      
      # timepointNAdead = timepoint1
      # for(k in 1:length(timepointNA)){timepointNAdead[k]=0}
      # for(k in 1:nrow(data)) {
      #   if(is.na(data[k,'TIMEPOINT'])) {
      #     
      #     if(i>1) {
      #       #weight = (data[k,'time']-timesplit[i-1])/(timesplit[i] - timesplit[i-1])
      #       weight = (predict(kmMod,data[k,'time'])-predict(kmMod,timesplit[i]))/(predict(kmMod,timesplit[i-1]) - predict(kmMod,timesplit[i]))
      #     }else {
      #       #weight = data[k,'time']/timesplit[i]
      #       weight = (predict(kmMod,data[k,'time'])-predict(kmMod,timesplit[i]))/(predict(kmMod,0) - predict(kmMod,timesplit[i]))
      #     }
      #     index = data[k,p[-1]]
      #     codeStringTotalNA = paste0('timepointNAdead[',toString(index),']=timepointNAdead[',toString(index),']+1-weight',collapse = '')
      #     eval(parse(text=codeStringTotalNA))
      #   }
      # }
      numData = nrow(data)
      priorWeight = 0
      
      alive = nrow(data[data$TIMEPOINT==0 & complete.cases(data),])
      dead = nrow(data[data$TIMEPOINT==1 & complete.cases(data),])
      censored = nrow(data[is.na(data$TIMEPOINT),])
      
      iss = 20
      iss_total = iss
      iss_dead = iss * (dead/(alive+dead+censored*0.5))
      
      timepointProb = timepoint1
      for(j in 1:length(timepoint1)) {
        #timepoint1[j] = rbeta(1,(numTimepoint-1)*weight*15+timepoint1[j]*5,1*weight*15+timepoint2[j]*5)
        survivalRate = 1 - (timepointNA[j]+iss_dead)/(timepointNATotal[j]+iss_total)
        #survivalRate = 1 - (timepoint2[j])/(timepoint1[j] + timepoint2[j] + 0.5*timepointNATotal[j])
        #survivalRate = 1-(iss_dead/iss_total)
        if(is.nan(survivalRate)) {
          survivalRate = 1-(iss_dead/iss_total)
          #print('divided by zero')
        }
        prob = survivalRate
        timepointProb[j] = prob
      }
      
      codeString1 = paste0('newcpt[',commaString[1],']=','timepointProb',collapse = '')
      codeString2 = paste0('newcpt[',commaString[2],']=1-','newcpt[',commaString[1],']',collapse = '')
      eval(parse(text=codeString1))
      eval(parse(text=codeString2))
      fitList[[i]]$TIMEPOINT = newcpt
    }
    
  }
  return(fitList)
}

weightedLearning <- function(dataList,dataListNA,fitList,timesplit,training,childrenList,first,originalTraining) {
  kmMod = prodlim(Surv(time,delta)~1, data = training)
  
  C1 = NULL
  if(is.null(C1)){
    C1 = mtlr_cv(Surv(time,delta)~.,data=originalTraining, loss= "conc", C1_vec = c(0.001,0.01,0.1,1,10,100,1000)
                 , train_biases = F,train_uncensored = F)$best_C1
  }
  mod = mtlr(Surv(time,delta)~., data = originalTraining, C1=C1, train_biases = F, train_uncensored = F)
  testCurvesToReturn = predict(mod)
  
  iss = 20
  print('Weighted fitting')
  for(i in 1:length(dataList)) {
    for(node in childrenList[[i]]) {
      data = dataList[[i]]
      #data = data[complete.cases(data),]
      #data = data[data['PREVTIMEPOINT']==0,]
      data = data[data['PREVTIMEPOINT']==0 | is.na(data[,'PREVTIMEPOINT']),]
      #data$PREVTIMEPOINT = NULL
      
      cpt = coef(fitList[[i]][[node]])
      p = colnames(as.data.frame(cpt))
      p = p[p!="Freq"]
      #iss = nrow(data)
      newcpt = cpt

      for(k in 1:length(newcpt)){newcpt[k]=0.0}
      for(k in 1:nrow(data)) {
        if(T) {
        #if(is.na(data[k,'TIMEPOINT'])) {
          
          
          id = data[k,'id']
          curveSpline = splinefun(testCurvesToReturn[,1],testCurvesToReturn[,id+1],method='hyman')
          
          if(i>1) {
            #weight = (data[k,'time']-timesplit[i-1])/(timesplit[i] - timesplit[i-1])
            #weight = (predict(kmMod,data[k,'time'])-predict(kmMod,timesplit[i]))/(predict(kmMod,timesplit[i-1]) - predict(kmMod,timesplit[i]))
            evidence = data[k,]
            weight = (curveSpline(timesplit[i])/curveSpline(timesplit[i-1]))
          }else {
            #weight = data[k,'time']/timesplit[i]
            #weight = (predict(kmMod,data[k,'time'])-predict(kmMod,timesplit[i]))/(predict(kmMod,0) - predict(kmMod,timesplit[i]))
            evidence = data[k,]
            weight = curveSpline(timesplit[i])
          }
          if(weight<0){weight = 0}
          if(curveSpline(timesplit[i])<0){weight = 0}
          negativeWeight = 1-weight
          # if(is.na(data[k,'PREVTIMEPOINT'])) {
          #   weight = weight*0.1
          #   negativeWeight = negativeWeight*0.1
          # }else {
          #   weight = weight*0.5
          #   negativeWeight = negativeWeight*0.5
          # }
          #weight = 0
          
          index = data[k,p]
          index['TIMEPOINT'] = 1
          flatIndex = convertToFlatIndex(index,dim(newcpt))
          newcpt[flatIndex] = newcpt[flatIndex] + weight
          index['TIMEPOINT'] = 2
          flatIndex = convertToFlatIndex(index,dim(newcpt))
          newcpt[flatIndex] = newcpt[flatIndex] + negativeWeight
        }else {
          weight = 1
          index = data[k,p]
          flatIndex = convertToFlatIndex(index,dim(newcpt))
          newcpt[flatIndex] = newcpt[flatIndex] + negativeWeight
        }
        
        
      }
      
      normal = as.table(array(integer(prod(dim(newcpt)[-1])),dim=dim(newcpt)[-1]))
      for(j in 1:dim(newcpt)[1]) {
        multipleFlateIndex = convertToFlatIndexMultiple(c(j,rep(NA,length(p)-1)),dim(newcpt))
        normal = normal + newcpt[multipleFlateIndex]
      }
      
      normal = normal + iss
      newcpt[!is.na(newcpt)] = newcpt[!is.na(newcpt)] + (iss/dim(newcpt)[1])
      
      for(j in 1:dim(newcpt)[1]) {
        multipleFlateIndex = convertToFlatIndexMultiple(c(j,rep(NA,length(p)-1)),dim(newcpt))
        newcpt[multipleFlateIndex] = newcpt[multipleFlateIndex]/normal
      }
      
      sumOther = as.table(array(integer(prod(dim(newcpt)[-1])),dim=dim(newcpt)[-1]))
      for(j in 1:dim(newcpt)[1]) {
        if(j != dim(newcpt)[1]) {
          multipleFlateIndex = convertToFlatIndexMultiple(c(j,rep(NA,length(p)-1)),dim(newcpt))
          sumOther = sumOther + newcpt[multipleFlateIndex]
        }else if(j == dim(newcpt)[1]) {
          multipleFlateIndex = convertToFlatIndexMultiple(c(j,rep(NA,length(p)-1)),dim(newcpt))
          newcpt[multipleFlateIndex] = 1 - sumOther
        }
      }
      newcpt[is.na(newcpt)] = 0.5
      fitList[[i]][[node]] = newcpt
    }
  }
return(fitList)
}

convertToFlatIndex = function(index,dimension) {
  flatIndex = 0
  index = as.integer(index)
  for(digit in 1:length(index)) {
    flatIndex = flatIndex + (index[digit]-1)*prod(dimension[0:(digit-1)])
  }
  flatIndex = flatIndex+1
}

convertToFlatIndexMultiple = function(index,dimension) {
  flatIndexVec = c()
  index = as.integer(index)
  naDigit = is.na(index)
  index[naDigit] = 1
  naDigitDim = dimension[naDigit]
  
  for(n in 0:(prod(dimension[naDigit])-1)) {
    for(digit in 1:length(naDigitDim)) {
      index[naDigit][digit] = (floor(n/prod(naDigitDim[0:(digit-1)])) %% naDigitDim[digit])+1
    }
    flatIndex = convertToFlatIndex(index,dimension)
    flatIndexVec = append(flatIndexVec,flatIndex)
  }
  
  return(flatIndexVec)
}
