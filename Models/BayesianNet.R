#------------------------------------------------------------
#Cumulative model with seperate network at each timepoint
#include all data in structure learning. Na if dead or censored
#mtlr parameter learning
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
source("Models/bayesianNetHelper.R")

BayesianNet = function(training, testing, BayesianC1 = 1, timePoints = 0,debug = FALSE){
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
  prior = T
  weighted = T
  cont = F
  fillAllMB = F
  
  variableList = variableTypes(originalTraining,10)
  
  print(timePoints)
  
  numTimepoint = length(timePoints)
  
  allData = rbind(testing,training)
  
  cleaned = cleanDataType(nrow(testing),allData,allDiscrete=T,discretize_level=6)
  testing = cleaned$test
  training = cleaned$train
  
  dataList = createTimeSplit(training,timePoints,assumption = F,includeNa = F)
  dataListNA = createTimeSplit(training,timePoints,includeNa = T)
  
  cleanedCont = cleanDataType(nrow(testing),allData,allDiscrete=F,discretize_level=10)
  trainingCont = cleanedCont$train
  dataListCont = createTimeSplit(trainingCont,timePoints,assumption = F,includeNa = F)
  
  print(length(dataList))
  
  #mLess = floor(length(timePoints)/2)
  mLess = 2
  quantileVals = seq(0,1,length.out = mLess+2)[-c(1)]
  #quantileVals = seq(0,1,length.out = m+2)[-1]
  timesplitLess = unname(quantile(training$time, quantileVals))
  timesplitLess = timesplitLess[!duplicated(timesplitLess)]
  dataListLess = createTimeSplit(training,timesplitLess,assumption = F,includeNa = T)
  
  numTimepoint = length(dataList)
  
  fitList <- vector("list", numTimepoint)
  fitListParents <- vector("list", numTimepoint)
  dagList <- vector("list", numTimepoint)
  nbList <- vector("list", numTimepoint)
  contFittingList <- vector("list", numTimepoint)
  print('learn starting graph')
  
  allmb = c()
  
  
  # structureData = dataListLess[[2]]
  # #structureData = dataList[[i]]
  # structureData$PREVTIMEPOINT = NULL
  # structureData$time = NULL
  # structureData$delta = NULL
  # structureData$id = NULL
  # 
  # dag = structural.em(structureData, maximize = "hc",maximize.args = list(restart=15,blacklist=NULL,whitelist=NULL), fit = "bayes",fit.args=list(iss=20),impute='bayes-lw',return.all = T,start = NULL, max.iter = 2, debug = FALSE)
  # # blackWhiteList = blacklistFunctionFixParents(colnames(structureData),parents(dag$dag,'TIMEPOINT'))
  
  #whiteList = c('PREVTIMEPOINT','TIMEPOINT')
  for(iter in 1:1) {
    iLess = 1
    structureData = dataListLess[[1]]
    #structureData = dataListLess[[2]]
    #structureData = dataList[[i]]
    structureData$PREVTIMEPOINT = NULL
    structureData$time = NULL
    structureData$delta = NULL
    structureData$id = NULL
    
    #blackList = blacklistFunction(names(variableList[variableList==T&names(variableList)!='time']))
    #blackList = blacklistFunction(colnames(structureData))
    start = structural.em(structureData, maximize = "hc",maximize.args = list(restart=25,blacklist=NULL,whitelist=NULL,score='bic'), fit = "bayes",fit.args=list(iss=20),impute='bayes-lw',return.all = T,start = NULL, max.iter = 2, debug = FALSE)
    
    #plotDag(start$dag)
    # print(parents(dag$dag,'TIMEPOINT'))
    # print(children(dag$dag,'TIMEPOINT'))
    #blackWhiteList = blacklistFunctionFixParents(colnames(structureData),parents(dag$dag,'TIMEPOINT'))
    
    for(i in 1:numTimepoint) {
      fitListParents[[i]] = start$dag
    }
    
    for(i in 1:numTimepoint) {
      if(isTRUE(debug)){cat(i)
        cat(' ')}
      data = dataListNA[[i]]
      
      # for(k in 1:nrow(data)) {
      #   if(!is.na(data[k,'PREVTIMEPOINT'])) {
      #     if(data[k,'PREVTIMEPOINT']==1) {
      #       data[k,'TIMEPOINT'] = NA
      #     }
      #   }
      # }

      # }
      # if(i>1&i<numTimepoint) {
      #   data = rbind(dataList[[i-1]],dataList[[i]],dataList[[i+1]])
      # }else if(i==1){
      #   data = rbind(dataList[[i]],dataList[[i+1]])
      # }else if(i==numTimepoint) {
      #   data = rbind(dataList[[i]],dataList[[i-1]])
      # }
      
      #rcount = nrow(data)
      data = data[data$PREVTIMEPOINT == 0,]
      data$PREVTIMEPOINT = NULL
      data$time = NULL
      data$delta = NULL
      data$id = NULL

      #dataDie = data[data$TIMEPOINT == 1,]
      #datadownsample = data[sample(nrow(data), nrow(dataDie)), ]
      #datadownsample = rbind(datadownsample,dataDie)
      
      if(timePoints[i]>timesplitLess[iLess]) {
        iLess = iLess+1
        structureData = dataListLess[[iLess]]
        #structureData = dataList[[i]]
        structureData$PREVTIMEPOINT = NULL
        structureData$time = NULL
        structureData$delta = NULL
        structureData$id = NULL
        
        #blackList = blacklistFunction(names(variableList[variableList==T&names(variableList)!='time']))
        #start = structural.em(structureData, maximize = "hc",maximize.args = list(restart=20,blacklist=NULL,whitelist=NULL), fit = "bayes",fit.args=list(iss=20),impute='bayes-lw',return.all = T,start = NULL, max.iter = 1, debug = FALSE)
        #plotDag(start$dag)
        #print(parents(dag$dag,'TIMEPOINT'))
        #print(children(dag$dag,'TIMEPOINT'))
      }
      
      # structureData = dataList[[i]]
      # structureData$PREVTIMEPOINT = NULL
      # structureData$time = NULL
      # structureData$delta = NULL
      # structureData$id = NULL
      dataComplete = data[complete.cases(data),]
      #blackList = blacklistFunction(names(variableList[variableList==T&names(variableList)!='time']))
      #if(i==1) {dag = structural.em(data, maximize = "hc",maximize.args = list(restart=100,blacklist=NULL,whitelist=NULL), fit = "bayes",fit.args=list(iss=2),impute='bayes-lw',return.all = T,start = start$dag, max.iter = 1, debug = FALSE)}
      #else {dag = structural.em(data, maximize = "hc",maximize.args = list(restart=10,blacklist=NULL,whitelist=NULL), fit = "bayes",fit.args=list(iss=2),impute='bayes-lw',return.all = T,start = start$dag, max.iter = 1, debug = FALSE)}
      #dag = structural.em(data, maximize = "tabu",maximize.args = list(tabu=200,blacklist=NULL,whitelist=NULL), fit = "bayes",fit.args=list(iss=2),impute='bayes-lw',return.all = T,start = start$dag, max.iter = 1, debug = FALSE)
      dag = hc(dataComplete,restart=100,start = start$dag)
      #dag = start
      # if(length(mb(dag$dag,'TIMEPOINT'))<length(mb(start$dag,'TIMEPOINT'))) {
      #   dag = start
      # }
      #start = dag
      #dag = start
      
      print(parents(dag,'TIMEPOINT'))
      print(children(dag,'TIMEPOINT'))
      #whitelist = rbind(c('LDH_SERUM','TIMEPOINT'),c('GRANULOCYTES','TIMEPOINT'))
      #blacklist = blacklistFunction(colnames(data)[variableList])
      #dag = hc(structureData,restart=50,blacklist=NULL,whitelist = NULL,start=NULL)
      
      fit = bn.fit(dag, data, method='bayes',iss=1)
      fitParents = bn.fit(fitListParents[[i]], data, method='bayes',iss=10)
      fitListParents[[i]] = fitParents
      #fit = emResult$fitted
      #fit = bn.fit(dag, data, method='bayes')
      fitList[[i]] <- fit
      #data = dataListNA[[i]]
      nbdata=data
      #nbdata$PREVTIMEPOINT = NULL
      nb_model = naiveBayes(as.factor(TIMEPOINT) ~., nbdata, laplace = 20)
      nbList[[i]] <- nb_model
      #print(parents(fitList[[i]],'TIMEPOINT'))
      allmb = unique(c(allmb,mb(fitList[[i]],'TIMEPOINT')))
      #plotDag(dag$dag)
    }
    
    prevFitList=NULL
    bnCurveList=NULL
    if(prior) {
      fitList = timepointPrior(dataListNA,fitList,prevFitList,originalTraining,mod,timePoints,kmMod,bnCurveList)
    }
    if(weighted) {
      fitList = weightedLearning(dataListNA,fitList,prevFitList,timePoints,originalTraining,mod,kmMod,bnCurveList)
    }
    prevFitList = fitList
    
    # bnCurveList = vector("list", nrow(training))
    # for(k in 1:nrow(training)){
    #   bnCurveList[[k]] = BnExactInferenceCumulated(prevFitList,training[k,])
    # }
    
    # for(para_iter in 1:5){
    #   if(prior) {
    #     fitList = timepointPrior(dataListNA,fitList,prevFitList,originalTraining,mod,timePoints,kmMod,bnCurveList)
    #   }
    #   if(weighted) {
    #     fitList = weightedLearning(dataListNA,fitList,prevFitList,timePoints,originalTraining,mod,kmMod,bnCurveList)
    #   }
    #   prevFitList = fitList
    # }
    
    
    if(cont) {
      contFittingList = contFitting(dataListCont,fitList,variableList,timePoints)
    }
    
  }
  
  #plot with bigger fontsize
  #plotDag(dagList[[1]])
  print('start predict')
  #prediction
  survivalFunctionTesting = predictFunction(fitList,fitListParents,nbList,testing,originalTesting,timePoints,queryMethod,dataList,contFittingList,cont,variableList,kmMod)
  survivalFunctionTraining = predictFunction(fitList,fitListParents,nbList,training,originalTraining,timePoints,queryMethod,dataList,contFittingList,cont,variableList,kmMod)
  
  testCurvesToReturn = survivalFunctionTesting
  timesAndCensTest = cbind.data.frame(time = originalTesting$time, delta = originalTesting$delta)
  timesAndCensTrain = cbind.data.frame(time = originalTraining$time, delta = originalTraining$delta)
  trainingCurvesToReturn = survivalFunctionTraining
  
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}

predictFunction <- function(fitList,fitListParents,nbList,testing,originalTesting,timePoints,queryMethod,dataList,contFittingList,cont,variableList,kmMod) {
  testing$time = NULL
  testing$delta = NULL
  #testing$PREVTIMEPOINT = factor(integer(nrow(testing)),levels = c('0','1'))
  numTimepoint = length(timePoints)
  numReturnNA = 0
  numNotDecreasing = 0
  threshold = 0.000001
  
  #first value in cpt
  previousTimepointProb = rep(1,nrow(testing))
  previousProb = rep(1-(1/numTimepoint),nrow(testing))
  previousProb2 = rep(1-(1/numTimepoint),nrow(testing))
  previousProb3 = rep(1-(1/numTimepoint),nrow(testing))
  
  survivalFunction <- data.frame(matrix(ncol = nrow(testing), nrow = numTimepoint))
  for(i in 1:numTimepoint) {
    tempFit = fitList[[i]]

    tempFitParents = fitListParents[[i]]
    tempNb = nbList[[i]]
    cat(i)
    cat(' ')
    threshold = 0.001
    
    for(j in 1:nrow(testing)) {
      #print(nrow(testing))
      #print(ncol(testCurvesToReturn))
      
      if(queryMethod == 'lw') {
        eviList = as.list(testing[j,])
        prob = cpquery(tempFit, event = (TIMEPOINT == 0), evidence = as.list(testing[j,]),method = 'lw')
      }else if(queryMethod == 'ls') {
        #evidence = testing[j,c('TIMEPOINT',children(dag$dag,'TIMEPOINT'))]
        #prob = predict(tempFit, evidence)
        #iss_total = (alive+dead+censored*0.5)*0.1
        #iss_dead = dead*0.1
        #prob = 1-(iss_dead/iss_total)
      }else if(queryMethod == 'predict') {
        # cpt = tempFit[["TIMEPOINT"]][['prob']]
        # p = colnames(as.data.frame(cpt))
        # p = p[p!="Freq"]
        # if(p=='Var1') {p=c('TIMEPOINT')}
        evidence = testing[j,]
        predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "bayes-lw", prob = TRUE)
        attr(predicted, "prob")
        prob = attr(predicted, "prob")[1]
      }else if(queryMethod == 'exact') {
        evidence = testing[j,]
        if(i>1){kmprob = predict(kmMod,timePoints[i])/predict(kmMod,timePoints[i-1])}else{kmprob = predict(kmMod,timePoints[i])}
        prob = BnExactInference(tempFit,evidence,kmprob=kmprob)
      }else if(queryMethod == 'naivebayes') {
        evidence = testing[j,]
        cpt = tempFit[["TIMEPOINT"]][['prob']]
        p = colnames(as.data.frame(cpt))
        p = p[p!="Freq"]
        if(p=='Var1') {p=c('TIMEPOINT')}
        temp_evidence = evidence
        temp_evidence$TIMEPOINT = factor(0,levels = c('0','1'))
        ap1 = cpt[convertToFlatIndex(temp_evidence[p],dim(cpt))]
        #temp_evidence$TIMEPOINT = factor(1,levels = c('0','1'))
        #ap2 = cpt[convertToFlatIndex(temp_evidence[p],dim(cpt))]
        #ap1 = ap1*previousTimepointProb[j]
        ap2 = 1-ap1
        
        p0 = 0
        p1 = 0
        variables = children(tempFit,'TIMEPOINT')
        for(v in variables) {
          if(cont==T & variableList[v]){
            par = contFittingList[[i]][[v]]
            t1 = dnorm(originalTesting[j,v], mean = par['mean0'], sd = par['sd0'])
            t2 = dnorm(originalTesting[j,v], mean = par['mean1'], sd = par['sd1'])
          }else {
            cpt = tempFit[[v]][['prob']]
            p = colnames(as.data.frame(cpt))
            p = p[p!="Freq"]
            temp_evidence = evidence
            temp_evidence$TIMEPOINT = factor(0,levels = c('0','1'))
            t1 = cpt[convertToFlatIndex(temp_evidence[p],dim(cpt))]
            temp_evidence$TIMEPOINT = factor(1,levels = c('0','1'))
            t2 = cpt[convertToFlatIndex(temp_evidence[p],dim(cpt))]
          }
          
          
          if(t1>threshold) {rawp0=t1} else {rawp0=threshold}
          if(t2>threshold) {rawp1=t2} else {rawp1=threshold}
          p0 = p0 + log(rawp0)
          p1 = p1 + log(rawp1)
          
        }
        #print(t)
        #p0 = p0
        #p1 = p1
        p0_Prob = exp(p0 + log(ap1))/(exp(p0 + log(ap1))+exp(p1 + log(ap2)))
        prob = p0_Prob
      }else{
        print('No such query method')
      }
      
      previousProb3[j] = previousProb2[j]
      previousProb2[j] = previousProb[j]
      previousProb[j] = prob
      # if(j==1) {prob = previousProb[j]}
      # else if(j==2) {prob = (previousProb3[j]+previousProb2[j])/2}
      # else {prob = (previousProb3[j]+previousProb2[j]+previousProb[j])/3}

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
  #survivalFunction = cbind(time = c(0,timePoints), survivalFunction) 
  return(survivalFunction)
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
  blackList = rbind(blackList, blackList3)
  
  return(blackList2)
}

blacklistFunctionFixParents <- function(nodeNames,parents){
  nodeNames = nodeNames[nodeNames!="TIMEPOINT"&nodeNames!="PREVTIMEPOINT"]
  notParents = setdiff(nodeNames,parents)
  blackList = matrix(ncol = 2, nrow = length(notParents))
  for(i in 1:length(notParents)) {
    blackList[i,] = c(notParents[i],"TIMEPOINT")
  }
  if(length(parents)>0) {
    whiteList = matrix(ncol = 2, nrow = length(parents))
    for(i in 1:length(parents)) {
      whiteList[i,] = c(parents[i],"TIMEPOINT")
    }
  }else {
    whiteList = NULL
  }
  
  
  return(list(blackList=blackList,whiteList=whiteList))
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


timepointPrior <- function(dataList,fitList,prevFitList,originalTraining,mod,timePoints,kmMod,bnCurveList) {
  print('Impose prior')
  #testCurvesToReturn = predict(mod)
  
  for(i in 1:length(dataList)) {
    data = dataList[[i]]
    #data = data[complete.cases(data),]
    dataiss = data[!is.na(data[,'PREVTIMEPOINT']) & data[,'PREVTIMEPOINT']==0,]
    #data = data[data[,'PREVTIMEPOINT']==0 | is.na(data[,'PREVTIMEPOINT']),]
    
    #data$id = NULL
    data$time = NULL
    data$delta = NULL
    dataComplete = data[complete.cases(data),]
    dataComplete = dataComplete[dataComplete[,'PREVTIMEPOINT']==0,]
    #data$PREVTIMEPOINT = NULL
    dataComplete$PREVTIMEPOINT = NULL
    cpt = coef(fitList[[i]]$TIMEPOINT)
    p = colnames(as.data.frame(cpt))
    p = p[p!="Freq"]
    if(p=='Var1') {p=c('TIMEPOINT')}
    newcpt = table(dataComplete[,p])
    
    multipleFlateIndex = convertToFlatIndexMultiple(c(1,rep(NA,length(p)-1)),dim(newcpt))
    timepoint1 = newcpt[multipleFlateIndex]
    multipleFlateIndex = convertToFlatIndexMultiple(c(2,rep(NA,length(p)-1)),dim(newcpt))
    timepoint2 = newcpt[multipleFlateIndex]
    
    timepointNA = timepoint1
    timepointNATotal = timepoint1
    timepointNACount = timepoint1
    timepointNANACount = timepoint1
    for(k in 1:length(timepointNA)){timepointNA[k]=0}
    for(k in 1:length(timepointNATotal)){timepointNATotal[k]=0}
    for(k in 1:length(timepointNACount)){timepointNACount[k]=0}
    for(k in 1:length(timepointNANACount)){timepointNANACount[k]=0}
    for(k in 1:nrow(data)) {
      if(F) {
      #if(is.na(data[k,'TIMEPOINT'])) {
      #if(is.na(data[k,'TIMEPOINT'])&is.na(data[k,'PREVTIMEPOINT'])) {
        id = data[k,'id']
        curveSpline = splinefun(testCurvesToReturn[,1],testCurvesToReturn[,id+1],method='hyman')
        if(!is.null(prevFitList)){bnCurve = bnCurveList[[id]]}
        
        if(i>1) {
          if(is.null(prevFitList)) {
            weight = curveSpline(timePoints[i-1])-curveSpline(timePoints[i])
            prior = curveSpline(timePoints[i-1])
          }else {
            weight = bnCurve[i-1]-bnCurve[i]
            prior = bnCurve[i-1]
          }
          
        }else {
          if(is.null(prevFitList)) {
            weight = 1-curveSpline(timePoints[i])
            prior = 1
          }else {
            weight = 1 - bnCurve[i]
            prior = 1
          }
          
        }
        if(weight<0 | is.na(weight)){weight = 0}
        if(curveSpline(timePoints[i])<0){weight = 0}
        deadWeight = 1 - weight
        
        index = data[k,p[-1]]
        if(length(index)==0) {
          timepointNA = timepointNA + weight
          timepointNATotal = timepointNATotal+prior
        }else {
          flatIndex = convertToFlatIndex(index,dim(timepointNA))
          timepointNA[flatIndex] = timepointNA[flatIndex] + weight
          flatIndex = convertToFlatIndex(index,dim(timepointNATotal))
          timepointNATotal[flatIndex] = timepointNATotal[flatIndex] + prior
        }
      }
      
      if(is.na(data[k,'TIMEPOINT'])&!is.na(data[k,'PREVTIMEPOINT'])&data[k,'PREVTIMEPOINT']==0) {
        index = data[k,p[-1]]
        if(length(index)==0) {
          timepointNACount = timepointNACount + 1
        }else {
          flatIndex = convertToFlatIndex(index,dim(timepointNACount))
          timepointNACount[flatIndex] = timepointNACount[flatIndex] + 1
        }
      }
      
      if(is.na(data[k,'TIMEPOINT'])&is.na(data[k,'PREVTIMEPOINT'])) {
        index = data[k,p[-1]]
        if(length(index)==0) {
          timepointNANACount = timepointNANACount + 1
        }else {
          flatIndex = convertToFlatIndex(index,dim(timepointNANACount))
          timepointNANACount[flatIndex] = timepointNANACount[flatIndex] + 1
        }
      }
      
    } 
    
    alive = sum(dataiss$TIMEPOINT==0,na.rm=T)
    dead = sum(dataiss$TIMEPOINT==1,na.rm=T)
    censored = sum(is.na(dataiss$TIMEPOINT))
    
    iss = 1
    iss_dead = iss * (dead/(alive+dead+censored*0.5))

    timepointProb = timepoint1
    
    for(k in 1:length(timepointProb)){timepointProb[k]=0}
    
    for(j in 1:length(timepointProb)) {
      #timepoint1[j] = rbeta(1,(numTimepoint-1)*weight*15+timepoint1[j]*5,1*weight*15+timepoint2[j]*5)
      #survivalRate = 1 - (timepoint2[j]+iss_dead)/(timepoint1[j]+timepoint2[j]+0.5*timepointNACount[j]+iss_total)
      #survivalRate = (timepointNA[j]+iss_total-iss_dead)/(timepointNATotal[j]+iss_total)
      #iss = (timepoint1[j]+timepoint2[j]+timepointNACount[j])
      # for(randomflip in 1:timepoint1[j]+timepoint2[j]) {
      #   if(runif(1, 0.0, 1.0)<0.5) {
      #     if(runif(1, 0.0, 1.0)<1-(iss_dead/iss)) {
      #       timepoint1[j] = timepoint1[j]+1
      #     }else{
      #       timepoint2[j] = timepoint2[j]+1
      #     }
      #   }
      # }
      #iss = timepointNANACount[j]
      iss=2
      if(iss<5) {iss=5}
      iss_dead = iss * (dead/(alive+dead+censored*0.5))
      # if(iss_dead<2){
      #   iss_dead = 2
      #   iss = (iss_dead*(alive+dead+censored*0.5))/dead
      # }
     
      # iss_dead = 3
      # iss = (iss_dead*(alive+dead+censored*0.5))/dead
      if(T) {
        if(length(timepoint2)>1) {
          survivalRate = 1 - (timepoint2[j]+iss_dead)/(timepoint1[j]+timepoint2[j]+0.5*timepointNACount[j]+iss)
        }else {
          survivalRate = 1 - (timepoint2+iss_dead)/(timepoint1+timepoint2+0.5*timepointNACount+iss)
        }
      }else{
        if(length(timepoint2)>1) {
          survivalRate = 1 - (timepointNA[j]+timepoint2[j]+iss_dead)/(timepoint1[j]+timepoint2[j]+timepointNATotal[j]+iss)
        }else {
          survivalRate = 1 - (timepointNA+timepoint2+iss_dead)/(timepoint1+timepoint2+timepointNATotal+iss)
        }
      }
      #survivalRate = 1-(iss_dead/iss_total)
      #survivalRate = 1 - sum(timepoint2)/(sum(timepoint1)+sum(timepoint2)+0.5*sum(timepointNATotal))
      #print(sum(timepointNA)/sum(timepointNATotal))
      #print(iss_dead/iss_total)
      if(is.nan(survivalRate)) {
        survivalRate = 1-(iss_dead/iss)
        print('divided by zero')
      }

      timepointProb[j] = survivalRate
    }
    
    # for(j in 1:length(timepointProb)) {
    #   if((timepoint1[j]+timepoint2[j]+timepointNACount[j])<20) {
    #     survivalRate = 1-(iss_dead/iss)
    #   }
    # }

    
    multipleFlateIndex = convertToFlatIndexMultiple(c(1,rep(NA,length(p)-1)),dim(newcpt))
    newcpt[multipleFlateIndex] = timepointProb
    multipleFlateIndex = convertToFlatIndexMultiple(c(2,rep(NA,length(p)-1)),dim(newcpt))
    newcpt[multipleFlateIndex] = 1 - timepointProb
    
    fitList[[i]]$TIMEPOINT = newcpt

  }
  return(fitList)
}

weightedLearning <- function(dataList,fitList,prevFitList,timePoints,originalTraining,mod,kmMod,bnCurveList) {
  
  #testCurvesToReturn = predict(mod)
  
  print('Weighted fitting')
  for(i in 1:length(dataList)) {
    for(node in children(fitList[[i]],'TIMEPOINT')) {
      data = dataList[[i]]
      #data = data[complete.cases(data),]
      #data = data[data['PREVTIMEPOINT']==0,]
      dataComplete = data[complete.cases(data),]
      dataComplete = dataComplete[dataComplete[,'PREVTIMEPOINT']==0,]
      #data = data[data['PREVTIMEPOINT']==0 | is.na(data[,'PREVTIMEPOINT']),]
      #data$PREVTIMEPOINT = NULL
      dataComplete$PREVTIMEPOINT = NULL
      
      
      cpt = coef(fitList[[i]][[node]])
      p = colnames(as.data.frame(cpt))
      p = p[p!="Freq"]
      #iss = nrow(data)
      newcpt = cpt
      for(k in 1:length(newcpt)){newcpt[k]=0}
      
      for(k in 1:nrow(data)) {
        if(F) {
        #if(is.na(data[k,'TIMEPOINT'])) {
          
          id = data[k,'id']
          curveSpline = splinefun(testCurvesToReturn[,1],testCurvesToReturn[,id+1],method='hyman')
          if(!is.null(prevFitList)){bnCurve = bnCurveList[[id]]}
          
          if(i>1) {
            if(is.null(prevFitList)) {
              weight = curveSpline(timePoints[i])
              negativeWeight = curveSpline(timePoints[i-1])-curveSpline(timePoints[i])
            }else {
              weight = bnCurve[i]
              negativeWeight = bnCurve[i-1]-bnCurve[i]
            }
          }else {
            if(is.null(prevFitList)) {
              weight = curveSpline(timePoints[i])
              negativeWeight = 1-curveSpline(timePoints[i])
            }else {
              weight = bnCurve[i]
              negativeWeight = 1-bnCurve[i]
            }
          }
          if(weight<0){weight = 0}
          if(curveSpline(timePoints[i])<0){weight = 0}
          #negativeWeight = 1-weight
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
          #}else if(!is.na(data[k,'TIMEPOINT'])) {
        }else if(F) {
          weight = 1
          index = data[k,p]
          flatIndex = convertToFlatIndex(index,dim(newcpt))
          newcpt[flatIndex] = newcpt[flatIndex] + weight
        }
      }
      
      for(k in 1:nrow(dataComplete)) {
        index = dataComplete[k,p]
        flatIndex = convertToFlatIndex(index,dim(newcpt))
        newcpt[flatIndex] = newcpt[flatIndex] + 1
      }
      
      # for(k in 1:nrow(data)) {
      #   if(is.na(data[k,'TIMEPOINT'])&!is.na(data[k,'PREVTIMEPOINT'])) {
      #     index = data[k,p]
      #     index['TIMEPOINT'] = 1
      #     flatIndex = convertToFlatIndex(index,dim(newcpt))
      #     newcpt[flatIndex] = newcpt[flatIndex] + 0.5
      #     index['TIMEPOINT'] = 2
      #     flatIndex = convertToFlatIndex(index,dim(newcpt))
      #     newcpt[flatIndex] = newcpt[flatIndex] + 0.5
      #   }
      # }
      
      newcpt[is.na(newcpt)] = 0
      iss = 5
      newcpt = newcpt + iss
      
      # if(i>6){
      #   for(k in 1:length(newcpt)){newcpt[k]=5}
      # }
      
      normal = as.table(array(integer(prod(dim(newcpt)[-1])),dim=dim(newcpt)[-1]))
      for(z in 1:length(normal)){normal[z]=0}
      # for(j in 1:dim(newcpt)[1]) {
      #   multipleFlateIndex = convertToFlatIndexMultiple(c(j,rep(NA,length(p)-1)),dim(newcpt))
      #   normal = normal + newcpt[multipleFlateIndex]
      # }
      # 
      # normal = normal/2
      # 
      # for(j in 1:dim(newcpt)[1]) {
      #   multipleFlateIndex = convertToFlatIndexMultiple(c(j,rep(NA,length(p)-1)),dim(newcpt))
      #   newcpt[multipleFlateIndex] = newcpt[multipleFlateIndex]+normal/dim(newcpt)[1]
      # }
      
      for(z in 1:length(normal)){normal[z]=0}
      for(j in 1:dim(newcpt)[1]) {
        multipleFlateIndex = convertToFlatIndexMultiple(c(j,rep(NA,length(p)-1)),dim(newcpt))
        normal = normal + newcpt[multipleFlateIndex]
      }
      
      for(j in 1:dim(newcpt)[1]) {
        multipleFlateIndex = convertToFlatIndexMultiple(c(j,rep(NA,length(p)-1)),dim(newcpt))
        newcpt[multipleFlateIndex] = newcpt[multipleFlateIndex]/normal
      }
      
      sumOther = as.table(array(integer(prod(dim(newcpt)[-1])),dim=dim(newcpt)[-1]))
      for(z in 1:length(sumOther)){sumOther[z]=0}
      for(j in 1:dim(newcpt)[1]) {
        if(j != dim(newcpt)[1]) {
          multipleFlateIndex = convertToFlatIndexMultiple(c(j,rep(NA,length(p)-1)),dim(newcpt))
          sumOther = sumOther + newcpt[multipleFlateIndex]
        }else if(j == dim(newcpt)[1]) {
          multipleFlateIndex = convertToFlatIndexMultiple(c(j,rep(NA,length(p)-1)),dim(newcpt))
          newcpt[multipleFlateIndex] = 1 - sumOther
        }
      }
      
      fitList[[i]][[node]] = newcpt
    }
  }
  return(fitList)
}

contFitting = function(dataListCont,fitList,variableList,timePoints) {
  print('fitting continuous data')
  contFittingList = vector("list", length(fitList))
  for(fit in 1:length(fitList)) {
    data = dataListCont[[fit]]
    data = data[data$PREVTIMEPOINT==0,]
    children = children(fitList[[fit]],'TIMEPOINT')
    children = children[variableList[children]]
    parList = c()
    for(child in children) {
      par = c()
      data0 = data[data$TIMEPOINT==0,][[child]]
      data1 = data[data$TIMEPOINT==1,][[child]]
      #bar0 = c(quantile(data0,0.01,na.rm=T),quantile(data0,0.99,na.rm=T))
      #bar1 = c(quantile(data1,0.01,na.rm=T),quantile(data1,0.99,na.rm=T))
      #data0 = data0[data0>bar0[1]&data0<bar0[2]]
      #data1 = data1[data1>bar1[1]&data1<bar1[2]]
      par['mean0'] = mean(data0,na.rm=T)
      par['mean1'] = mean(data1,na.rm=T)
      par['sd0'] = sd(data0,na.rm=T)
      par['sd1'] = sd(data1,na.rm=T)
      #print(sd(data0,na.rm=T))
      #print(sd(data1,na.rm=T))
      parList[[child]] = par
    }
    contFittingList[[fit]] = parList
  }
  return(contFittingList)
}



