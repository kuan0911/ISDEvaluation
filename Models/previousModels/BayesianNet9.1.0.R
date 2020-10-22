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
source("Models/bayesianNetHelper.R")

BayesianNet = function(training, testing, BayesianC1 = 1, timePoints = 0,debug = FALSE){
  originalTraining = training
  originalTesting = testing
  
  C1 = NULL
  if(is.null(C1)){
    C1 = mtlr_cv(Surv(time,delta)~.,data=originalTraining, loss= "conc", C1_vec = c(0.001,0.01,0.1,1,10,100,1000)
                 , train_biases = F,train_uncensored = F)$best_C1
  }
  mod = mtlr(Surv(time,delta)~., data = originalTraining, C1=C1, train_biases = F, train_uncensored = F)
  
  
  queryMethod = 'naivebayes'
  prior = T
  weighted = T
  
  print(timePoints)
  
  numTimepoint = length(timePoints)
  
  allData = rbind(testing,training)
  
  cleaned = cleanDataType(nrow(testing),allData,discretize_level=10)
  testing = cleaned$test
  training = cleaned$train
  noCensorTraining = training[training$delta == 1,]

  dataList = createTimeSplit(training,timePoints,assumption = F,includeNa = F)
  dataListNA = createTimeSplit(training,timePoints,includeNa = T)
  print(length(dataList))
  
  mLess = 5
  quantileVals = seq(0,1,length.out = mLess+2)[-c(1)]
  #quantileVals = seq(0,1,length.out = m+2)[-1]
  timesplitLess = unname(quantile(training$time, quantileVals))
  timesplitLess = timesplitLess[!duplicated(timesplitLess)]
  dataListLess = createTimeSplit(training,timesplitLess,assumption = F,includeNa = T)

  numTimepoint = length(dataList)
  
  fitList <- vector("list", numTimepoint)
  dagList <- vector("list", numTimepoint)
  nbList <- vector("list", numTimepoint)
  childrenList <- vector("list", numTimepoint)
  print('learn starting graph')

  allmb = c()
  
  
  #whiteList = c('PREVTIMEPOINT','TIMEPOINT')
  
  iLess = 1
  for(i in 1:numTimepoint) {
    if(isTRUE(debug)){cat(i)
      cat(' ')}
    data = dataList[[i]]
    # if(i>1&i<numTimepoint) {
    #   data = rbind(dataList[[i-1]],dataList[[i]],dataList[[i+1]])
    # }else if(i==1){
    #   data = rbind(dataList[[i]],dataList[[i+1]])
    # }else if(i==numTimepoint) {
    #   data = rbind(dataList[[i]],dataList[[i-1]])
    # }
    
    rcount = nrow(data[data$PREVTIMEPOINT == 0,])
    #rcount = nrow(data)
    data = data[data$PREVTIMEPOINT == 0 & !is.na(data$PREVTIMEPOINT),]
    data$PREVTIMEPOINT = NULL
    data$time = NULL
    data$delta = NULL
    data$id = NULL
    #data[is.na(data$TIMEPOINT),'TIMEPOINT'] = 0
    
    if(timePoints[i]>timesplitLess[iLess]) {iLess = iLess+1}
    structureData = dataListLess[[iLess]]
    #structureData = dataList[[i]]
    structureData$PREVTIMEPOINT = NULL
    structureData$time = NULL
    structureData$delta = NULL
    structureData$id = NULL
    
    #structureData = structureData[complete.cases(structureData),]
    #structureData = rbind(structureData,data)
    #structureData = data
    
    #whitelist = rbind(c('LDH_SERUM','TIMEPOINT'),c('GRANULOCYTES','TIMEPOINT'))
    blacklist = blacklistFunction(colnames(data))
    #dag = hc(structureData,restart=50,blacklist=NULL,whitelist = NULL,start=NULL,score='bic')
    dag = structural.em(structureData, maximize = "hc",maximize.args = list(restart=40,blacklist=NULL,whitelist=NULL), fit = "bayes",fit.args=list(iss=20),impute='bayes-lw',return.all = T,start = NULL, max.iter = 15, debug = FALSE)
    
    # data = data[complete.cases(data),]
    # if(length(mb(dag$dag,'TIMEPOINT'))==0){p=c('TIMEPOINT','LDH_SERUM',mb(dag$dag,'TIMEPOINT'))}
    # else{p=c('TIMEPOINT',mb(dag$dag,'TIMEPOINT'))}
    # selectedData = data[,p]
    # selectedStructureData = structureData[,p]
    # dag = structural.em(selectedStructureData, maximize = "hc",maximize.args = list(restart=40,blacklist=NULL,whitelist=NULL), fit = "bayes",fit.args=list(iss=20),impute='bayes-lw',return.all = T,start = NULL, max.iter = 20, debug = FALSE)

    #tan = tree.bayes(selectedData, "TIMEPOINT")
    #emResult = structural.em(data, maximize = "hc",maximize.args = list(restart=20,blacklist=NULL,whitelist=NULL), fit = "bayes",fit.args=list(iss=10),impute='bayes-lw',return.all = T,start = fitList[[i]], max.iter = 20, debug = FALSE)
    #dag = emResult$dag
    #dag = cextend(dag, strict = TRUE, debug = FALSE)
    #dag = start
    #dagList[[i]] = dag
    #start = dag
    fit = bn.fit(dag$dag, data,method='bayes',iss=10)
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
  mList = c(5,10,20)
  priorFitList = NULL
  priorTimePoints = NULL
  for(iter in 1:3) {
    m = mList[iter]
    quantileVals = seq(0,1,length.out = m+2)[-c(1,m+2)]
    newTimePoints = unname(quantile(training$time, quantileVals))
    newTimePoints = newTimePoints[!duplicated(newTimePoints)]
    
    dataList = createTimeSplit(training,newTimePoints,assumption = F,includeNa = F)
    newFitList <- vector("list", numTimepoint)
    kLess = 1
    for(k in 1:length(newTimePoints)) {
      if(newTimePoints[k]>timePoints[kLess]){kLess+1}
      newFitList[k] = fitList[kLess]
    }
    
    if(prior) {
      newFitList = timepointPrior(dataList,newFitList,training,mod,newTimePoints,priorFitList,priorTimePoints)
    }
    
  
    if(weighted) {
      newFitList = weightedLearning(dataList,dataListNA,newFitList,newTimePoints,training,mod,priorFitList,priorTimePoints)
    }
    
    priorFitList = newFitList
    priorTimePoints = newTimePoints
    
  }

  #plot with bigger fontsize
  #plotDag(dagList[[1]])
  print('start predict')
  #prediction
  survivalFunctionTesting = predictFunction(newFitList,testing,newTimePoints,method = queryMethod)
  survivalFunctionTraining = predictFunction(newFitList,training,newTimePoints,method = queryMethod)
  
  testCurvesToReturn = survivalFunctionTesting
  timesAndCensTest = cbind.data.frame(time = originalTesting$time, delta = originalTesting$delta)
  timesAndCensTrain = cbind.data.frame(time = originalTraining$time, delta = originalTraining$delta)
  trainingCurvesToReturn = survivalFunctionTraining
  
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}

predictFunction <- function(fitList,testing,timePoints,method) {
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
    cat(i)
    cat(' ')
    threshold = 0.001
    
    
    for(j in 1:nrow(testing)) {
      #print(nrow(testing))
      #print(ncol(testCurvesToReturn))

      if(method == 'lw') {
        eviList = as.list(testing[j,])
        prob = cpquery(tempFit, event = (TIMEPOINT == 0), evidence = as.list(testing[j,]),method = 'lw')
      }else if(method == 'ls') {
        #evidence = testing[j,c('TIMEPOINT',children(dag$dag,'TIMEPOINT'))]
        #prob = predict(tempFit, evidence)
      }else if(method == 'predict') {
        # cpt = tempFit[["TIMEPOINT"]][['prob']]
        # p = colnames(as.data.frame(cpt))
        # p = p[p!="Freq"]
        # if(p=='Var1') {p=c('TIMEPOINT')}
        evidence = testing[j,]
        predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "parents", prob = TRUE)
        attr(predicted, "prob")
        prob = attr(predicted, "prob")[1]
      }else if(method == 'naivebayes') {
        evidence = testing[j,]
        cpt = tempFit[["TIMEPOINT"]][['prob']]
        p = colnames(as.data.frame(cpt))
        p = p[p!="Freq"]
        if(p=='Var1') {p=c('TIMEPOINT')}
        temp_evidence = evidence
        temp_evidence$TIMEPOINT = factor(0,levels = c('0','1'))
        ap1 = cpt[convertToFlatIndex(temp_evidence[p],dim(cpt))]
        temp_evidence$TIMEPOINT = factor(1,levels = c('0','1'))
        ap2 = cpt[convertToFlatIndex(temp_evidence[p],dim(cpt))]

        p0 = 0
        p1 = 0
        variables = children(fitList[[i]],'TIMEPOINT')
        for(v in variables) {
          e = as.integer(evidence[[v]])
          cpt = tempFit[[v]][['prob']]
          p = colnames(as.data.frame(cpt))
          p = p[p!="Freq"]
          temp_evidence = evidence
          temp_evidence$TIMEPOINT = factor(0,levels = c('0','1'))
          t1 = cpt[convertToFlatIndex(temp_evidence[p],dim(cpt))]
          temp_evidence$TIMEPOINT = factor(1,levels = c('0','1'))
          t2 = cpt[convertToFlatIndex(temp_evidence[p],dim(cpt))]
          
          if(t1>threshold) {rawp0=t1} else {rawp0=threshold}
          if(t2>threshold) {rawp1=t2} else {rawp1=threshold}
          p0 = p0 + log(rawp0)
          p1 = p1 + log(rawp1)
        }
        #print(t)
        p0 = p0 + log(ap1)
        p1 = p1 + log(ap2)
        p0_Prob = exp(p0)/(exp(p0)+exp(p1))
        
        prob = p0_Prob
      }else{
        print('No such query method')
      }

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
  #survivalFunction = cbind(time = c(0,timePoints), survivalFunction) 
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
  blackList = rbind(blackList, blackList3)
  
  return(blackList)
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


timepointPrior <- function(dataList,fitList,training,mod,timePoints,priorFitList,priorTimePoints) {
  print('Impose prior to TIMEPOINT')
  kmMod = prodlim(Surv(time,delta)~1, data = training)
  
  testCurvesToReturn = predict(mod)
  
  for(i in 1:length(dataList)) {
    data = dataList[[i]]
    #data = data[complete.cases(data),]
    data = data[data[,'PREVTIMEPOINT']==0,]
    #data = data[data[,'PREVTIMEPOINT']==0 | is.na(data[,'PREVTIMEPOINT']),]
    #data$PREVTIMEPOINT = NULL
    data$time = NULL
    data$delta = NULL
    #data$id = NULL
    dataComplete = data[complete.cases(data),]
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
    for(k in 1:length(timepointNA)){timepointNA[k]=0}
    for(k in 1:length(timepointNATotal)){timepointNATotal[k]=0}
    for(k in 1:nrow(data)) {
      if(is.na(data[k,'TIMEPOINT'])) {
      #if(T) {
        id = data[k,'id']
        curveSpline = splinefun(testCurvesToReturn[,1],testCurvesToReturn[,id+1],method='hyman')
        
        if(i>1) {
          #weight = (data[k,'time']-timePoints[i-1])/(timePoints[i] - timePoints[i-1])
          weight = curveSpline(timePoints[i])/curveSpline(timePoints[i-1])
          #weight = predict(kmMod,timePoints[i])/predict(kmMod,timePoints[i-1])
        }else {
          #weight = data[k,'time']/timePoints[i]
          weight = curveSpline(timePoints[i])
          #weight = predict(kmMod,timePoints[i])
        }
        if(weight<0 | is.na(weight)){weight = 0}
        if(curveSpline(timePoints[i])<0){weight = 0}
        deadWeight = 1 - weight
        # if(is.na(data[k,'PREVTIMEPOINT'])) {
        #   deadWeight = deadWeight*0.1
        # }else {
        #   deadWeight = deadWeight*0.5
        # }
        
        index = data[k,p[-1]]
        flatIndex = convertToFlatIndex(index,dim(timepointNA))
        timepointNA[flatIndex] = timepointNA[flatIndex] + deadWeight
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
    #       #weight = (data[k,'time']-timePoints[i-1])/(timePoints[i] - timePoints[i-1])
    #       weight = (predict(kmMod,data[k,'time'])-predict(kmMod,timePoints[i]))/(predict(kmMod,timePoints[i-1]) - predict(kmMod,timePoints[i]))
    #     }else {
    #       #weight = data[k,'time']/timePoints[i]
    #       weight = (predict(kmMod,data[k,'time'])-predict(kmMod,timePoints[i]))/(predict(kmMod,0) - predict(kmMod,timePoints[i]))
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
    if(!is.null(priorFitList)) {
      prior_i=1
      for(k in 1:length(priorTimePoints)) {
        if(priorTimePoints[k]>timePoints[i]) {
          prior_i = k
          break
        }
      }
      priorCpt = coef(priorFitList[[prior_i]]$TIMEPOINT)
      multipleFlateIndex = convertToFlatIndexMultiple(c(1,rep(NA,length(p)-1)),dim(newcpt))
      priorTimepointProb = priorCpt[multipleFlateIndex]
    }
    
    
    iss = 20
    iss_total = iss
    iss_dead = iss * (dead/(alive+dead+censored*0.5))
    
    timepointProb = timepoint1
    for(j in 1:length(timepoint1)) {
      if(!is.null(priorFitList)) {
        iss_dead = iss * (1-sqrt(priorTimepointProb[j]))
      }
      #timepoint1[j] = rbeta(1,(numTimepoint-1)*weight*15+timepoint1[j]*5,1*weight*15+timepoint2[j]*5)
      survivalRate = 1 - (timepoint2[j]+iss_dead)/(timepoint1[j]+timepoint2[j]+0.5*timepointNATotal[j]+iss_total)
      #survivalRate = 1- (timepointNA[j]+iss_dead)/(timepointNATotal[j]+iss_total)
      #survivalRate = 1-(iss_dead/iss_total)
      #print(sum(timepointNA)/sum(timepointNATotal))
      #print(iss_dead/iss_total)
      if(is.nan(survivalRate)) {
        survivalRate = 1-(iss_dead/iss_total)
        print('divided by zero')
      }
      prob = survivalRate
      timepointProb[j] = prob
    }
    
    multipleFlateIndex = convertToFlatIndexMultiple(c(1,rep(NA,length(p)-1)),dim(newcpt))
    newcpt[multipleFlateIndex] = timepointProb
    multipleFlateIndex = convertToFlatIndexMultiple(c(2,rep(NA,length(p)-1)),dim(newcpt))
    newcpt[multipleFlateIndex] = 1 - timepointProb
    
    fitList[[i]]$TIMEPOINT = newcpt
  }
  return(fitList)
}

weightedLearning <- function(dataList,dataListNA,fitList,timePoints,training,mod,priorFitList,priorTimePoints) {
  kmMod = prodlim(Surv(time,delta)~1, data = training)
  
  
  testCurvesToReturn = predict(mod)
  
  
  print('Weighted fitting')
  for(i in 1:length(dataList)) {
    for(node in children(fitList[[i]],'TIMEPOINT')) {
      data = dataList[[i]]
      #data = data[complete.cases(data),]
      data = data[data['PREVTIMEPOINT']==0,]
      #data = data[data['PREVTIMEPOINT']==0 | is.na(data[,'PREVTIMEPOINT']),]
      #data$PREVTIMEPOINT = NULL
      
      cpt = coef(fitList[[i]][[node]])
      p = colnames(as.data.frame(cpt))
      p = p[p!="Freq"]
      #iss = nrow(data)
      newcpt = cpt
      for(k in 1:length(newcpt)){newcpt[k]=0}
      
      for(k in 1:nrow(data)) {
        #if(T) {
        if(is.na(data[k,'TIMEPOINT'])) {
          
          
          id = data[k,'id']
          curveSpline = splinefun(testCurvesToReturn[,1],testCurvesToReturn[,id+1],method='hyman')
          
          if(i>1) {
            #weight = (data[k,'time']-timePoints[i-1])/(timePoints[i] - timePoints[i-1])
            #weight = (predict(kmMod,data[k,'time'])-predict(kmMod,timePoints[i]))/(predict(kmMod,timePoints[i-1]) - predict(kmMod,timePoints[i]))
            evidence = data[k,]
            weight = (curveSpline(timePoints[i])/curveSpline(timePoints[i-1]))
          }else {
            #weight = data[k,'time']/timePoints[i]
            #weight = (predict(kmMod,data[k,'time'])-predict(kmMod,timePoints[i]))/(predict(kmMod,0) - predict(kmMod,timePoints[i]))
            evidence = data[k,]
            weight = curveSpline(timePoints[i])
          }
          if(weight<0){weight = 0}
          if(curveSpline(timePoints[i])<0){weight = 0}
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
          prior = sum(as.integer(data[,node])==as.integer(index[node]))/length(data[,node])
          if(prior==0){
            print(prior)
            prior=0.1
            }
          # index['TIMEPOINT'] = 1
          # flatIndex = convertToFlatIndex(index,dim(newcpt))
          # newcpt[flatIndex] = newcpt[flatIndex] + (weight)
          # index['TIMEPOINT'] = 2
          # flatIndex = convertToFlatIndex(index,dim(newcpt))
          # newcpt[flatIndex] = newcpt[flatIndex] + negativeWeight
        }else {
          weight = 1
          index = data[k,p]
          flatIndex = convertToFlatIndex(index,dim(newcpt))
          newcpt[flatIndex] = newcpt[flatIndex] + weight
        }
      }
      
      if(!is.null(priorFitList)) {
        prior_i=1
        for(k in 1:length(priorTimePoints)) {
          if(priorTimePoints[k]>timePoints[i]) {
            prior_i = k
            break
          }
        }
        priorCpt = coef(priorFitList[[prior_i]][[node]])
        #multipleFlateIndex1 = convertToFlatIndexMultiple(c(1,rep(NA,length(p)-1)),dim(newcpt))
        #multipleFlateIndex2 = convertToFlatIndexMultiple(c(2,rep(NA,length(p)-1)),dim(newcpt))
        #priorSum = priorCpt[multipleFlateIndex1] + priorCpt[multipleFlateIndex2]
        #priorCpt[multipleFlateIndex1] = sqrt(priorCpt[multipleFlateIndex1])
        #priorCpt[multipleFlateIndex2] = sqrt(priorSum) - priorCpt[multipleFlateIndex1]
        
      }
      
      
      iss = 20
      
      newcpt[is.na(newcpt)] = 0
      if(!is.null(priorFitList)) {
        newcpt = newcpt + iss*priorCpt
      }else{
        newcpt = newcpt + (iss/dim(newcpt)[1])
      }
      
      #newcpt = sqrt(newcpt)
      
      normal = as.table(array(integer(prod(dim(newcpt)[-1])),dim=dim(newcpt)[-1]))
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
