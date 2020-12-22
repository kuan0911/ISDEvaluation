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
#library(e1071)
library(prodlim)
library(MTLR)
library(Rcpp)
source("Models/cleanDataForBayesNet.R")
source("Models/bayesianNetHelper.R")
source("Models/changesource.R")

BayesianNet = function(training, testing, timePoints = 0,debug = FALSE){
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
  cont = F
  fillAllMB = F
  
  variableList = variableTypes(originalTraining,10)
  
  print(timePoints)
  
  numTimepoint = length(timePoints)
  
  allData = rbind(testing,training)
  
  cleaned = cleanDataType(nrow(testing),allData,allDiscrete=T,discretize_level=3)
  testing = cleaned$test
  training = cleaned$train
  
  trainingList = vector("list", numTimepoint)
  testingList = vector("list", numTimepoint)
  for(i in 1:numTimepoint) {
    levelnum = numTimepoint - i
    if(levelnum>10) {levelnum=10}
    else if(levelnum<3) {levelnum=3}
    cleaned = cleanDataType(nrow(testing),allData,allDiscrete=T,discretize_level=levelnum)
    testingList[[i]] = cleaned$test
    trainingList[[i]] = cleaned$train
  }
  
  cleanedLevels = cleanDataType(nrow(testing),allData,allDiscrete=T,discretize_level=3)
  testingLevels = cleanedLevels$test
  trainingLevels = cleanedLevels$train

  dataList = createTimeSplit(training,timePoints,assumption = F,includeNa = F)
  dataListNA = createTimeSplit(training,timePoints,includeNa = T)

  dataListNAadapt = vector("list", numTimepoint)
  
  for(i in 1:numTimepoint) {
    dataListNAList = createTimeSplit(trainingList[[i]],timePoints,includeNa = T)
    dataListNAadapt[[i]] = dataListNAList[[i]]
  }
  
  cleanedCont = cleanDataType(nrow(testing),allData,allDiscrete=F,discretize_level=10)
  trainingCont = cleanedCont$train
  dataListCont = createTimeSplit(trainingCont,timePoints,assumption = F,includeNa = F)
  
  print(length(dataList))
  
  #mLess = floor(length(timePoints)/2)
  mLess = 3
  quantileVals = seq(0,1,length.out = mLess+2)[-c(1)]
  #quantileVals = seq(0,1,length.out = m+2)[-1]
  timesplitLess = unname(quantile(training[training$delta==1,]$time, quantileVals))
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
  
  prevLikelihood = 0
  
  for(iter in 1:1) {
    if(iter==1) {
      iLess = 1
      structureData = dataListLess[[2]]
  
      #structureData = do.call("rbind", dataListLess)
      #structureData = structureData[structureData$PREVTIMEPOINT==0,]
      structureData[,c('PREVTIMEPOINT','time','delta','id')] = NULL
      allstart = structural.em(structureData, maximize = "hc",maximize.args = list(restart=500,blacklist=NULL,whitelist=NULL,score='bic'), fit = "bayes",fit.args=list(iss=2),impute='bayes-lw',return.all = T,start = NULL, max.iter = 2, debug = FALSE)
      prevFitS = NULL
      prevDag = NULL
      prevFit = NULL
      
      allstart$fitted = timepointFixer(allstart$fitted,structureData)
      #plotDag(start$dag)
      print(mb(allstart$dag,'TIMEPOINT'))
    }
    for(i in 1:numTimepoint) {
      if(isTRUE(debug)){cat(i)
        cat(' ')}

      data = dataListNAadapt[[i]]
      
      datapara2 = data
      datapara2 = datapara2[datapara2$PREVTIMEPOINT == 0|is.na(datapara2$PREVTIMEPOINT),]

      data[,c('time','delta','id')] = NULL
      
      # dataStart = data
      # dataStart$PREVTIMEPOINT = NULL
      
      datapara = data
      datapara = datapara[!is.na(datapara$PREVTIMEPOINT),]
      datapara = datapara[datapara$PREVTIMEPOINT == 0,]
      datapara$PREVTIMEPOINT = NULL
      
      data$PREVTIMEPOINT = NULL
      
      #dataComplete = dataStart[complete.cases(dataStart),]
      #internaltest = datapara[complete.cases(datapara),]
      
      #imputedRes = imputeFunctionSuvival(fit,prevFitS,datapara2)
      #imputed = imputedRes$data
      #weight = imputedRes$weight
      #dag = hc(internaltest,restart=100,start=allstart$dag,score='bic')
      # #dag = structural.em(datapara, maximize = "hc",maximize.args = list(restart=100,score='bic'), fit = "bayes",impute='bayes-lw',return.all = T,start = allstart$dag, max.iter = 5, debug = FALSE)
      # #dag = dag$dag
      # fit = bn.fit(dag, datapara, method='bayes')
      # fit = timepointFixer(fit,datapara)
      if(iter>1) {
        fit = fitList[[i]]
        imputeRes = weightedImpute(fitList,fit,dataListNAadapt,datapara2,i)
      }else {
        internaltest = datapara[complete.cases(datapara),]
        dag = hc(internaltest,restart=100,start=allstart$dag,score='bic')
        fit = bn.fit(dag, datapara, method='bayes')
        fit = timepointFixer(fit,datapara)
        #imputeRes = weightedImputeKM(kmMod,datapara2,timePoints,i)
      }
      #imputed = imputeRes$data
      #weight = imputeRes$weight
      # prevLogLike = 0
      # maxfit = fit
      #dag = hc(internaltest,start=NULL,score="custom",fun = customScoreFunction, args=list(weight=weight,prevFit=NULL,nextFit=NULL))
      #fit = bn.fit(dag, datapara, method='bayes')
      #fit = timepointFixer(fit,imputed,weight)
      #fit = timepointFixer(fit,datapara)
      #fit = childrenFixer(fit,imputed,weight)
      
        
      #prevDag = dag
      fitList[[i]] <- fit
      print(mb(fitList[[i]],'TIMEPOINT'))
      allmb = unique(c(allmb,mb(fitList[[i]],'TIMEPOINT')))
    }
    
    # imputenum = 0
    # imputedDataListNa = dataListNA
    # for(dataListIter in 1:length(fitList)) {
    #   for(k in 1:nrow(imputedDataListNa[[dataListIter]])) {
    #     if(is.na(imputedDataListNa[[dataListIter]][k,'TIMEPOINT'])) {
    #       imputenum = imputenum+1
    #       evidence = imputedDataListNa[[dataListIter]][k,]
    #       #prob = BnExactInferenceCumulated(fitList,evidence)[dataListIter]
    #       prob = BnExactInference(fitList[[dataListIter]],evidence)
    #       if(dataListIter>1) {
    #         probS = BnExactInferenceCumulated(fitList,evidence)[(dataListIter-1)]
    #       }else {
    #         probS=1
    #       }
    #       id = imputedDataListNa[[dataListIter]][k,'id']
    #       #if(runif(1,0,1)>prob) {
    #       if(0.5>prob) {
    #         imputedDataListNa[[dataListIter]][k,'TIMEPOINT'] = 1
    #         if(dataListIter<length(fitList)) {
    #           for(paddingIter in (dataListIter+1):length(fitList)) {
    #             imputedDataListNa[[paddingIter]][imputedDataListNa[[paddingIter]][,'id']==id,'TIMEPOINT'] = 1
    #             imputedDataListNa[[paddingIter]][imputedDataListNa[[paddingIter]][,'id']==id,'PREVTIMEPOINT'] = 1
    #           }
    #           #paddingIter = dataListIter+1
    #           #imputedDataListNa[[paddingIter]][imputedDataListNa[[paddingIter]][,'id']==id,'PREVTIMEPOINT'] = 0
    #         }
    #       }else {
    #         imputedDataListNa[[dataListIter]][k,'TIMEPOINT'] = 0
    #         if(dataListIter<length(fitList)) {
    #           paddingIter = dataListIter+1
    #           imputedDataListNa[[paddingIter]][imputedDataListNa[[paddingIter]][,'id']==id,'PREVTIMEPOINT'] = 0
    #         }
    #         # for(paddingIter in dataListIter:1) {
    #         #   imputedDataListNa[[paddingIter]][imputedDataListNa[[paddingIter]][,'id']==id,'TIMEPOINT'] = 0
    #         #   imputedDataListNa[[paddingIter]][imputedDataListNa[[paddingIter]][,'id']==id,'PREVTIMEPOINT'] = 0
    #         # }
    #       }
    #       # if(runif(1,0,1)>probS) {
    #       #   imputedDataListNa[[dataListIter]][k,'PREVTIMEPOINT'] = 1
    #       # }else {
    #       #   imputedDataListNa[[dataListIter]][k,'PREVTIMEPOINT'] = 0
    #       # }
    #     }
    #   }
    # }
    #print(imputenum)
    
    # likelihoodTestdata = dataListNA
    # likelihood = 0
    # for(dataListIter in 1:length(fitList)) {
    #   for(k in 1:nrow(likelihoodTestdata[[dataListIter]])) {
    #     instance = likelihoodTestdata[[dataListIter]][k,]
    #     if(instance['TIMEPOINT']==1 & instance['PREVTIMEPOINT']==0 &instance['delta']==1) {
    #       prob = 1-BnExactInference(fitList[[dataListIter]],instance)
    #       likelihood = likelihood + prob
    #     }else if(is.na(instance['TIMEPOINT'])&!is.na(instance['PREVTIMEPOINT'])&instance['delta']==0) {
    #       prob = BnExactInferenceCumulated(fitList,instance)[dataListIter]
    #       likelihood = likelihood + prob
    #     }
    #   }
    # }
    # print(likelihood)
    
    if(weighted) {
      fitList = weightedLearning(dataListNA,fitList,prevFitList,timePoints,originalTraining,mod,kmMod,bnCurveList)
    }
    prevFitList = fitList

    
    if(cont) {
      contFittingList = contFitting(dataListCont,fitList,variableList,timePoints)
    }
    
  }
  
  print('start predict')
  #prediction
  survivalFunctionTesting = predictFunction(fitList,fitListParents,nbList,testingList,originalTesting,timePoints,queryMethod,dataList,contFittingList,cont,variableList,kmMod)
  survivalFunctionTraining = predictFunction(fitList,fitListParents,nbList,trainingList,originalTraining,timePoints,queryMethod,dataList,contFittingList,cont,variableList,kmMod)
  
  testCurvesToReturn = survivalFunctionTesting
  timesAndCensTest = cbind.data.frame(time = originalTesting$time, delta = originalTesting$delta)
  timesAndCensTrain = cbind.data.frame(time = originalTraining$time, delta = originalTraining$delta)
  trainingCurvesToReturn = survivalFunctionTraining
  
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}

predictFunction <- function(fitList,fitListParents,nbList,testingList,originalTesting,timePoints,queryMethod,dataList,contFittingList,cont,variableList,kmMod) {
  numTimepoint = length(timePoints)
  numReturnNA = 0
  numNotDecreasing = 0

  #first value in cpt
  previousTimepointProb = rep(1,nrow(testingList[[1]]))

  survivalFunction <- data.frame(matrix(ncol = nrow(testingList[[1]]), nrow = numTimepoint))
  for(i in 1:numTimepoint) {
    testing = testingList[[i]]
    testing$time = NULL
    testing$delta = NULL
    
    tempFit = fitList[[i]]

    tempNb = nbList[[i]]
    cat(i)
    cat(' ')

    for(j in 1:nrow(testing)) {
      #print(nrow(testing))
      #print(ncol(testCurvesToReturn))
      
      if(queryMethod == 'lw') {
        eviList = as.list(testing[j,])
        prob = cpquery(tempFit, event = (TIMEPOINT == 0), evidence = as.list(testing[j,]),method = 'lw')
      }else if(queryMethod == 'predict') {
        evidence = testing[j,]
        predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "bayes-lw", prob = TRUE)
        attr(predicted, "prob")
        prob = attr(predicted, "prob")[1]
      }else if(queryMethod == 'exact') {
        evidence = testing[j,]
        if(i>1){kmprob = predict(kmMod,timePoints[i])/predict(kmMod,timePoints[i-1])}else{kmprob = predict(kmMod,timePoints[i])}
        prob = BnExactInference(tempFit,evidence,kmprob=kmprob)
      }else{
        print('No such query method')
      }
      
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

  colnames(survivalFunction) = 1:nrow(testing)
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
    
    #iss = 5
    #iss_dead = iss * (dead/(alive+dead+0.5*censored))

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
      iss_dead = iss * (dead/(alive+dead+0.5*censored))
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

weightedLearning <- function(dataList,fitList,prevFitList,timePoints,originalTraining,mod,kmMod,bnCurveList,weight=NULL) {
  
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
        newcpt[flatIndex] = newcpt[flatIndex] + weight[k]
      }
      
      newcpt[is.na(newcpt)] = 0
      iss = 5
      newcpt = newcpt + iss
      
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

imputeFunction = function(fitted1,fitted2,inputdata) {
  #n=floor(50*nrow(inputdata[!is.na(inputdata),])/nrow(inputdata[is.na(inputdata),]))
  n=1
  if(n<1){n=1}
  if(n>20){n=20}
  imputeCount = 0
  ouputdata = inputdata[complete.cases(inputdata),]
  ouputdata = ouputdata[1,]
  for(k in 1:nrow(inputdata)) {
    if(is.na(inputdata[k,'TIMEPOINT'])) {
    #if(T) {
      imputeCount = imputeCount+1
      evidence = inputdata[k,]
      datainstance = inputdata[k,]
      evidence$TIMEPOINT = NULL
      if(is.null(fitted1)) {
        prob1 = 1
      }else {
        prob1 = BnExactInference(fitted1,evidence)
      }
      prob2 = BnExactInference(fitted2,evidence)
      imputeProb = prob2/prob1
      if(is.nan(prob1)){prob1=0}
      if(is.nan(imputeProb)){imputeProb=1}
      if(imputeProb>1){imputeProb=1}
      if(imputeProb<0){imputeProb=0}
      for(w in 1:n) {
        if(runif(1,0,1)<imputeProb) {
          datainstance['TIMEPOINT'] = 0
          ouputdata = rbind(ouputdata,datainstance)
        }else {
          datainstance['TIMEPOINT'] = 1
          ouputdata = rbind(ouputdata,datainstance)
        }
      }
    }else {
      for(w in 1:n) {
        ouputdata = rbind(ouputdata,inputdata[k,])
      }
    }
    # if(!is.na(inputdata[k,'TIMEPOINT'])) {
    #   for(w in 1:n) {
    #     ouputdata = rbind(ouputdata,inputdata[k,])
    #   }
    # }
  }
  print(imputeCount)
  return(ouputdata)
}
