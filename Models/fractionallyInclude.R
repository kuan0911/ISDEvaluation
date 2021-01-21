source('Evaluations/EvaluationHelperFunctions.R')

weightedImputeKM = function(kmMod,inputdata,timePoints,currentTime) {
  if(is.null(inputdata$PREVTIMEPOINT)) {
    print('PREVTIMEPOINT missing when imputing weight')
    inputdata$PREVTIMEPOINT = rep(0,nrow(inputdata))
  }
  imputeCount = 0
  outputdata = inputdata
  weight = rep(1,nrow(outputdata))
  for(k in 1:nrow(inputdata)) {
    if(is.na(inputdata[k,'TIMEPOINT'])) {
      imputeCount = imputeCount+1
      datainstance = inputdata[k,]
      evidence = inputdata[k,]
      evidence$TIMEPOINT = NULL
      if(!is.na(inputdata[k,'PREVTIMEPOINT'])|currentTime==1) {
        probSurvive = 1
      }else {
        probSurvive = predict(kmMod,timePoints[currentTime-1])/predict(kmMod,datainstance$time)
        #probSurvive = predict(kmMod,timePoints[currentTime-1])
        #probSurvive = runif(1,0,1)
      }
      if(is.na(inputdata[k,'PREVTIMEPOINT'])) {
        prob = 1 - ((predict(kmMod,timePoints[currentTime-1])/predict(kmMod,datainstance$time))-(predict(kmMod,timePoints[currentTime])/predict(kmMod,datainstance$time)))
      }else {
        prob = predict(kmMod,timePoints[currentTime])/predict(kmMod,datainstance$time)
      }
      #prob = runif(1,0,1)
      if(is.nan(probSurvive) | is.nan(prob)) {print('weightedImputeKM: Error imputing: prob is nan')}
      if(probSurvive>1 |probSurvive<0){print('weightedImputeKM: probSurvive error')}
      if(prob>1 |prob<0){print('weightedImputeKM: prob error')}
      
      outputdata[k,'TIMEPOINT'] = 0
      weight[k] = probSurvive*prob
      
      datainstance['TIMEPOINT'] = 1
      outputdata = rbind(outputdata,datainstance)
      weight = c(weight,probSurvive*(1-prob))
    }
  }
  #print(imputeCount)
  #outputdata[,c('PREVTIMEPOINT','time','delta','id')] = NULL
  if(anyNA(outputdata$TIMEPOINT)) {print('Warning: data imputed contain NA value')}
  return(list(data=outputdata,weight=weight))
}

weightedImputeLR = function(fitList,inputdata,timePoints,currentTime,curveCache) {
  if(is.null(fitList[[length(fitList)]])){
    print('weightedImputeLR: fitList is null. Return input data')
    return(inputdata)
  }
  if(is.null(inputdata$PREVTIMEPOINT)) {
    print('PREVTIMEPOINT missing when imputing weight')
    inputdata$PREVTIMEPOINT = rep(0,nrow(inputdata))
  }
  imputeCount = 0
  outputdata = inputdata
  weight = rep(1,nrow(outputdata))
  
  for(k in 1:nrow(inputdata)) {
    if(is.na(inputdata[k,'TIMEPOINT'])) {
      imputeCount = imputeCount+1
      datainstance = inputdata[k,]
      evidence = inputdata[k,]
      evidence[c('PREVTIMEPOINT','TIMEPOINT','time','delta','id')] = NULL
      
      

      # if(!is.na(inputdata[k,'PREVTIMEPOINT'])|currentTime==1) {probSurvive = 1}
      # else {probSurvive = predictLRCurve(fitList,evidence,currentTime-1)}
      # prob = 1 - predict(fitList[[currentTime]],newx=as.matrix(evidence),type='response',s="lambda.min")

      survivalCurve = c(1,predictLRCurveAll(fitList,evidence))
      survivalCurveTime = c(0,timePoints)
      if(!is.na(inputdata[k,'PREVTIMEPOINT'])|currentTime==1) {
        probSurvive = 1
      }else {
        a = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,timePoints[currentTime-1])
        c = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,datainstance$time)
        probSurvive = a/c
        #probSurvive = predictLRCurve(fitList,evidence,currentTime-1)/predictLRCurve(fitList,evidence,censorTime)
      }
      if(currentTime>1) {
        if(is.na(inputdata[k,'PREVTIMEPOINT'])) {
          a = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,timePoints[currentTime-1])
          b = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,timePoints[currentTime])
          c = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,datainstance$time)
          prob = 1- ((a/c)-(b/c))
          #prob = 1 - ((predictLRCurve(fitList,evidence,currentTime-1)/predictLRCurve(fitList,evidence,censorTime))-(predictLRCurve(fitList,evidence,currentTime)/predictLRCurve(fitList,evidence,censorTime)))
        }else {
          b = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,timePoints[currentTime])
          c = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,datainstance$time)
          prob = b/c
          #prob = predictLRCurve(fitList,evidence,currentTime)/predictLRCurve(fitList,evidence,censorTime)
        }
      }else {
        b = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,timePoints[currentTime])
        c = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,datainstance$time)
        prob = b/c
        #prob = predictLRCurve(fitList,evidence,currentTime)
      }

      #prob = runif(1,0,1)
      if(is.nan(probSurvive) | is.nan(prob)) {print('weightedImputeLR: Error imputing: prob is nan')}
      if(probSurvive>1 |probSurvive<0){print('weightedImputeLR: probSurvive error')}
      if(prob>1 |prob<0){print('weightedImputeLR: prob error')}
      
      outputdata[k,'TIMEPOINT'] = 0
      weight[k] = probSurvive*prob
      
      if(is.na(inputdata[k,'PREVTIMEPOINT'])) {
        datainstance['TIMEPOINT'] = 1
        outputdata = rbind(outputdata,datainstance)
        weight = c(weight,probSurvive*(1-prob))
      }
    }
  }
  #print(imputeCount)
  #outputdata[,c('PREVTIMEPOINT','time','delta','id')] = NULL
  if(anyNA(outputdata$TIMEPOINT)) {print('Warning: data imputed contain NA value')}
  return(list(data=outputdata,weight=weight))
}

predictLRCurve = function(fitList,evidence,currentTime) {
  if(currentTime==0){return(1)}
  cumulatedProb = 1
  curve = rep(1,currentTime)
  evidence$PREVTIMEPOINT = NULL
  evidence$TIMEPOINT = NULL
  evidence$time = NULL
  evidence$delta = NULL
  evidence$id = NULL
  
  for(i in 1:currentTime) {
    fitted = fitList[[i]]
    if(is.null(fitted)|is.null(evidence)) {print('predictLRCurve: error. lenth not match')}
    #evidence[c('PREVTIMEPOINT','TIMEPOINT','time','delta','id')] = NULL
    evidence = as.matrix(evidence)
    prob = 1 - predict(fitted,newx=evidence,type='response',s="lambda.min")
    cumulatedProb = cumulatedProb * prob
    curve[i] = cumulatedProb
  }
  return(curve[currentTime])
}

predictLRCurveAll = function(fitList,evidence) {
  cumulatedProb = 1
  curve = rep(1,length(fitList))
  evidence$PREVTIMEPOINT = NULL
  evidence$TIMEPOINT = NULL
  evidence$time = NULL
  evidence$delta = NULL
  evidence$id = NULL
  
  for(i in 1:length(fitList)) {
    fitted = fitList[[i]]
    if(is.null(fitted)|is.null(evidence)) {print('predictLRCurve: error. lenth not match')}
    #evidence[c('PREVTIMEPOINT','TIMEPOINT','time','delta','id')] = NULL
    evidence = as.matrix(evidence)
    prob = 1 - predict(fitted,newx=evidence,type='response',s="lambda.min")
    cumulatedProb = cumulatedProb * prob
    curve[i] = cumulatedProb
  }
  return(curve)
}