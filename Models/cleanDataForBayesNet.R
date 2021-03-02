timeSplitFunction <- function(numTimepoint = 10, method = 'quantile', data, debug = FALSE){
  if(method == 'even') {
    timesplit <- seq(0,quantile(data$time,0.85),(quantile(data$time,0.85)-0)/(numTimepoint-1) )
  }else if(method == 'quantile'){
    timeBreaks = seq(0, 100,100/numTimepoint)
    timesplit <- quantile(data[data$delta == 1,]$time,timeBreaks[-1]/100)
  }else if(method == 'quantile2'){
    timeBreaks = seq(0, 100,100/numTimepoint)
    timesplit <- quantile(data$time,timeBreaks[-1]/100)
  }else if(method == 'fixrate'){
    timesplit<-vector()
    remainData = data[data$delta == 1,]
    #remainData = data
    while(nrow(remainData)>50) {
      res = quantile(remainData[remainData$delta==1,]$time,1/numTimepoint)
      timesplit = c(timesplit, res)
      remainData = remainData[remainData$time > res,]
    }
    print('time split num remain data')
    print(nrow(remainData))
    print(length(timesplit))
  }else {
    print('No such time split method')
  }
  if(debug) {
    #print(method)
    #print(timesplit)
  }
  return(timesplit)
}

cleanDataType <- function(numTesting, allData, allDiscrete=T, discretize_level = 4){
  maxFactorLevel = 10
  keepOriginalLevels = T
  
  if(keepOriginalLevels == F){
    print('Do not keep original levels')
  }
  
  #Change datatype to factor
  if(allDiscrete == F) {
    for(c in names(allData)) {
      if(dim(table(allData[[c]])) < maxFactorLevel) {
        allData[[c]] = factor(allData[[c]])
      }else {
        allData[[c]] = as.numeric(allData[[c]])
      }
    }
  }else {
    print('Turn all data into descrete')
    timepreserved = allData$time
    deltapreserved = allData$delta
    #allData$delta = NULL
    for(c in names(allData)) {
      if(keepOriginalLevels) {
        if(dim(table(allData[[c]])) < maxFactorLevel) {
          allData[[c]] = factor(allData[[c]])
        }else {
          breaks = seq(0,100,100/discretize_level)
          quantileVals = seq(0,1,length.out = discretize_level+1)[-c(1,discretize_level+1)]
          breakList = c(-Inf, unique(quantile(allData[[c]],quantileVals)), Inf)
          
          #breakList = seq(quantile(allData[[c]],0.10),quantile(allData[[c]],0.90),(quantile(allData[[c]],0.90)-quantile(allData[[c]],0.10))/(discretize_level-2) )
          #breakList = c(-Inf,breakList,Inf)
          
          allData[[c]] = cut(allData[[c]],breaks = breakList,labels = 1:(length(breakList)-1))
        }
      }else{
        if(dim(table(allData[[c]])) <= discretize_level) {
          allData[[c]] = factor(allData[[c]])
        }else {
          breaks = seq(0,100,100/discretize_level)
          breakList = c(-Inf, unique(quantile(allData[[c]], probs = breaks[2:discretize_level]/100)), Inf)
          allData[[c]] = cut(allData[[c]],breaks = breakList,labels = 1:(length(breakList)-1))
        }
      }
    }
    allData$time = timepreserved
    allData$delta = factor(deltapreserved)
  }
  
  #drop column
  #for(c in names(allData)) {
  #  if(dim(table(allData[[c]])) == 1) {
  #    allData[[c]] = NULL
  #  }
  #}
  
  testing = allData[1:numTesting,]
  training = allData[(numTesting+1):nrow(allData),]
  
  return(list(test = testing, train = training))
}

createTimeSplit <- function(data, timesplit,assumption = F,includeNa = F) {
  
  data$id = 1:nrow(data)
  
  dataList <- vector("list", length(timesplit))
  
  previousTimepoint = factor(integer(nrow(data)))
  
  #TIMEPOINT = one means dead
  for(i in 1:length(timesplit)) {
    tempData <- data
    timepoint = integer(nrow(data))
    if(assumption == F & includeNa == F) {
      timepoint[data$time < timesplit[i]] <- 1
    }else if(includeNa == T){
      timepoint[data$time < timesplit[i] & data$delta == 1] <- 1
      timepoint[data$time < timesplit[i] & data$delta == 0] <- NA
    }else {
      timepoint[data$time < timesplit[i] & data$delta == 1] <- 1
    }
    
    tempData$TIMEPOINT <- factor(timepoint)
    tempData$PREVTIMEPOINT <- previousTimepoint
    previousTimepoint = tempData$TIMEPOINT
    tempData$TIMEPOINT[tempData$PREVTIMEPOINT == 0 & tempData$TIMEPOINT == 1 & tempData$delta == 0] = NA
    #CENSORED == 'false' or tempData$TIMEPOINT == not dead
    if(assumption == F & includeNa == F) {
      tempData <- tempData[tempData$delta == 1| is.na(tempData$TIMEPOINT) | tempData$TIMEPOINT == 0,]
    }
    
    #make sure TIMEPOINT has two levels
    levels(tempData$TIMEPOINT) = c('0','1')
    levels(tempData$PREVTIMEPOINT) = c('0','1')
    dataList[[i]] <- tempData
  }
  
  return(dataList)
}

createTimeSplitIntegrated <- function(data, timesplit) {
  
  newData = data
  
  #TIMEPOINT = one means dead
  for(i in 1:length(timesplit)) {
    timepoint = integer(nrow(data))
    timepoint[data$time < timesplit[i] & data$delta == 1] <- 1
    timepoint[data$time < timesplit[i] & data$delta == 0] <- NA
    tempData = factor(timepoint,levels = c('0','1'))
    newData[,paste0("t_",i)] =  tempData
    
    #CENSORED == 'false' or tempData$TIMEPOINT == not dead
    #make sure TIMEPOINT has two levels
    #levels(timeData[,i]) = c('0','1')
  }
  
  return(newData)
}

variableTypes = function(data,maxFactorLevel) {
  variableList = c()
  for(c in names(data)) {
    if(dim(table(data[[c]])) < maxFactorLevel) {
      variableList[c] = FALSE
    }else {
      variableList[c] = TRUE
    }
  }
  return(variableList)
}

kmTimesplit = function(m,kmMod,data,step=1) {
  timePoints = rep(0,m)
  curves = rep(0,nrow(data))
  windowMax = max(data$time)
  summation = 0
  for(k in 1:nrow(data)) {
    if(data[k,'delta']==1) {
      summation = summation + 1
      curves[k] = 1
    }else {
      summation = summation + 1 - predict(kmMod,windowMax)/predict(kmMod,data[k,'time'])
      curves[k] = 0
    }
  }
  avesum = summation/m
  prevTime = 0
  sortData = data[order(data$time),]
  for(i in 1:m) {
    for(s in seq(prevTime, windowMax, step)) {
      localsum = 0
      timesetflag = F
      for(k in 1:nrow(sortData)) {
        if(data[k,'delta']==1) {
          if(data[k,'time']>=prevTime && data[k,'time']<s) {localsum = localsum + 1}
          #print(localsum)
        }else {
          if(data[k,'time']<s) {
            localsum = localsum + 1 - predict(kmMod,s)/predict(kmMod,prevTime)
          }
        }
        if(localsum >= avesum) {
          timePoints[i] = s
          prevTime = s
          timesetflag = T
          print(i)
          print('time set')
          #print(timePoints[i])
          break
        }
      }
      if(timesetflag) {break}
      #print(localsum)
    }
  }
  return(timePoints)
}

kmTimesplitV2 = function(m,kmMod,data,step=0.1) {
  timePoints = rep(0,m)
  windowMax = max(data$time)
  
  summation = 0
  for(k in 1:nrow(data)) {
    if(data[k,'delta']==1) {
      summation = summation + 1
    }else {
      summation = summation + predict(kmMod,data[k,'time']) - predict(kmMod,windowMax)
    }
  }
  
  avesum = summation/m
  prevTime = 0
  sortData = data[order(data$time),]
  for(i in 1:m) {
    for(s in seq(prevTime, windowMax, step)) {
      localsum = 0
      timesetflag = F
      for(k in 1:nrow(sortData)) {
        if(data[k,'delta']==1) {
          if(data[k,'time']>=prevTime && data[k,'time']<s) {localsum = localsum + 1}
        }else {
          if(data[k,'time']<=prevTime) {
            localsum = localsum + (predict(kmMod,prevTime) - predict(kmMod,s))/predict(kmMod,data[k,'time'])
          }else if(data[k,'time']<s & data[k,'time']>prevTime) {
            localsum = localsum + (predict(kmMod,data[k,'time']) - predict(kmMod,s))/predict(kmMod,data[k,'time'])
          }
        }
        if(localsum >= avesum) {
          timePoints[i] = s
          prevTime = s
          timesetflag = T
          #print(timePoints[i])
          break
        }
      }
      if(timesetflag) {break}
      #print(localsum)
    }
  }
  timePoints = timePoints[timePoints!=0]
  return(timePoints)
}
  

kmTimeSplitSimple = function(data, p=0.1,step=0.2, minSurvProb = 0.01) {
  timePoints <- vector()
  kmMod = prodlim(Surv(time,delta)~1, data = data)
  targetProb = 1-p
  preProb=1

  for(t in seq(0, max(data$time), step)) {
    if(predict(kmMod,abs(t-step))>targetProb & predict(kmMod,t)<targetProb) {
      timePoints = c(timePoints, t)
      targetProb = targetProb * (1-p)
      if(targetProb<minSurvProb) {break}
    }
  }
  return(timePoints)
}

