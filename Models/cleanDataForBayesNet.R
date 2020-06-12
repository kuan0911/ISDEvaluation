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
    while(nrow(remainData)>10) {
      res = quantile(remainData$time,1/numTimepoint)
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
    print(method)
    print(timesplit)
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
          breakList = c(-Inf, unique(quantile(allData[[c]], probs = breaks[2:discretize_level]/100)), Inf)
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

createTimeSplit <- function(data, timesplit,assumption = F) {
  
  dataList <- vector("list", length(timesplit))
  
  previousTimepoint = factor(integer(nrow(data)))
  
  #TIMEPOINT = one means dead
  for(i in 1:length(timesplit)) {
    tempData <- data
    timepoint = integer(nrow(data))
    if(assumption == F) {
      timepoint[data$time < timesplit[i]] <- 1
    }else {
      timepoint[data$time < timesplit[i] & data$delta == 1] <- 1
    }
    
    tempData$TIMEPOINT <- factor(timepoint)
    tempData$PREVTIMEPOINT <- previousTimepoint
    previousTimepoint = tempData$TIMEPOINT
    #CENSORED == 'false' or tempData$TIMEPOINT == not dead
    if(assumption == F) {
      tempData <- tempData[tempData$delta == 1 | tempData$TIMEPOINT == 0,]
    }
    tempData$time <- NULL
    tempData$delta <- NULL
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
    tempData = factor(timepoint,levels = c('0','1'))
    newData[,paste0("t_",i)] =  tempData
    
    #CENSORED == 'false' or tempData$TIMEPOINT == not dead
    #make sure TIMEPOINT has two levels
    #levels(timeData[,i]) = c('0','1')
  }
  
  return(newData)
}

