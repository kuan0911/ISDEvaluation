source("Models/cleanDataForBayesNet.R")
source("ValidateCleanCV/validateAndClean.R")

sythesize = function(n=5000) {
  survivalDataset <- read.csv(file = "Data/All_Data_updated_may2011_CLEANED.csv")
  names(survivalDataset)[c(1,2)] = c("time", "delta")
  survivalDataset$delta = 1 - survivalDataset$delta 
  
  survivalDataset = validateAndClean(survivalDataset, T)
  
  m = floor(sqrt(nrow(survivalDataset))+1)
  #if(m>20) {m=20}
  m = 5
  #timePoints = seq(0,quantile(survivalDataset$time,0.85),(quantile(survivalDataset$time,0.85)-0)/m )
  quantileVals = seq(0,1,length.out = m+2)[-c(1,m+2)]
  timePoints = unname(quantile(survivalDataset[survivalDataset$delta==1,]$time, quantileVals))
  timePoints = timePoints[!duplicated(timePoints)]
  fillup = seq(max(timePoints),max(survivalDataset$time),tail(timePoints, n=2)[2]-tail(timePoints, n=2)[1])
  timePoints = c(timePoints,fillup)
  timePoints = timePoints[!duplicated(timePoints)]
  m = length(timePoints)
  
  cleanedsurvivalDataset = cleanDataType(0,survivalDataset,allDiscrete=T,discretize_level=6)
  survivalDatasetList = createTimeSplit(cleanedsurvivalDataset$train,timePoints,assumption = F,includeNa = F)
  
  fitList <- vector("list", length(timePoints))
  dagList <- vector("list", length(timePoints))
  
  for(i in 1:length(timePoints)) {
    cat(i)
    cat(' ')
    data = survivalDatasetList[[i]]
    
    # for(k in 1:nrow(data)) {
    #   if(!is.na(data[k,'PREVTIMEPOINT'])) {
    #     if(data[k,'PREVTIMEPOINT']==1) {
    #       data[k,'TIMEPOINT'] = NA
    #     }
    #   }
    # }
    
    
    data$time = NULL
    data$delta = NULL
    data$id = NULL
    dataFit = data[data$PREVTIMEPOINT == 0,]
    data = data[data$PREVTIMEPOINT == 0,]
    data$PREVTIMEPOINT = NULL
    dataFit$PREVTIMEPOINT = NULL
    
    dag = structural.em(data, maximize = "hc",maximize.args = list(restart=100,blacklist=NULL,whitelist=NULL), fit = "bayes",fit.args=list(iss=2),impute='bayes-lw',return.all = T,start = NULL, max.iter = 1, debug = FALSE)
    
    print(parents(dag$dag,'TIMEPOINT'))
    print(children(dag$dag,'TIMEPOINT'))
    
    data = data[data$PREVTIMEPOINT == 0,]
    
    fit = bn.fit(dag$dag, dataFit, method='bayes',iss=10)
    fitList[[i]] <- fit
    #plotDag(dag$dag)
  }
  
  covariateData = survivalDatasetList[[1]]
  covariateData$TIMEPOINT = NULL
  covariateData$PREVTIMEPOINT = NULL
  covariateData$time = NULL
  covariateData$delta = NULL
  covariateData$id = NULL
  
  covariateDag = hc(covariateData)
  covariateFit = bn.fit(covariateDag, covariateData, method='bayes',iss=5)
  covariateSim = rbn(covariateFit, n, covariateData)
  
  time = rep(NA,nrow(covariateSim))
  timePointsWithZero = c(0,timePoints)
  for(i in 1:length(timePoints)) {
    cat(i)
    cat(' ')
    for(k in 1:nrow(covariateSim)) {
      if(is.na(time[k])) {
        covariate = covariateSim[k,]
        predicted = predict(fitList[[i]], "TIMEPOINT", covariate, method = "bayes-lw", prob = TRUE)
        if(runif(1,0,1)>attr(predicted, "prob")[1]) {
          time[k] = runif(1,timePointsWithZero[i], timePointsWithZero[i+1])
        }
      }
    }
  }
  delta = rep(as.integer(1),nrow(covariateSim))
  syntheticData = cbind(time,delta,covariateSim)
  
  for(k in 1:nrow(syntheticData)) {
    if(is.na(syntheticData[k,'time'])) {
      syntheticData[k,'delta'] = as.integer(0)
      syntheticData[k,'time'] = runif(1,timePoints[m],2*timePoints[m]-timePoints[m-1])
    }
  }
  
  for(k in 1:nrow(syntheticData)) {
    if(runif(1,0,1)<0.3) {
      censoreTime = runif(1,0.1,timePoints[m])
      if(censoreTime<syntheticData[k,'time']) {
        syntheticData[k,'delta'] = as.integer(0)
        syntheticData[k,'time'] = censoreTime
      }
    }
  }
  
  write.csv(syntheticData,'Data/syntheticData.csv')
  
}


