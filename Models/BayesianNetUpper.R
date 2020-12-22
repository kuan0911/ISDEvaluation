#-------------------
#cross validation
#-------------------
source("Models/BayesianNet.R")
source("Models/cleanDataForBayesNet.R")

BayesianNetUpper = function(training,testing){
  print('start internal cross validation')
  #bestC1 = internalCV_BayesianNet(training,5,10)
  bestC1 = 1
  print('bestC1')
  print(bestC1)
  
  #m = floor(sqrt(nrow(training))+1)
  #m = floor(nrow(training)^(1/3)+1)
  #if(m>20) {m=20}
  m = 10
  quantileVals = seq(0,1,length.out = m+2)[-c(1,m+2)]
  #quantileVals = seq(0,1,length.out = m+2)[-1]
  timePoints = unname(quantile(training[training$delta==1,]$time, quantileVals))
  timePoints = timePoints[!duplicated(timePoints)]
  fillup = seq(max(timePoints),max(training$time),(tail(timePoints, n=2)[2]-tail(timePoints, n=2)[1])*2)
  #timePoints = c(timePoints,fillup)
  #timePoints = c(timePoints,quantile(training$time,0.99))
  timePoints = timePoints[!duplicated(timePoints)]
  #timePoints = timeSplitFunction(numTimepoint = 10, method = 'fixrate', data = training, debug=T)
  #timePoints = timePoints[-1]
  #kmMod = prodlim(Surv(time,delta)~1, data = training)
  #timePoints = kmTimesplit(m,kmMod,training,step=1)
  m = length(timePoints)
 
  
  #timePoints = timeSplitFunction(numTimepoint = m, method = 'quantile', data = training, debug=T)
  
  mod = BayesianNet(training, testing, timePoints,debug=T)
  
  # m=20
  # quantileVals = seq(0,1,length.out = m+2)[-c(1,m+2)]
  # #quantileVals = seq(0,1,length.out = m+2)[-1]
  # timePoints = unname(quantile(training$time, quantileVals))
  # timePoints = timePoints[!duplicated(timePoints)]
  
  #timePoints = fixtime(timePoints)
  survivalFunctionTesting = mod$TestCurves
  survivalFunctionTesting= rbind(rep(1,nrow(testing)),survivalFunctionTesting)
  testCurvesToReturn = cbind(time = c(0,timePoints), survivalFunctionTesting)
  
  #testCurvesToReturn = cbind.data.frame(time = timePoints, survivalProbabilitiesTest) 
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)

  survivalFunctionTraining = mod$TrainCurves
  survivalFunctionTraining= rbind(rep(1,nrow(training)),survivalFunctionTraining)
  trainingCurvesToReturn = cbind(time = c(0,timePoints), survivalFunctionTraining)
  #trainingCurvesToReturn = cbind.data.frame(time = timePoints, survivalProbabilitiesTrain) 
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}



internalCV_BayesianNet = function(training, numFolds, m=0){
  
  C1Vec = c(0.05,0.1,0.3,0.5,0.8,1,2,4,6,10)
  
  foldedData = createFoldsOfData(training, numFolds)[[2]] 
  resultsMatrix = matrix(rep(0,numFolds*length(C1Vec)), ncol = length(C1Vec),nrow =numFolds) #7 vals of C1
  
  #Much of this code is duplicated from MTLR(). See comments in that function.
  for(i in 1:numFolds){
    trainingFold = foldedData[[1]][[i]]
    testingFold = foldedData[[2]][[i]]
    
    trainingFold = trainingFold[order(trainingFold$delta),]
    testingFold = testingFold[order(testingFold$delta),]
    
    if(m==0) {
      m = floor(sqrt(nrow(training))+1)
    }
    quantileVals = seq(0,1,length.out = m+2)[-c(1,m+2)]
    timePoints = unname(quantile(trainingFold$time, quantileVals))
    timePoints = timePoints[!duplicated(timePoints)]
    d = as.matrix(trainingFold[,-c(1,2)])
    dAsZero = matrix(0,ncol = ncol(d), nrow = nrow(d))
    yval = matrix(1 -Reduce(c,Map(function(x) trainingFold$time > timePoints[x],1:length(timePoints))), ncol = nrow(trainingFold), byrow = T)
    
    #Train biases first and then all parameters.
    resultVec = c()
    for(C1 in C1Vec){
      survivalCurves = BayesianNet(trainingFold, testingFold,C1,timePoints,debug=F)$TestCurves
      
      #print(C1)
      #print(avgLogLikLossBayesian(survivalCurves,testingFold,timePoints))
      resultVec = c(resultVec,avgLogLikLossBayesian(survivalCurves,testingFold,timePoints))
    }
    resultsMatrix[i,] = resultVec
  }
  meanResults = apply(resultsMatrix, 2, mean)
  bestC1 = C1Vec[which.min(meanResults)]
  return(bestC1)
}


avgLogLikLossBayesian = function(survivalCurves, dat, timePoints){
  #For the log-likelihood loss we need to compute losses differently for censored and uncensored patients.
  #For censored patients the loss will correspond the the (log) survival probability assigned by the model at the time of censoring.
  #For uncensored patients, we will consider the log of the probability assigned to the time interval when the patient died. 
  #Then we take the negative of this loss and thus would like to minimize the loss.
  NCens = sum(1-  dat$delta)
  
  logloss = 0
  #Censored patients
  censorTimes = dat$time[1:NCens]
  probAtCensorTime = sapply(seq_along(censorTimes),
                            function(index) predictProbabilityFromCurve(survivalCurves[,index],
                                                                        timePoints,
                                                                        censorTimes[index]))
  logloss = logloss - sum(log(probAtCensorTime + 1e-5))
  
  #Uncensored patients
  deathTimes = dat$time[(NCens+1):nrow(dat)]
  uncenSurvival = as.matrix(survivalCurves[,(NCens+1):nrow(dat),drop=FALSE])
  uncenSurvival = rbind(1,uncenSurvival,0)
  pmfProbs = -diff(uncenSurvival)
  indexRow = sapply(deathTimes, function(x) findInterval(x, timePoints)) + 1
  indexCol = 1:(nrow(dat) - NCens)
  indexMat = matrix(c(indexRow,indexCol),ncol = 2)
  probs = pmfProbs[indexMat]
  logloss = logloss  - sum(log(probs))
  
  return(logloss/nrow(dat))
}

fixtime <- function(timePoints) {
  newTimePoints = timePoints
  newTimePoints[1] = timePoints[1]/2
  for(i in 2:length(timePoints)) {
    newTimePoints[i] = (timePoints[i-1] + timePoints[i])/2
  }
  return(newTimePoints)
}

