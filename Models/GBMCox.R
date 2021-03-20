#### File Information #####################################################################################################################
#File Name: coxPH_KP.R
#Date Created: May 26, 2018
#Author: Humza Haider
#Email: hshaider@ualberta.ca

### General Comments ######################################################################################################################
#This file is used to run the Cox proportional hazards model and then use the  Kalbfleisch-Prentice method to attain survival probabilities.
#Additionally, the elastic net cox is also implemented here.

### Functions #############################################################################################################################

## Function 1: CoxPH_KP(training, testing,ElasticNet =F, numFolds =5)

#Inputs:
#   training:   The training dataset (after normalization and imputation).
#   testing:    The testing dataset (after normalization and imputation).
#   ElasticNet: Boolean indicating to use or not use elastic net cox.
#   numFolds:   Number of folds for internal cross validation (elastic-net cox only).

# Output: A list of 4 items:(1) TestCurves - The survival curves for the testing set.
#                           (2) TestData - The censor/death indicator and event time for the testing set. 
#                           (3) TrainData - The censor/death indicator and event time for the training set. 
#                           (4) TrainCurves - The survival curves for the training set.

# Usage: Train and evaluate the Cox-KP/CoxEN-KP model.

## Function 2: KPEstimator(lp,lpTime,censorStatus){

# Inputs:
#   lp:           The linear predictors from the coxEN model (for training data).
#   lpTime:       The training times.
#   censorStatus: A vector indicating the censor status of the patients.

# Output: The baseline survival function via the kalbfleisch-prentice estimator.

# Usage: This will build the baseline survival function for patients (because cox only considers hazards we need to use a method
#         of estimating the baseline survival function). (Helpher function for CoxPH_KP).

### Code #############################################################################################################################
#Library Dependencies:
#survival is needed to get survfit and the implimentation of the KP estimator.
library(survival)
#For EN-Cox
library(fastcox)
#For sindex
library(prodlim)

library(gbm)

GBMCox_KP = function(training, testing,ElasticNet=F, numFolds = 5){
  print('Gradien boost')
  timeInd = which(names(training) == "time")
  deltaInd = which(names(training) == "delta")

  coxModel = gbm(Surv(time,delta)~.,       # formula
                 data=training,                 # dataset
                 #weights=w,
                 #var.monotone=c(0,0,0),     # -1: monotone decrease, +1: monotone increase, 0: no monotone restrictions
                 distribution="coxph",
                 n.trees=3000,              # number of trees
                 shrinkage=0.001,           # shrinkage or learning rate, 0.001 to 0.1 usually work
                 interaction.depth=3,       # 1: additive model, 2: two-way interactions, etc
                 bag.fraction = 0.5,        # subsampling fraction, 0.5 is probably best
                 train.fraction = 1,      # fraction of data for training, first train.fraction*N used for training
                 cv.folds = 5,              # do 5-fold cross-validation
                 n.minobsinnode = 10,       # minimum total weight needed in each node
                 keep.data = TRUE)
  best.iter <- gbm.perf(coxModel,method='cv')
  
  linearPredictionsTraining = predict(coxModel,training[,-c(timeInd, deltaInd)],best.iter)
  linearPredictionsTesting = predict(coxModel,testing[,-c(timeInd, deltaInd)],best.iter)
  survivalEstimate = KPEstimator(linearPredictionsTraining, training$time,training$delta)
  survCurvs = t(sapply(survivalEstimate[[2]], function(x) x^exp(linearPredictionsTesting)))
  survCurvsTraining = t(sapply(survivalEstimate[[2]], function(x) x^exp(linearPredictionsTraining)))
  
  survivalCurves = list(time = survivalEstimate[[1]], surv = survCurvs)
  survivalCurvesTrain = list(time = survivalEstimate[[1]], surv = survCurvsTraining)
  #If 0 wasnt included in the timepoints we would like to manually add it with a survival probability of 1.
  if(0 %in% survivalCurves$time){
    timePoints = survivalCurves$time
    probabilities = survivalCurves$surv
    
    probabilitiesTrain = survivalCurvesTrain$surv
  } else{
    timePoints = c(0,survivalCurves$time)
    probabilities = rbind(1,survivalCurves$surv)
    
    probabilitiesTrain = rbind(1,survivalCurvesTrain$surv)
  }
  curvesToReturn = cbind.data.frame(time = timePoints, probabilities)
  trainingCurvesToReturn = cbind.data.frame(time = timePoints, probabilitiesTrain)
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)
  return(list(TestCurves = curvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
  
}

#Currently this estimator uses parameters which assume new ties. Future implementations should adjust the code for this (though estimates)
#are nearly identical when using the ties-adjusted model.
KPEstimator = function(lp,lpTime,censorStatus){
  indexToKeep = sindex(sort(lpTime), unique(sort(lpTime)))
  orderLPTime = order(lpTime)
  cumHaz = rev(cumsum(rev(exp(lp[orderLPTime]))))
  alpha = ((1-(exp(lp[orderLPTime])/cumHaz)))^exp(-lp[orderLPTime])
  survivalFunc = cumprod(alpha^censorStatus[orderLPTime])
  return(list(time = lpTime[orderLPTime][indexToKeep],surv = survivalFunc[indexToKeep]))
}

