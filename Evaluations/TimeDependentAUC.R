library(survival)
#We use this for the prodlim function.
library(prodlim)
#Helper Functions: predictMeanSurvivalTimeSpline(survivalCurve,predictedTimes)
source("Evaluations/EvaluationHelperFunctions.R")
library(pROC)

TimeDependentAUC = function(survMod, timePoints, weighting="marginal", iid=TRUE){
  #Being passed an empty model.
  if(is.null(survMod)) return(NULL)
  #Being passed a model that failed.
  suppressWarnings(if(is.na(survMod[[1]])) return(NULL))
  predictedTimes = survMod[[1]]$time
  survivalCurves = survMod[[1]][-1]
  trueDeathTimes = survMod[[2]]$time
  censorStatus = survMod[[2]]$delta
  censorTimes = trueDeathTimes[as.logical(1-censorStatus)]
  trainingDeathTimes = survMod[[3]]$time
  trainingCensorStatus = survMod[[3]]$delta

  AUC_Vec = rep(0,length(timePoints))
  
  for(i in 1:length(timePoints)) {
    tk = timePoints[i]
    survivalCurvesAtRisk = survivalCurves[trueDeathTimes>tk | censorStatus==1]
    trueDeathTimesAtRisk = trueDeathTimes[trueDeathTimes>tk | censorStatus==1]
    D = trueDeathTimesAtRisk<=tk
    predicted = apply(survivalCurvesAtRisk,2, function(z) predictProbabilityFromCurve(z,predictedTimes,tk))
    roc_obj = roc(D, predicted)
    AUC_Vec[i] = auc(roc_obj)
  }
  
  return(AUC_Vec)
}
