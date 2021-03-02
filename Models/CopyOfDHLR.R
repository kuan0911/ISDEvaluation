#------------------------------------------------------------
#Cumulative model with seperate network at each timepoint
#include all data in structure learning. Na if dead or censored
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
library(glmnet)
source("Models/cleanDataForBayesNet.R")
source("Models/bayesianNetHelper.R")
source("Models/fractionallyInclude.R")
source("Models/changesource.R")
source("Models/kerasHelper.R")
source("ValidateCleanCV/createFoldsAndNormalize.R")
library(keras)
library(tensorflow)
library(magrittr)
library(dplyr)
library(kerastuneR)


DHLR = function(training, testing, timePoints = 0,debug = FALSE){
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
  
  variableList = variableTypes(originalTraining,10)
  
  print(timePoints)
  
  numTimepoint = length(timePoints)
  
  dataListNA = createTimeSplit(training,timePoints,includeNa = T)
  
  print(length(dataListNA))
  
  fitList <- vector("list", numTimepoint)
  print('learn starting graph')
  
  res = prepareDataKeras(training,timePoints)
  x=res$x
  y=res$y
  w=res$w
  
  for(iter in 1:1) {
    break
    for(i in 1:numTimepoint) {
      if(isTRUE(debug)){cat(i)
        cat(' ')}
      
      data = dataListNA[[i]]
      
      dataHazardNA = data
      dataHazardNA = dataHazardNA[is.na(dataHazardNA$PREVTIMEPOINT)|dataHazardNA$PREVTIMEPOINT==0,]
      resKM = weightedImputeKM(kmMod,dataHazardNA,timePoints,i)
      weightsNAKM = resKM$weight
      dataHazardNAKM = resKM$data
      resLR = weightedImputeLR(fitList,dataHazardNA,timePoints,i,curveCache)
      weightsNALR = resLR$weight
      dataHazardNALR = resLR$data
      
      dataHazard = data
      dataHazard = dataHazard[!is.na(dataHazard$PREVTIMEPOINT),]
      dataHazard = dataHazard[dataHazard$PREVTIMEPOINT == 0,]
      weights = rep(1,nrow(dataHazard))
      weights[is.na(dataHazard$TIMEPOINT)] = 0.5
      dataHazard[is.na(dataHazard$TIMEPOINT),'TIMEPOINT'] = 0
      
      if(iter==1) {
        trainingData = dataHazard
        weights = weights
      }else{
        trainingData = dataHazardNALR
        weights = weightsNALR
      }
      
      
      cvFoldIndex = createFoldsOfData(trainingData, numberOfFolds=5)[[1]]
      transCvFoldIndex = rep(1,nrow(trainingData))
      for(cvFoldIter in 1:length(cvFoldIndex)){transCvFoldIndex[cvFoldIndex[[cvFoldIter]]]=cvFoldIter}
      
      trainingData[,c('time','delta','id','PREVTIMEPOINT')] = NULL
      # y = trainingData$TIMEPOINT
      # y <- as.integer(y)-1
      # trainingData$TIMEPOINT = NULL
      # x = as.matrix(trainingData)
      
      
      
      if(sum(y == 0)<4 | sum(y == 1)<4) {
        fitList[[i]] = NULL
        cat('skip ')
      }else {
        #fitted = cv.glmnet(x,y,family='binomial', weights=weights,foldid=transCvFoldIndex)
        #fitted <- glmnet(x,y,family='binomial',weights=weights,lambda=0)
        
        #bestLambda = keras_cv(x,y,weights,cvFoldIndex)
        #customBCEIntegrated_wrapper <- custom_metric("wbce", function(y_true, y_pred){customBCE(y_true, y_pred, weights=weights)})
        
        
        model <- keras_model_sequential()
        model %>%
          layer_dense(units=1, activation='sigmoid', input_shape=c(ncol(x)),
                      kernel_regularizer = regularizer_l1(l = 0.5)
                      #bias_regularizer=my_bias_regularizer
          )
        #summary(model)
        optimizer = optimizer_sgd(
          lr = 0.0001,
          momentum = 0.8,
          nesterov = FALSE,
          clipnorm = NULL,
          clipvalue = NULL
        )
        model %>% compile(
          loss = 'binary_crossentropy',
          optimizer = 'sgd',
          metrics = "binary_accuracy"
        )
        fitted <- model %>% fit(
          x, y[,i],
          epochs = 100, batch_size = 8,
          validation_split = 0,
          shuffle=T,
          verbose = 0,
          sample_weight = w[,i]
        )
        
        
        fitList[[i]] <- model
      }
      
    }
  }
  
  oldx = x
  oldy = y
  oldw = w
  
  if(anyNA(x)|anyNA(y)|anyNA(w)) {print('Missing value in training data')}
  
  # cvFoldIndex = createFoldsOfData(training, numberOfFolds=5)[[1]]
  # transCvFoldIndex = rep(1,nrow(training))
  # for(cvFoldIter in 1:length(cvFoldIndex)){transCvFoldIndex[cvFoldIndex[[cvFoldIter]]]=cvFoldIter}
  # bestLambda = keras_cvInt(x,y,w,cvFoldIndex)
  
  bestLambda = 0.05
  cat('best lambda: ');cat(bestLambda);
  
  
  
  for(iter in 1:1) {
    print('EM')
    
    customBCEIntegrated_wrapper <- custom_metric("wbce", function(y_true, y_pred){customBCE(y_true, y_pred, weights=w)})
    aveAccuracy_wrapper <- custom_metric("ave_acc", function(y_true, y_pred){aveAccuracy(y_true, y_pred, weights=w)})
    my_regularizer_wrapper <- custom_metric("reg", function(x){my_regularizer(x, weights=w, lambda1=bestLambda)})
    my_bias_regularizer_wrapper <- custom_metric("regb", function(x){my_bias_regularizer(x, weights=w)})
    
    model <- keras_model_sequential()
    model %>%  layer_dense(units=2*ncol(y), activation='sigmoid', input_shape=c(ncol(x)),use_bias=TRUE,
                           #bias_initializer = initializer_random_normal(-2.19,0.5),
                           #kernel_regularizer = regularizer_l1(l = 0.5),
                           #kernel_constraint = constraint_maxnorm(max_value = 10*ncol(x), axis = 0),
                           #bias_constraint = constraint_maxnorm(max_value = 10*ncol(x), axis = 0),
                           kernel_regularizer = my_regularizer_wrapper
                           #bias_regularizer = my_bias_regularizer_wrapper
    )
    #summary(model)
    optimizer = optimizer_sgd(
      lr = 5.0,
      momentum = 0.5,
      decay = 0.0,
      nesterov = FALSE,
      clipnorm = 1.0,
      clipvalue = 1.0
    )
    model %>% compile(
      loss = customBCEIntegrated_wrapper,
      optimizer = optimizer,
      metrics = NULL
      #sample_weight_mode='temporal',
      #weighted_metrics = 'binary_accuracy'
    )
    
    early_stopping = callback_early_stopping(monitor='loss', patience=3000, verbose=2)
    reduceLearningRate = callback_reduce_lr_on_plateau(monitor='loss', patience=30, factor=0.9,min_lr=0.001, verbose=0)
    
    fitted <- model %>% fit(
      x, cbind(y,w),
      epochs =10000, batch_size = 256,
      validation_split = 0,
      shuffle=T,
      verbose = 0,
      callbacks = list(reduceLearningRate,early_stopping)
    )
    res = EMdata(oldx,oldy,oldw,training,model,timePoints)
    x=res$x
    y=res$y
    w=res$w
  }
  
  
  
  
  print('start predict')
  #prediction
  #survivalFunctionTesting = predictFunctionLR(fitList,originalTesting,timePoints)
  #survivalFunctionTraining = predictFunctionLR(fitList,training,timePoints)
  survivalFunctionTesting = predictFunctionLRInt(model,originalTesting,timePoints)
  survivalFunctionTraining = predictFunctionLRInt(model,originalTraining,timePoints)
  
  testCurvesToReturn = survivalFunctionTesting
  timesAndCensTest = cbind.data.frame(time = originalTesting$time, delta = originalTesting$delta)
  timesAndCensTrain = cbind.data.frame(time = originalTraining$time, delta = originalTraining$delta)
  trainingCurvesToReturn = survivalFunctionTraining
  
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
}

customBCE = function(y_true, y_pred, weights){
  K <- backend()
  numTimepoints = ncol(weights)
  
  weightsF = y_true[,(numTimepoints+1):(2*numTimepoints)]
  y_trueF = y_true[,1:numTimepoints]
  y_predF = y_pred[,1:numTimepoints]
  loss = (y_trueF)*K$log(y_predF+0.00001) + (1.0-y_trueF)*K$log(1.0-y_predF+0.00001)
  weighted_loss = loss*weightsF
  
  # y_trueDummy = y_true[,(numTimepoints+1):(2*numTimepoints)]
  # y_predDummy = y_pred[,(numTimepoints+1):(2*numTimepoints)]
  # lossDummy = K$mean((y_trueDummy)*K$log(y_predDummy) + (1.0-y_trueDummy)*K$log(1.0-y_predDummy))
  
  
  return(K$mean(-weighted_loss))
}

aveAccuracy = function(y_true, y_pred, weights){
  K <- backend()
  numTimepoints = ncol(weights)
  weightsF = y_true[,(numTimepoints+1):(2*numTimepoints)]
  y_trueF = y_true[,1:numTimepoints]
  y_predF = y_pred[,1:numTimepoints]
  
  acc_vec = K$equal(y_trueF,K$round(y_predF))
  acc_vec_w = K$cast_to_floatx(acc_vec)*weightsF
  acc_sum = K$sum(acc_vec_w)
  weights_sum = K$cast_to_floatx(K$sum(weightsF))
  
  acc = acc_sum/weights_sum
  return(acc)
}

my_regularizer = function(x,weights,lambda1){
  K <- backend()
  xF = x[,1:ncol(weights)]
  xW = x[,(ncol(weights)+1):ncol(x)]
  
  shiftx = xF[,2:(ncol(xF))]
  concat = xF[,ncol(xF)]
  concat2D = K$expand_dims(concat)
  shiftx_concat = K$concatenate(list(shiftx,concat2D))
  sub = K$abs(shiftx_concat - xF)
  
  # 
  # weights_col = K$sum(weights,axis=as.integer(0),keepdims = F)
  # weights_sum_all = K$cast_to_floatx(K$sum(weights))
  # 
  # sub_sum = K$sum(sub,axis=as.integer(0),keepdims = F)
  # weighted_sub = sub_sum*K$cast_to_floatx(weights_col)
  # 
  # x_abs = K$abs(xF)
  # x_sum = K$sum(x_abs,axis=as.integer(0),keepdims = F)
  # weighted_x = x_sum*K$cast_to_floatx(weights_col)
  
  #regabs = 0.05*(K$sum(weighted_sub)/weights_sum_all)
  regabs = 0.01*(K$mean(sub))
  #regweight = 0.0005*(K$sum(weighted_x)/weights_sum_all)
  
  regweight = lambda1*K$mean(K$abs(xF))
  regconstrain = K$mean(K$abs(xW))
  return(regweight+regconstrain+regabs) 
}

my_bias_regularizer = function(x,weights){
  K <- backend()
  xW = x[(ncol(weights)+1):length(x)]
  sub = xW + 0.5
  
  xF = x[1:ncol(weights)]
  # regabs = 0.0001*K$mean(K$abs(sub))
  
  regconstrain = 0.5*K$mean(K$abs(sub))
  reg = 0.0001*K$mean(K$abs(xF))
  return(regconstrain) 
}







