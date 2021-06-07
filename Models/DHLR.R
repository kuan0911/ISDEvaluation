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
source("Models/kerasLossFunction.R")
source("ValidateCleanCV/createFoldsAndNormalize.R")
library(keras)
library(tensorflow)
library(magrittr)
library(dplyr)


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
  
  # res = prepareDataKeras(training,timePoints)
  # x=res$x
  # y=res$y
  # w=res$w

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

  
  
  
  
  
  
  bestLambda = 0.1
  bestLambda = bestLambda * 0.66
  cvFoldIndex = createFoldsOfData(training, numberOfFolds=3)[[1]]
  #bestLambda = cvGeneral(training,timePoints,numberOfFolds=3)
  
  res = prepareDataKeras(training,timePoints)
  x=res$x
  y=res$y
  w=res$w
  bestLambda2 = 0.05
  #bestLambda2 = keras_cvInt_lambda2(x,y,w,cvFoldIndex,bestLambda)
  cat('best lambda2: ');cat(bestLambda2);
  
  #bestLambda = keras_cvInt(x,y,w,cvFoldIndex)
  
  cat('best lambda: ');cat(bestLambda);cat('  ');
  # bestLambda = cvGeneral(training,timePoints,numberOfFolds=3)
  # quantileVals = seq(0,1,length.out = bestLambda+2)[-c(1,bestLambda+2)]
  # timePoints = unname(quantile(training$time, quantileVals))
  # timePoints = timePoints[!duplicated(timePoints)]
  # print(timePoints)
  # k_clear_session()
  
  res = prepareDataKeras(training,timePoints)
  x=res$x
  y=res$y
  w=res$w
  oldx = x
  oldy = y
  oldw = w
  
  if(F) {
    print('start with KM imputation')
    resKM = EMdataKM(oldx,oldy,oldw,training,kmMod,timePoints)
    x=resKM$x
    y=resKM$y
    w=resKM$w
  }
  
  #y=1-y
  
  
  for(iter in 1:1) {
    print('EM')
    if(anyNA(x)|anyNA(y)|anyNA(w)) {print('Missing value in training data')}
    allHazard = colSums(oldy*oldw)/colSums(oldw)
    
    customBCEIntegrated_wrapper <- custom_metric("wbce", function(y_true, y_pred){customBCE(y_true, y_pred)})
    #aveAccuracy_wrapper <- custom_metric("ave_acc", function(y_true, y_pred){aveAccuracy(y_true, y_pred, weights=w)})
    my_regularizer_wrapper <- custom_metric("reg", function(x){my_regularizer(x, lambda1=bestLambda,lambda2=bestLambda2,weights=w)})
    my_bias_regularizer_wrapper <- custom_metric("regb", function(x){my_bias_regularizer(x, allHazard=log(allHazard))})
    
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
      lr = 0.1,
      momentum = 0.1,
      decay = 0.0,
      nesterov = FALSE,
      clipnorm = 5.0,
      clipvalue = 5.0
    )
    optimizerADAGRAD = optimizer_adagrad(
      lr = 5.0
    )
    model %>% compile(
      loss = customBCEIntegrated_wrapper,
      optimizer = optimizer,
      metrics = NULL
      #sample_weight_mode='temporal',
      #weighted_metrics = 'binary_accuracy'
    )
    
    if(iter==1) {
      patience = 20
      batch_size = 32
    }else {
      patience = 20
      batch_size = 32
    }
    
    early_stopping = callback_early_stopping(monitor='loss', patience=2000, verbose=2)
    reduceLearningRate = callback_reduce_lr_on_plateau(monitor='loss', patience=patience, factor=0.9,min_lr=0.01, verbose=0)
    
    history <- model %>% fit(
      x, cbind(y,w),
      epochs =10000, batch_size = batch_size,
      validation_split = 0,
      shuffle=T,
      verbose = 0,
      callbacks = list(early_stopping,reduceLearningRate)
    )
    
    probInt = model %>% predict(oldx)
    probInt = probInt[,1:ncol(oldy)]
    loss = sum((oldy*log(probInt+0.00001) + (1-oldy)*log(1-probInt+0.00001))*oldw)
    print('loss: ')
    print(loss)
    
    TestCurves = predictFunctionLRInt(model,originalTesting,timePoints)
    TrainCurves = predictFunctionLRInt(model,originalTraining,timePoints)
    EMMod = makeMod(TestCurves,TrainCurves, timePoints, training, testing)
    resEval = Evaluation(EMMod)
    print(resEval$Dcal)
    combinedBins =colSums(ldply(lapply(seq_along(list(EMMod)), function(x) getBinned(list(EMMod)[[x]], 10)), rbind))
    barplot(combinedBins,horiz=TRUE,xlab=iter)

    res = EMdata(oldx,oldy,oldw,training,model,timePoints)
    x=res$x
    y=res$y
    w=res$w
  }

  layerkernalweights = as.matrix(model$weights[[1]])
  layerkernalweights[abs(layerkernalweights)<0.001]=0
  layerbiasweights = as.matrix(model$weights[[2]])
  rownames(layerkernalweights) = colnames(x)
  write.csv(layerkernalweights,'layerkernalweights.csv')
  
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
  
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn,timePoints=timePoints))  
}


  

  
  
