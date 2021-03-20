source('Evaluations/EvaluationHelperFunctions.R')
keras_cv = function(x,y,weights,cvFoldIndex) {
  lambdaList = c(0.0001,0.001,0.01,0.1)
  lossList = rep(nrow(y),length(lambdaList))
  cat('internal cross validation: ')
  for(lambdaIter in 1:length(lambdaList)) {
    cat(lambdaIter);cat(' ');
    for(cvIter in 1:length(cvFoldIndex)) {
      cv_x <- x[-cvFoldIndex[[cvIter]],]
      cv_y <-  y[-cvFoldIndex[[cvIter]]]
      cv_w = weights[-cvFoldIndex[[cvIter]]]
      cv_x_v <- x[cvFoldIndex[[cvIter]],]
      cv_y_v <-  y[cvFoldIndex[[cvIter]]]
      cv_w_v = weights[cvFoldIndex[[cvIter]]]
      
      
      model <- keras_model_sequential()
      model %>%
        layer_dense(units=1, activation='sigmoid', input_shape=c(ncol(cv_x)),
                    kernel_regularizer = regularizer_l1(l = lambdaList[lambdaIter]))
      #summary(model)
      model %>% compile(
        loss = "binary_crossentropy",
        optimizer = 'sgd',
        metrics = "binary_accuracy"
      )
      fitted <- model %>% fit(
        cv_x, cv_y,
        sample_weight = cv_w,
        epochs = 20, batch_size = 32,
        validation_data = list(cv_x_v,cv_y_v,cv_w_v),
        verbose = 0
      )
      lossList[lambdaIter] = tail(fitted$metrics$val_loss,1)
    }
  }
  bestLambda = lambdaList[which.min(lossList)]
  
  return(bestLambda)
}

keras_cvInt = function(x,y,w,cvFoldIndex) {
  lambdaList = c(0.001,0.01,0.05,0.1,0.5,1)
  lossList = rep(0,length(lambdaList))
  cat('internal cross validation: ')
  for(lambdaIter in 1:length(lambdaList)) {
    cat(lambdaList[lambdaIter]);cat(' ');
    for(cvIter in 1:length(cvFoldIndex)) {
      cv_x = x[-cvFoldIndex[[cvIter]],]
      cv_y = y[-cvFoldIndex[[cvIter]],]
      cv_w = w[-cvFoldIndex[[cvIter]],]
      cv_x_v = x[cvFoldIndex[[cvIter]],]
      cv_y_v = y[cvFoldIndex[[cvIter]],]
      cv_w_v = w[cvFoldIndex[[cvIter]],]
      
      
      my_regularizer_wrapper_cv <- custom_metric("reg", function(x){my_regularizer(x, lambda1=lambdaList[lambdaIter],lambda2=0)})

      model_cv <- keras_model_sequential()
      model_cv %>%  layer_dense(units=2*ncol(cv_y), activation='sigmoid', input_shape=c(ncol(cv_x)),use_bias=TRUE,
                             #bias_initializer = initializer_random_normal(-2.19,0.5),
                             #kernel_regularizer = regularizer_l1(l = 0.5),
                             kernel_regularizer = my_regularizer_wrapper_cv
                             #bias_regularizer=my_bias_regularizer
      )
      #summary(model)
      optimizerSGD = optimizer_sgd(
        lr = 2.0,
        momentum = 0.5,
        decay = 0.0,
        nesterov = FALSE,
        clipnorm = 1,
        clipvalue = 1
      )
      model_cv %>% compile(
        loss = customBCE,
        optimizer = optimizerSGD,
        metrics = NULL
        #sample_weight_mode='temporal',
        #weighted_metrics = 'binary_accuracy'
      )
      early_stopping = callback_early_stopping(monitor='loss', patience=300, verbose=2)
      reduceLearningRate = callback_reduce_lr_on_plateau(monitor='loss', patience=30, factor=0.9,min_lr=0.01, verbose=0)
      
      fitted <- model_cv %>% fit(
        cv_x, cbind(cv_y,cv_w),
        epochs = 6000, batch_size = 16,
        validation_split = 0,
        #validation_data = list(cv_x_v,cbind(cv_y_v,cv_w_v)),
        shuffle=T,
        verbose = 0,
        callbacks = list(early_stopping,reduceLearningRate)
      )
      
      # probInt = model_cv %>% predict(cv_x_v)
      # probInt = probInt[,1:ncol(cv_y_v)]
      # survival = hazard2survival(1-probInt)
      # pmf = survival
      # for(i in 1:(ncol(survival)-1)) {
      #   pmf[,i] = survival[,i] - survival[,i+1]
      # }
      # loglike = 0
      # for(k in 1:nrow(cv_x_v)) {
      #   if(sum(cv_y_v[k,]*cv_w_v[k,])==1) {
      #     loglike = loglike + sum(pmf[k,]*cv_y_v[k,]*cv_w_v[k,])+0.00001
      #   }else if(sum(cv_y[k,]*cv_w[k,])==0.5) {
      #     loglike = loglike + sum(survival[k,]*cv_y_v[k,]*cv_w_v[k,])+0.00001
      #     print('hi')
      #   }else {
      #     loglike = loglike + sum(survival[k,][cv_w_v[k,] == 0.5])+0.00001
      #   }
      # }

      probInt = model_cv %>% predict(cv_x_v)
      probInt = probInt[,1:ncol(cv_y_v)]
      loss = sum((cv_y_v*log(probInt+0.00001) + (1-cv_y_v)*log(1-probInt+0.00001))*cv_w_v)
      #lossList[lambdaIter] = lossList[lambdaIter]+tail(fitted$metrics$val_loss,1)
      lossList[lambdaIter] = lossList[lambdaIter]+loss
      #print(loss)
    }
    print(lossList[lambdaIter])
  }
  bestLambda = lambdaList[which.max(lossList)]
  
  return(bestLambda)
}

EMdata = function(x,y,w,data,model,timePoints) {
  survivalFunction <- data.frame(matrix(nrow = nrow(x), ncol = ncol(y)))
  hazardFunction <- data.frame(matrix(nrow = nrow(x), ncol = ncol(y)))
  previousTimepointProb = rep(1,nrow(x))
  probInt = 1 - model %>% predict(x)
  probInt = probInt[,1:ncol(y)]
  for(i in 1:ncol(y)) {
    prob = probInt[,i]
    hazardFunction[,i] = prob
    survivalFunction[,i] = prob*previousTimepointProb
    previousTimepointProb = survivalFunction[,i]
  }

  oldy = y
  copy_y = y
  weight = w
  copyWeight = matrix(0,nrow(w),ncol(w))
  for(k in 1:nrow(y)) {
    for(i in 1:ncol(y)) {
      if(w[k,i]<1&data[k,'delta']==0) {
        oldy[k,i] = 0
        copy_y[k,i] = 1
        survivalCurve = c(1,survivalFunction[k,])
        survivalCurveTime = c(0,timePoints)
        if(i>1) {
          a = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,timePoints[i-1])
        }else {
          a = 1
        }
        b = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,timePoints[i])
        c = predictProbabilityFromCurve(survivalCurve,survivalCurveTime,data[k,'time'])+0.0001
        if(i>1) {
          if(w[k,i-1]<1) {
            survprob = survivalFunction[k,i-1]
            #survprob = a/c
          }else {
            survprob = 1
          }
        }else {
          survprob = 1
        }
        prob = hazardFunction[k,i]
        # if(i>1) {
        #   prob = 1- ((a/c)-(b/c))
        # }else {
        #   prob = b/c
        # }
        weight[k,i] = survprob*prob
        copyWeight[k,i] = survprob*(1-prob)
      }
    }
  }
  newx = rbind(x,x)
  newy = rbind(oldy,copy_y)
  neww = rbind(weight,copyWeight)
  
  return(list(x=newx,y=newy,w=neww))
}

EMdataKM = function(x,y,w,data,kmMod,timePoints) {
  
  copy_y = y
  weight = w
  copyWeight = matrix(0,nrow(w),ncol(w))
  for(k in 1:nrow(y)) {
    for(i in 1:ncol(y)) {
      if(w[k,i]<1&data[k,'delta']==0) {
        y[k,i] = 0
        copy_y[k,i] = 1
        if(i>1) {
          a = predict(kmMod,timePoints[i-1])
        }else {
          a = 1
        }
        b = predict(kmMod,timePoints[i])
        c = predict(kmMod,data[k,'time'])
        if(i>1) {
          if(w[k,i-1]<1) {
            #survprob = survivalFunction[k,i-1]
            survprob = a/c
          }else {
            survprob = 1
          }
        }else {
          survprob = 1
        }
        #prob = hazardFunction[k,i]
        if(i>1) {
          prob = 1- ((a/c)-(b/c))
        }else {
          prob = b/c
        }
        weight[k,i] = survprob*prob
        copyWeight[k,i] = survprob*(1-prob)
      }
    }
  }
  newx = rbind(x,x)
  newy = rbind(y,copy_y)
  neww = rbind(weight,copyWeight)
  
  return(list(x=newx,y=newy,w=neww))
}
predictFunctionLR <- function(fitList,testing,timePoints) {
  numTimepoint = length(timePoints)
  numReturnNA = 0
  numNotDecreasing = 0
  testing[,c('time','delta')] = NULL
  
  previousTimepointProb = rep(1,nrow(testing))
  
  survivalFunction <- data.frame(matrix(ncol = nrow(testing), nrow = numTimepoint))
  for(i in 1:numTimepoint) {
    if(i>length(fitList)){
      survivalFunction[i,] = previousTimepointProb
      break
    }
    fitted = fitList[[i]]
    if(is.null(fitted)) {
      prob = 1
    }else {
      testing = as.matrix(testing)
      prob = 1 - fitted %>% predict(testing)
      #prob = 1 - predict(fitted,newx=testing,type='response',s="lambda.min")
    }
    survivalFunction[i,] = prob*previousTimepointProb
    
    previousTimepointProb = survivalFunction[i,]
  }
  if(numReturnNA>0) {cat('return NA: ',numReturnNA)}
  if(numNotDecreasing>0) {cat('Not decreasing: ',numNotDecreasing)}
  
  colnames(survivalFunction) = 1:nrow(testing)
  return(survivalFunction)
}

predictFunctionLRInt <- function(model,testing,timePoints) {
  numTimepoint = length(timePoints)
  testing[,c('time','delta')] = NULL
  
  survivalFunction <- data.frame(matrix(ncol = nrow(testing), nrow = numTimepoint))
  previousTimepointProb = rep(1,nrow(testing))
  testing = as.matrix(testing)
  probInt = 1 - model %>% predict(testing)
  if(anyNA(probInt)) {print('NaN in keras prediction')}
  #probInt[is.na(probInt)] = 1
  probInt = probInt[,1:length(timePoints)]
  for(i in 1:numTimepoint) {
    prob = probInt[,i]
    survivalFunction[i,] = prob*previousTimepointProb
    previousTimepointProb = survivalFunction[i,]
  }
  colnames(survivalFunction) = 1:nrow(testing)
  return(survivalFunction)
}



prepareDataKeras = function(data,timePoints) {
  
  y = matrix(0,nrow(data),length(timePoints))
  
  for(i in 1:length(timePoints)) {
    singleTime = rep(0,nrow(data))
    singleTime[data$time <= timePoints[i] & data$delta == 1] <- 1
    singleTime[data$time <= timePoints[i] & data$delta == 0] <- NA
    y[,i] =  singleTime
  }
  x = data
  x[,c('time','delta','id')] = NULL
  x = as.matrix(x)
  
  w = matrix(1,nrow(y),ncol(y))
  for(i in 2:length(timePoints)) {
    singleweight = rep(1,nrow(y))
    singleweight[y[,i]==1 & y[,i-1]==0 ] <- 1
    singleweight[y[,i]==1 & y[,i-1]==1 ] <- 0
    singleweight[is.na(y[,i]) & y[,i-1]==0] <- 0.5
    singleweight[is.na(y[,i]) & is.na(y[,i-1])] <- 0
    w[,i] =  singleweight
  }
  y[is.na(y)] = 0
  
  return(list(x=x,y=y,w=w))
}

hazard2survival = function(hazard) {
  survivalFunction <- matrix(ncol = ncol(hazard), nrow = nrow(hazard))
  previousTimepointProb = rep(1,nrow(hazard))
  for(i in 1:ncol(hazard)) {
    survivalFunction[,i] = hazard[,i]*previousTimepointProb
    previousTimepointProb = survivalFunction[,i]
  }
  return(survivalFunction)
}

makeMod = function(TestCurves,TrainCurves, timePoints, training, testing) {
  
  #timePoints = fixtime(timePoints)
  survivalFunctionTesting = TestCurves
  survivalFunctionTesting= rbind(rep(1,nrow(testing)),survivalFunctionTesting)
  testCurvesToReturn = cbind(time = c(0,timePoints), survivalFunctionTesting)
  
  #testCurvesToReturn = cbind.data.frame(time = timePoints, survivalProbabilitiesTest) 
  timesAndCensTest = cbind.data.frame(time = testing$time, delta = testing$delta)
  timesAndCensTrain = cbind.data.frame(time = training$time, delta = training$delta)
  
  survivalFunctionTraining = TrainCurves
  survivalFunctionTraining= rbind(rep(1,nrow(training)),survivalFunctionTraining)
  trainingCurvesToReturn = cbind(time = c(0,timePoints), survivalFunctionTraining)
  #trainingCurvesToReturn = cbind.data.frame(time = timePoints, survivalProbabilitiesTrain) 
  
  curveCheck(testCurvesToReturn)
  curveCheck(trainingCurvesToReturn)
  
  return(list(TestCurves = testCurvesToReturn, TestData = timesAndCensTest,TrainData = timesAndCensTrain,TrainCurves= trainingCurvesToReturn))  
  
}

Evaluation = function(EMMod) {
  survivalPredictionMethod = 'Median'
  bayesConc = Concordance(EMMod, 'None',survivalPredictionMethod)
  bayesBrierInt = BrierScore(EMMod, type = "Integrated", numPoints =  1000, integratedBrierTimes = NULL)
  bayesL1 = L1(EMMod, 'Margin', F,survivalPredictionMethod)
  bayesDcal = DCalibrationCumulative(list(EMMod),10)
  return(list(Conc=bayesConc,BrierInt=bayesBrierInt,L1=bayesL1,Dcal=bayesDcal))
}
