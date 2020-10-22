predictSingle <- function(fitList,childrenList,timesplit,evidence,time,method='naivebayes') {
  evidence$time = NULL
  evidence$TIMEPOINT = NULL
  evidence$PREVTIMEPOINT = NULL
  evidence$delta = NULL
  #testing$PREVTIMEPOINT = factor(integer(nrow(testing)),levels = c('0','1'))
  threshold = 0.000001
  
  #first value in cpt
  previousTimepointProb = 1
  
  survivalFunction = rep(0,length(fitList))
  for(i in 1:length(fitList)) {
    tempFit = fitList[[i]]
    
    if(method == 'lw') {
      eviList = as.list(evidence)
      prob = cpquery(tempFit, event = (TIMEPOINT == 0), evidence = as.list(testing[j,]),method = 'lw')
    }else if(method == 'predict') {
      predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "parents", prob = TRUE, n=500)
      attr(predicted, "prob")
      prob = attr(predicted, "prob")[1]
    }else if(method == 'naivebayes') {
      predicted = predict(tempFit, node="TIMEPOINT", evidence,method = "parents", prob = TRUE)
      ap = attr(predicted, "prob")
      
      p0 = 1
      p1 = 1
      variables = childrenList[[i]]
      for(v in variables) {
        e = as.integer(evidence[[v]])
        t = tempFit[[v]][['prob']]
        temp_evidence = evidence
        temp_evidence[[v]] = NULL
        temp_evidence$TIMEPOINT = factor(0,levels = c('0','1'))
        t1 = predict(tempFit, node=v, temp_evidence,method = "parents", prob = TRUE)
        t1 = attr(t1, "prob")
        temp_evidence$TIMEPOINT = factor(1,levels = c('0','1'))
        t2 = predict(tempFit, node=v, temp_evidence,method = "parents", prob = TRUE)
        t2 = attr(t2, "prob")
        
        if(t1[e]>threshold) {rawp0=t1[e]} else {rawp0=threshold}
        if(t2[e]>threshold) {rawp1=t2[e]} else {rawp1=threshold}
        p0 = p0 * rawp0
        p1 = p1 * rawp1
      }
      #print(t)
      p0 = p0 * ap[1]
      p1 = p1 * ap[2]
      p0_Prob = p0/(p0+p1)
      
      prob = p0_Prob
    }else{
      print('No such query method')
    }
    
    survivalFunction[i] = prob*previousTimepointProb
    if(time<timesplit[i]) {
      toReturn = previousTimepointProb + (survivalFunction[i]-previousTimepointProb) * (time - timesplit[i-1])/(timesplit[i] - timesplit[i-1])
      return(toReturn)
    }
    
    
    if(is.na(survivalFunction[i])) {
      survivalFunction[i] = previousTimepointProb-0.01
      if(survivalFunction[i] <0){survivalFunction[i] = 0}
    }
    
    if(previousTimepointProb<survivalFunction[i]){
      survivalFunction[i] = previousTimepointProb[j]-0.0001
      if(survivalFunction[i] <0){survivalFunction[i] = 0}
      #print('probability not decreasing')
    }
    
    previousTimepointProb = survivalFunction[i]
  }
  return(survivalFunction[length(fitList)])
}

convertToFlatIndex = function(index,dimension) {
  flatIndex = 0
  index = as.integer(index)
  for(digit in 1:length(index)) {
    flatIndex = flatIndex + (index[digit]-1)*prod(dimension[0:(digit-1)])
  }
  flatIndex = flatIndex+1
}

convertToFlatIndexMultiple = function(index,dimension) {
  flatIndexVec = c()
  index = as.integer(index)
  naDigit = is.na(index)
  index[naDigit] = 1
  naDigitDim = dimension[naDigit]
  
  for(n in 0:(prod(dimension[naDigit])-1)) {
    for(digit in 1:length(naDigitDim)) {
      index[naDigit][digit] = (floor(n/prod(naDigitDim[0:(digit-1)])) %% naDigitDim[digit])+1
    }
    flatIndex = convertToFlatIndex(index,dimension)
    flatIndexVec = append(flatIndexVec,flatIndex)
  }
  
  return(flatIndexVec)
}

BnExactInference = function(tempFit,evidence,threshold=0.0005,kmprob=NULL){
  evidence$time = NULL
  evidence$TIMEPOINT = NULL
  evidence$PREVTIMEPOINT = NULL
  evidence$delta = NULL
  evidence$id = NULL
  
  cpt = tempFit[["TIMEPOINT"]][['prob']]
  p = colnames(as.data.frame(cpt))
  p = p[p!="Freq"]
  if(p=='Var1') {p=c('TIMEPOINT')}
  temp_evidence = evidence
  temp_evidence$TIMEPOINT = factor(0,levels = c('0','1'))
  ap1 = cpt[convertToFlatIndex(temp_evidence[p],dim(cpt))]
  #temp_evidence$TIMEPOINT = factor(1,levels = c('0','1'))
  #ap2 = cpt[convertToFlatIndex(temp_evidence[p],dim(cpt))]
  ap2 = 1-ap1
  p0 = 0
  p1 = 0
  variables = children(tempFit,'TIMEPOINT')
  for(v in variables) {
    cpt = tempFit[[v]][['prob']]
    p = colnames(as.data.frame(cpt))
    p = p[p!="Freq"]
    temp_evidence = evidence
    temp_evidence$TIMEPOINT = factor(0,levels = c('0','1'))
    t1 = cpt[convertToFlatIndex(temp_evidence[p],dim(cpt))]
    temp_evidence$TIMEPOINT = factor(1,levels = c('0','1'))
    t2 = cpt[convertToFlatIndex(temp_evidence[p],dim(cpt))]
    
    
    if(t1>threshold) {rawp0=t1} else {rawp0=threshold}
    if(t2>threshold) {rawp1=t2} else {rawp1=threshold}
    p0 = p0 + log(rawp0)
    p1 = p1 + log(rawp1)
    
  }
  #ap1=kmprob
  p0_Prob = exp(p0 + log(ap1))/(exp(p0 + log(ap1))+exp(p1 + log(1-ap1)))

  return(p0_Prob)
}

BnExactInferenceCumulated <- function(fitList,evidence) {
  evidence$time = NULL
  evidence$TIMEPOINT = NULL
  evidence$PREVTIMEPOINT = NULL
  evidence$delta = NULL
  evidence$id = NULL
  #testing$PREVTIMEPOINT = factor(integer(nrow(testing)),levels = c('0','1'))

  cumulatedProb = 1
  curve = rep(1,length(fitList))
  
  for(i in 1:length(fitList)) {
    tempFit = fitList[[i]]
    
    cumulatedProb = cumulatedProb * BnExactInference(tempFit,evidence)
    curve[i]=cumulatedProb
    
  }
  return(curve)
}

plotDag <-function(dag){
  #plot with bigger fontsize
  g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(dag))
  graph::nodeRenderInfo(g) <- list(fontsize=80)
  Rgraphviz::renderGraph(g)
}
