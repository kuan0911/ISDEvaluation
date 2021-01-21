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

BnExactInference = function(tempFit,evidence,kmprob=NULL,noise=F,lev=NULL,noiserate = 0.1){
  #threshold=0.0005
  evidence[c('PREVTIMEPOINT','TIMEPOINT','time','delta','id')] = NULL
  
  cpt = tempFit[["TIMEPOINT"]][['prob']]
  p = colnames(as.data.frame(cpt))
  p = p[p!="Freq"]
  if(length(p)==1){
    if(p=='Var1') {p=c('TIMEPOINT')}
  }
  if(noise) {
    if(runif(1,0,1)<noiserate) {
      index = sample(1:length(evidence), 1)
      evidence[index] = sample(lev[[index]],1)
    }
  }
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
    temp_evidence$TIMEPOINT = factor(0,levels = c(0,1))
    t1 = cpt[convertToFlatIndex(temp_evidence[p],dim(cpt))]
    temp_evidence$TIMEPOINT = factor(1,levels = c(0,1))
    t2 = cpt[convertToFlatIndex(temp_evidence[p],dim(cpt))]
    if(is.na(t2)|is.na(t2)) {'BnExactInference: t1 or t2 is NA'}
    
    if(t1>0.0005) {rawp0=t1} else {rawp0=0.0005}
    if(t2>0.0005) {rawp1=t2} else {rawp1=0.0005}
    p0 = p0 + log(rawp0)
    p1 = p1 + log(rawp1)
    
  }
  #ap1=kmprob
  p0_Prob = exp(p0 + log(ap1))/(exp(p0 + log(ap1))+exp(p1 + log(1-ap1)))

  return(p0_Prob)
}

BnExactInferenceCumulated <- function(fitList,evidenceList,currentTime) {
  if(currentTime==0){return(1)}
  #testing$PREVTIMEPOINT = factor(integer(nrow(testing)),levels = c('0','1'))
  cumulatedProb = 1
  curve = rep(1,currentTime)

  for(i in 1:currentTime) {
    tempFit = fitList[[i]]
    evidence = evidenceList[[i]]
    if(is.null(tempFit)|is.null(evidence)) {print('BnExactInferenceCumulated: error. lenth not match')}
    evidence[c('PREVTIMEPOINT','TIMEPOINT','time','delta','id')] = NULL
    cumulatedProb = cumulatedProb * BnExactInference(tempFit,evidence)
    curve[i] = cumulatedProb
  }
  return(curve[currentTime])
}

plotDag <-function(dag){
  #plot with bigger fontsize
  g <- Rgraphviz::layoutGraph(bnlearn::as.graphNEL(dag))
  graph::nodeRenderInfo(g) <- list(fontsize=80)
  Rgraphviz::renderGraph(g)
}

timepointFixer <- function(fit,data,weight=NULL,iss=5) {
  if(!is.null(data$PREVTIMEPOINT)) {print('Warning : timepointFixer data include missing data in PREVTIMEPOINT')}
  if(!is.null(weight)) {if(anyNA(data)) {print('Warning : timepointFixer include missing data')}}
  data[,c('time','delta','id')] = NULL
  dataComplete = data[complete.cases(data),]
  dataComplete$PREVTIMEPOINT = NULL
  
  cpt = coef(fit$TIMEPOINT)
  p = colnames(as.data.frame(cpt))
  p = p[p!="Freq"]
  if(p=='Var1') {p=c('TIMEPOINT')}
  newcpt = table(dataComplete[,p])
  
  if(!is.null(weight)) {
    for(k in 1:nrow(data)) {
      index = data[k,p]
      flatIndex = convertToFlatIndex(index,dim(newcpt))
      newcpt[flatIndex] = newcpt[flatIndex] - (1-weight[k])
    }
  }
  
  multipleFlateIndex = convertToFlatIndexMultiple(c(1,rep(NA,length(p)-1)),dim(newcpt))
  timepoint1 = newcpt[multipleFlateIndex]
  multipleFlateIndex = convertToFlatIndexMultiple(c(2,rep(NA,length(p)-1)),dim(newcpt))
  timepoint2 = newcpt[multipleFlateIndex]
  
  timepointNACount = timepoint1
  for(k in 1:length(timepointNACount)){timepointNACount[k]=0}
  
  for(k in 1:nrow(data)) {
    if(is.na(data[k,'TIMEPOINT'])) {
      index = data[k,p[-1]]
      if(length(index)==0) {
        timepointNACount = timepointNACount + 1
      }else {
        flatIndex = convertToFlatIndex(index,dim(timepointNACount))
        timepointNACount[flatIndex] = timepointNACount[flatIndex] + 1
      }
    }
  }
  
  alive = sum(data$TIMEPOINT==0,na.rm=T)
  dead = sum(data$TIMEPOINT==1,na.rm=T)
  censored = sum(is.na(data$TIMEPOINT))
  iss_dead = iss * (dead/(alive+dead+censored))
    
  timepointProb = timepoint1
  for(k in 1:length(timepointProb)){timepointProb[k]=0}
    
  for(j in 1:length(timepointProb)) {
    if(length(timepoint2)>1) {
      survivalRate = 1 - (timepoint2[j]+iss_dead)/(timepoint1[j]+timepoint2[j]+0.5*timepointNACount[j]+iss)
    }else {
      survivalRate = 1 - (timepoint2+iss_dead)/(timepoint1+timepoint2+0.5*timepointNACount+iss)
    }
    if(is.nan(survivalRate)) {
      survivalRate = 1-(iss_dead/iss)
      print('divided by zero')
    }
    timepointProb[j] = survivalRate
  }

  multipleFlateIndex = convertToFlatIndexMultiple(c(1,rep(NA,length(p)-1)),dim(newcpt))
  newcpt[multipleFlateIndex] = timepointProb
  multipleFlateIndex = convertToFlatIndexMultiple(c(2,rep(NA,length(p)-1)),dim(newcpt))
  newcpt[multipleFlateIndex] = 1 - timepointProb
  
  fit$TIMEPOINT = newcpt
  return(fit)
}

childrenFixer <- function(fit,data,weight=NULL,iss=1) {
  if(!is.null(data$PREVTIMEPOINT)) {print('Warning : timepointFixer data include missing data in PREVTIMEPOINT')}
  if(!is.null(weight)) {if(anyNA(data)) {print('Warning : timepointFixer include missing data')}}
  data[,c('time','delta','id')] = NULL
  dataComplete = data[complete.cases(data),]
  dataComplete$PREVTIMEPOINT = NULL
  
  #print('Weighted fitting')
  for(node in children(fit,'TIMEPOINT')) {
    cpt = coef(fit[[node]])
    p = colnames(as.data.frame(cpt))
    p = p[p!="Freq"]
    newcpt = cpt
    for(k in 1:length(newcpt)){newcpt[k]=0}
    if(is.null(weight)) {weight = rep(1,nrow(dataComplete))}
    
    for(k in 1:nrow(dataComplete)) {
      index = dataComplete[k,p]
      flatIndex = convertToFlatIndex(index,dim(newcpt))
      newcpt[flatIndex] = newcpt[flatIndex] + weight[k]
      
    }

    newcpt = newcpt + iss
    
    normal = as.table(array(integer(prod(dim(newcpt)[-1])),dim=dim(newcpt)[-1]))
    for(z in 1:length(normal)){normal[z]=0}
    for(j in 1:dim(newcpt)[1]) {
      multipleFlateIndex = convertToFlatIndexMultiple(c(j,rep(NA,length(p)-1)),dim(newcpt))
      normal = normal + newcpt[multipleFlateIndex]
    }
    
    for(j in 1:dim(newcpt)[1]) {
      multipleFlateIndex = convertToFlatIndexMultiple(c(j,rep(NA,length(p)-1)),dim(newcpt))
      newcpt[multipleFlateIndex] = newcpt[multipleFlateIndex]/normal
    }
    
    sumOther = as.table(array(integer(prod(dim(newcpt)[-1])),dim=dim(newcpt)[-1]))
    for(z in 1:length(sumOther)){sumOther[z]=0}
    for(j in 1:dim(newcpt)[1]) {
      if(j != dim(newcpt)[1]) {
        multipleFlateIndex = convertToFlatIndexMultiple(c(j,rep(NA,length(p)-1)),dim(newcpt))
        sumOther = sumOther + newcpt[multipleFlateIndex]
      }else if(j == dim(newcpt)[1]) {
        multipleFlateIndex = convertToFlatIndexMultiple(c(j,rep(NA,length(p)-1)),dim(newcpt))
        newcpt[multipleFlateIndex] = 1 - sumOther
      }
    }
    
    fit[[node]] = newcpt
  }
  return(fit)
}

customScoreFunction = function(node, parents, data, args) {
  if(!is.null(args$weight)) {weight = args$weight}
  else {weight = rep(1,nrow(data))}
  kpenalty = (log(sum(weight))/2)
  #kpenalty = 1
  df = as.data.frame(data[,node])
  if(length(parents)==0){
    t = table(df)
    for(k in 1:nrow(df)) {
      index = as.factor(df[k,node])
      t[index] = t[index] - (1-weight)
    }
    td = table(df)/nrow(df)
    totalloglike = sum(log(td)*t) - kpenalty*(length(t)-1)
  }else {
    e = empty.graph(c(node,parents))
    m = matrix(parents,nrow=length(parents),ncol=2)
    m[,2] = node
    arcs(e) = m
    df = as.data.frame(data[,c(node,parents),drop=FALSE])
    fit = bn.fit(e,df)
    totalloglike = 0
    cpt = coef(fit[[node]])
    p = colnames(as.data.frame(cpt))
    p = p[p!="Freq"]
    indexarray = df[,p]
    for(k in 1:nrow(indexarray)) {
      flatIndex = convertToFlatIndex(indexarray[k,],dim(cpt))
      loglike = log(cpt[flatIndex])*weight[k]
      #evidence = as.data.frame(data[k,parents,drop=FALSE])
      #predicted = predict(fit, node=node, data=evidence,method = "parents", prob = TRUE)
      #index = as.factor(df[k,node])
      totalloglike = totalloglike + loglike
    }
    parentparams = 0
    for(p in parents) {
      parentparams = parentparams + length(levels(df[,p]))-1
    }
    parapenalty = nparams(fit)-parentparams
    totalloglike = totalloglike - kpenalty*parapenalty
  }
  if(!is.null(args$prevFit)) {
    prevParents = parents(args$prevFit,node)
    totalloglike = totalloglike - kpenalty*4*length(setdiff(parents,prevParents))
  }
  if(!is.null(args$nextFit)) {
    nextParents = parents(args$nextFit,node)
    totalloglike = totalloglike - kpenalty*4*length(setdiff(parents,nextParents))
  }
  return(totalloglike)
}

weightedImpute = function(fitList,fit,dataListNAadapt,inputdata,currentTime) {
  if(is.null(inputdata$PREVTIMEPOINT)) {
    print('PREVTIMEPOINT missing when imputing weight')
    inputdata$PREVTIMEPOINT = rep(0,nrow(inputdata))
  }
  imputeCount = 0
  outputdata = inputdata
  weight = rep(1,nrow(outputdata))
  for(k in 1:nrow(inputdata)) {
    if(is.na(outputdata[k,'TIMEPOINT'])) {
      imputeCount = imputeCount+1
      datainstance = inputdata[k,]
      evidenceList = getEvidenceList(dataListNAadapt,datainstance$id,currentTime)
      if(!is.na(inputdata[k,'PREVTIMEPOINT'])) {
        probSurvive = 1
        fitList[[currentTime]] = fit
        prob = BnExactInferenceCumulated(fitList,evidenceList,currentTime)
      }else {
        fitList[[currentTime]] = fit
        probSurvive = BnExactInferenceCumulated(fitList,evidenceList,currentTime-1)
        prob = 1-(BnExactInferenceCumulated(fitList,evidenceList,currentTime-1)-BnExactInferenceCumulated(fitList,evidenceList,currentTime))
        
      }
      #prob = BnExactInference(fit,evidenceList[[currentTime]])
      fitList[[currentTime]] = fit
      #prob = 1-(BnExactInferenceCumulated(fitList,evidenceList,currentTime-1)-BnExactInferenceCumulated(fitList,evidenceList,currentTime))
      if(is.nan(probSurvive) | is.nan(prob)) {print('Error imputing: prob is nan')}
      if(prob>1 |prob<0){print('weightedImpute: prob error')}
      
      outputdata[k,'TIMEPOINT'] = 0
      weight[k] = probSurvive*prob
      
      datainstance['TIMEPOINT'] = 1
      outputdata = rbind(outputdata,datainstance)
      weight = c(weight,probSurvive*(1-prob))
    }
  }
  #print(imputeCount)
  outputdata[,c('PREVTIMEPOINT','time','delta','id')] = NULL
  if(anyNA(outputdata)) {print('Warning: data imputed contain NA value')}
  return(list(data=outputdata,weight=weight))
}

randomflip = function(data,rate=0.4) {
  if(anyNA(data)){print('randomflip: Error. missing data')}
  for(k in 1:nrow(data)) {
    if(runif(1,0,1)<rate) {
      if(data[k,'TIMEPOINT']==1){
        data[k,'TIMEPOINT']=0
      }else {
        data[k,'TIMEPOINT']=1
      }
      if(is.na(data[k,'TIMEPOINT'])){print('randomflip: Error, produce missing data')}
    }
  }
  return(data)
}
getEvidenceList = function(dataListNAadapt,id,currentTime) {
  evidenceList = vector("list", currentTime)
  for(i in 1:currentTime) {
    data = dataListNAadapt[[i]]
    evidence = data[data$id==id,]
    evidence[c('PREVTIMEPOINT','TIMEPOINT','time','delta','id')] = NULL
    evidenceList[[i]] = evidence
  }
  return(evidenceList)
}




  
  
  
  
  
  