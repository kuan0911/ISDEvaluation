source("Models/cleanDataForBayesNet.R")
source("Models/bayesianNetHelper.R")

sythesize2 = function(n=1000) {
  
  
  timePoints = c(5, 13, 24, 32, 44, 57)
  m = length(timePoints)
  
  # data generation.
  LV3 = c(1, 2, 3)
  LorD = as.integer(c(0,1))
  
  a = sample(LV3, 5000, prob = c(0.4,0.4,0.2), replace = TRUE)
  b = sample(c(1, 2), 5000, prob = c(0.75, 0.25), replace = TRUE)
  c = sample(LV3, 5000, prob = c(0.75, 0.2, 0.05), replace = TRUE)
  d = sample(LV3, 5000, prob = c(0.4, 0.3, 0.3), replace = TRUE)
  e = sample(c(1, 2), 5000, prob = c(0.75, 0.25), replace = TRUE)
  g = sample(LV3, 5000, prob = c(0.4, 0.3, 0.3), replace = TRUE)
  
  TIMEPOINT = a
  TIMEPOINT[TIMEPOINT == 1] = sample(LorD, length(which(TIMEPOINT == 1)), prob = c(0.95, 0.05), replace = TRUE)
  TIMEPOINT[TIMEPOINT == 2] = sample(LorD, length(which(TIMEPOINT == 2)), prob = c(0.85, 0.15), replace = TRUE)
  TIMEPOINT[TIMEPOINT == 3] = sample(LorD, length(which(TIMEPOINT == 3)), prob = c(0.7, 0.3), replace = TRUE)
  
  # b = TIMEPOINT
  # b[b == 0] = sample(LV3, length(which(b == 0)), prob = c(0.6, 0.2, 0.2), replace = TRUE)
  # b[b == 1] = sample(LV3, length(which(b == 1)), prob = c(0.15, 0.15, 0.7), replace = TRUE)
  
  # g = TIMEPOINT
  # g[g == 0] = sample(LV3, length(which(g == 0)), prob = c(0.5, 0.4, 0.1), replace = TRUE)
  # g[g == 1] = sample(LV3, length(which(g == 1)), prob = c(0.1, 0.4, 0.5), replace = TRUE)
  
  syndata1 = data.frame(
    TIMEPOINT = factor(TIMEPOINT, levels = LorD),
    A = factor(a, levels = LV3),
    B = factor(b, levels = c(1, 2)),
    C = factor(c, levels = LV3),
    D = factor(d, levels = LV3),
    E = factor(e, levels = c(1, 2)),
    G = factor(g, levels = LV3)
  )
  dag1 = empty.graph(c('TIMEPOINT','A','B','C','D','E','G'))
  #arc.set1 = matrix(c("A", "TIMEPOINT", "TIMEPOINT", "B", "TIMEPOINT", "G"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
  arc.set1 = matrix(c("A", "TIMEPOINT"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
  arcs(dag1) = arc.set1
  fit1 = bn.fit(dag1,syndata1)
  #graphviz.plot(fit1)
  
  TIMEPOINT = d
  TIMEPOINT[TIMEPOINT == 1] = sample(LorD, length(which(TIMEPOINT == 1)), prob = c(0.95, 0.05), replace = TRUE)
  TIMEPOINT[TIMEPOINT == 2] = sample(LorD, length(which(TIMEPOINT == 2)), prob = c(0.85, 0.15), replace = TRUE)
  TIMEPOINT[TIMEPOINT == 3] = sample(LorD, length(which(TIMEPOINT == 3)), prob = c(0.7, 0.3), replace = TRUE)
  
  # e = TIMEPOINT
  # e[e == 0] = sample(c(1, 2), length(which(e == 0)), prob = c(0.7, 0.3), replace = TRUE)
  # e[e == 1] = sample(c(1, 2), length(which(e == 1)), prob = c(0.3, 0.7), replace = TRUE)
  
  # g = TIMEPOINT
  # g[g == 0] = sample(LV3, length(which(g == 0)), prob = c(0.5, 0.4, 0.1), replace = TRUE)
  # g[g == 1] = sample(LV3, length(which(g == 1)), prob = c(0.1, 0.4, 0.5), replace = TRUE)
  
  syndata2 = data.frame(
    TIMEPOINT = factor(TIMEPOINT, levels = LorD),
    A = factor(a, levels = LV3),
    B = factor(b, levels = c(1, 2)),
    C = factor(c, levels = LV3),
    D = factor(d, levels = LV3),
    E = factor(e, levels = c(1, 2)),
    G = factor(g, levels = LV3)
  )
  
  dag2 = empty.graph(c('TIMEPOINT','A','B','C','D','E','G'))
  #arc.set2 = matrix(c("D", "TIMEPOINT", "TIMEPOINT", "E", "TIMEPOINT", "G"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
  arc.set2 = matrix(c("D", "TIMEPOINT"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
  arcs(dag2) = arc.set2
  fit2 = bn.fit(dag2,syndata2)
  #graphviz.plot(fit2)
  
  TIMEPOINT = g
  TIMEPOINT[TIMEPOINT == 1] = sample(LorD, length(which(TIMEPOINT == 1)), prob = c(0.8, 0.2), replace = TRUE)
  TIMEPOINT[TIMEPOINT == 2] = sample(LorD, length(which(TIMEPOINT == 2)), prob = c(0.7, 0.3), replace = TRUE)
  TIMEPOINT[TIMEPOINT == 3] = sample(LorD, length(which(TIMEPOINT == 3)), prob = c(0.6, 0.4), replace = TRUE)
  
  syndata3 = data.frame(
    TIMEPOINT = factor(TIMEPOINT, levels = LorD),
    A = factor(a, levels = LV3),
    B = factor(b, levels = c(1, 2)),
    C = factor(c, levels = LV3),
    D = factor(d, levels = LV3),
    E = factor(e, levels = c(1, 2)),
    G = factor(g, levels = LV3)
  )
  dag3 = empty.graph(c('TIMEPOINT','A','B','C','D','E','G'))
  arc.set3 = matrix(c("G", "TIMEPOINT"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
  arcs(dag3) = arc.set3
  fit3 = bn.fit(dag3,syndata3)
  #graphviz.plot(fit3)
  
  fitList <- vector("list", length(timePoints))
  
  fitList[[1]] = fit1
  fitList[[2]] = fit1
  fitList[[3]] = fit2
  fitList[[4]] = fit2
  fitList[[5]] = fit3
  fitList[[6]] = fit3
  #fitList[[7]] = fit2
  #fitList[[8]] = fit3
  
  a = sample(LV3, 5000, prob = c(0.4,0.4,0.2), replace = TRUE)
  b = sample(c(1, 2), 5000, prob = c(0.5, 0.5), replace = TRUE)
  c = sample(LV3, 5000, prob = c(0.75, 0.2, 0.05), replace = TRUE)
  d = sample(LV3, 5000, prob = c(0.4, 0.3, 0.3), replace = TRUE)
  e = sample(c(1, 2), 5000, prob = c(0.75, 0.25), replace = TRUE)
  g = sample(LV3, 5000, prob = c(0.3, 0.3, 0.4), replace = TRUE)
  
  covariateData = data.frame(
    A = factor(a, levels = LV3),
    B = factor(b, levels = c(1, 2)),
    C = factor(c, levels = LV3),
    D = factor(d, levels = LV3),
    E = factor(e, levels = c(1, 2)),
    G = factor(g, levels = LV3)
  )
  
  covariateDag = empty.graph(c('A','B','C','D','E','G'))
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
        #predicted = predict(object = fitList[[i]], node = "TIMEPOINT", data = covariate, method = "bayes-lw", n=2000, prob = TRUE)
        prob = BnExactInference(fitList[[i]],covariate,kmprob=NULL,noise=F,lev=list(LV3,c(1, 2),LV3,LV3,c(1, 2),LV3))
        if(runif(1,0,1)>prob) {
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
    if(runif(1,0,1)<0.4) {
      censoreTime = runif(1,0.1,timePoints[m])
      if(censoreTime<syntheticData[k,'time']) {
        syntheticData[k,'delta'] = as.integer(0)
        syntheticData[k,'time'] = censoreTime
      }
    }
  }
  
  write.csv(syntheticData,'Data/syntheticData.csv')
  
}