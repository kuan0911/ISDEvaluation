survivalDataset <- read.csv(file = "Data/All_Data_updated_may2011_CLEANED.csv")
# names(survivalDataset)[c(1,2)] = c("time", "delta")
# survivalDataset$delta = 1 - survivalDataset$delta 

#survivalDataset <- read.csv(file = "Data/NACD_colorectal.csv")
names(survivalDataset)[c(1,2)] = c("time", "delta")
survivalDataset$delta = 1 - survivalDataset$delta 

survivalDataset <- read.csv(file = "Data/NSBCD.csv")
survivalDataset <- read.csv(file = "Data/EffectiveISD_data/GBM.csv")

survivalDataset <- read.csv(file = "Data/apcari/Neeraj_filtered_APCaRI_biopsydiagnosedcases_14Dec2020.csv")
survivalDataset = survivalDataset[survivalDataset$Biopsy.Date!=''&survivalDataset$Last.Follow.Up.Date!='',]
survivalDataset$time = as.Date(survivalDataset$Last.Follow.Up.Date,format="%m/%d/%Y")-as.Date(survivalDataset$Biopsy.Date,format="%m/%d/%Y")
survivalDataset$delta = as.integer(survivalDataset$Date.of.Death!="")
survivalDataset$weight = survivalDataset$weight..kgs.
survivalDataset$height = survivalDataset$height..cm.

for(k in 1:nrow(survivalDataset)) {
  if(survivalDataset[k,'delta']==1) {
    survivalDataset[k,'time'] = as.Date(survivalDataset[k,'Date.of.Death'],format="%m/%d/%Y")-as.Date(survivalDataset[k,'Biopsy.Date'],format="%m/%d/%Y")
  }
  if(survivalDataset[k,'time']<0) {
    survivalDataset[k,'time'] = 0
  }
  if(is.na(survivalDataset[k,'weight'])&!is.na(survivalDataset[k,'weight..lbs.'])) {
    survivalDataset[k,'weight'] = survivalDataset[k,'weight..lbs.']
  }
  if(is.na(survivalDataset[k,'height'])&!is.na(survivalDataset[k,'height..inches.'])) {
    survivalDataset[k,'height'] = survivalDataset[k,'height..inches.']
  }
}




source('Data/synthesizeData.R')
source('Data/synthesizeData2.R')
sythesize2(n=1000)
survivalDataset <- read.csv(file = "Data/syntheticData.csv")
survivalDataset$X = NULL 
#lapply(survivalDataset, as.double)

source('analysisMaster.R')

ISD = analysisMaster(survivalDataset, CoxKP = F, MTLRModel=F,BayesianNetModel=T,KaplanMeier = F,RSFModel=F, FS = T, numberOfFolds = 5)

ISD = analysisMaster(survivalDataset, CoxKP = F, MTLRModel=T,BayesianNetModel=F,KaplanMeier = F,RSFModel=F, FS = F, numberOfFolds = 5, foldIndex = foldIndex)

res = aveMetrics(ISD)

plotSurvivalCurves(ISD$survivalCurves$Bayes, 1:10)

binnames = c("[0.9,1]","[0.8,0.9)","[0.7,0.8)","[0.6,0.7)","[0.5,0.6)","[0.4,0.5)","[0.3,0.4)","[0.2,0.3)","[0.1,0.2)","[0,0.1)")
barplot(ISD$DcalHistogram$BayesianNet/sum(ISD$DcalHistogram$BayesianNet),horiz=TRUE,xlab='Proportion in bin',names.arg=binnames,las=1,main='')

barplot(100*ISD$DcalHistogram$Cox/sum(ISD$DcalHistogram$Cox),horiz=TRUE,xlab='Percentage in bin',names.arg=binnames,las=1,main='CoxKP',cex.axis=1.6,cex.names=1.6,cex.lab=1.6,cex.main=2)

delta0curveX = delta0curve$time
delta0curveY = delta0curve[,2]
delta1curveX = delta1curve$time
delta1curveY = delta1curve[,2]

imputeZero=T
FS = T
verbose = T
numberOfFolds =5
i = 1

aveMetrics = function(ISD) {
  print(paste(ISD$results$Model[1],' N =',ISD$results$N[1],sep=' '))
  C = mean(ISD$results$Concordance)
  Cstd = sd(ISD$results$Concordance)
  print(paste('Conc:',trunc(C*1000)/1000,'+/-',trunc(Cstd*1000)/1000,sep=' '))
  B = mean(ISD$results$BrierInt)
  Bstd = sd(ISD$results$BrierInt)
  print(paste('Brier:',trunc(B*1000)/1000,'+/-',trunc(Bstd*1000)/1000,sep=' '))
  L = mean(ISD$results$L1Loss)
  Lstd = sd(ISD$results$L1Loss)
  print(paste('L1Loss:',trunc(L*1000)/1000,'+/-',trunc(Lstd*1000)/1000,sep=' '))
  D = ISD$results$DCalibration[1]
  print(paste('DCal:',trunc(D*1000)/1000,sep=' '))
  return(list(Concordance=C,Concordance_std=Cstd,BrierInt=B,BrierInt_std=Bstd,L1Loss=L,L1Loss_std=Lstd,DCalibration=D))
}



ISD = analysisMaster(survivalDataset, CoxKP = T, MTLRModel=F,BayesianNetModel=F,KaplanMeier = F,RSFModel=F, FS = F, numberOfFolds = 5,verbose = F)
res = aveMetrics(ISD)
ISD = analysisMaster(survivalDataset, CoxKP = F, MTLRModel=T,BayesianNetModel=F,KaplanMeier = F,RSFModel=F, FS = F, numberOfFolds = 5,verbose = F)
res = aveMetrics(ISD)
ISD = analysisMaster(survivalDataset, CoxKP = F, MTLRModel=F,BayesianNetModel=F,KaplanMeier = F,RSFModel=T, FS = F, numberOfFolds = 5,verbose = F)
res = aveMetrics(ISD)


