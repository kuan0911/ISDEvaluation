survivalDataset <- read.csv(file = "Data/All_Data_updated_may2011_CLEANED.csv")
names(survivalDataset)[c(1,2)] = c("time", "delta")
survivalDataset$delta = 1 - survivalDataset$delta 

survivalDataset <- read.csv(file = "Data/NACD_colorectal.csv")
names(survivalDataset)[c(1,2)] = c("time", "delta")
survivalDataset$delta = 1 - survivalDataset$delta 

survivalDataset <- read.csv(file = "Data/NACD2.csv")
names(survivalDataset)[c(1,2)] = c("time", "delta")
survivalDataset$delta = 1 - survivalDataset$delta 

survivalDataset <- read.csv(file = "Data/DBCD.csv")
names(survivalDataset)[c(1,2)] = c("time", "delta")
survivalDataset$delta = 1 - survivalDataset$delta 

survivalDataset <- read.csv(file = "Data/NSBCD.csv")
names(survivalDataset)[c(1,2)] = c("time", "delta")
survivalDataset$delta = 1 - survivalDataset$delta 

survivalDataset = survival::lung
names(survivalDataset)[c(3)] = c("delta")
survivalDataset$delta = survivalDataset$delta - 1

survivalDataset = read.csv(file = "http://pssp.srv.ualberta.ca/system/predictors/datasets/000/000/076/original/All_Data_updated_may2011_Stage4_Stomach.csv?1354708504")
names(survivalDataset)[c(1,2)] = c("time", "delta")
survivalDataset$delta = 1 - survivalDataset$delta 

survivalDataset = loadCOADREAD()

rm(list = ls())

source('analysisMaster.R')

extra = survivalDataset[c(2,4,6,13,14),]
survivalDataset = survivalDataset[survivalDataset$delta==1,]
survivalDataset = rbind(survivalDataset,extra)

C1Vec = c(0.1,0.5,1,1.5,2,3,4,5,6,7,8,9,10,14)
histogramList <- c()
DCalList = c()
ConcList = c()
L1List = c()
curveList = c()
resultList = c()
for(C1 in C1Vec) {
  print('Start collect data for')
  print(C1)
  ISD = analysisMaster(survivalDataset,BayesianC1=C1, MTLRModel=F,BayesianNetModel=T,KaplanMeier = F, FS = T, numberOfFolds = 5)
  histogramList = c(histogramList,ISD$DcalHistogram)
  DCalList = c(DCalList,ISD$results$DCalibration[1])
  ConcList = c(ConcList,mean(ISD$results$Concordance))
  L1List = c(L1List,mean(ISD$results$L1Loss))
  curveList = c(curveList,ISD$survivalCurves)
  resultList = c(resultList,ISD$results)
  print(DCalList)
  print(ConcList)
  print(L1List)
}
ISD = analysisMaster(survivalDataset, MTLRModel=F,BayesianNetModel=T,KaplanMeier = F, FS = T, numberOfFolds = 5)

ISD = analysisMaster(survivalDataset, CoxKP = T,CoxKPEN = F, KaplanMeier = T, RSFModel = F, AFTModel = T, MTLRModel = T, BayesianNetModel = F, numberOfFolds = 5)

curve = ISD$survivalCurves$Bayes
for(j in 2:ncol(curve)) {
  for(i in 1:(nrow(curve)-1)) {
    if(curve[i+1,j] > curve[i,j]){
      #print('fix curve')
      curve[i+1,j] = curve[i,j]
    }
  }
}

plotSurvivalCurves(curve, 1:10)

plotSurvivalCurves(ISD$survivalCurves$Bayes, 1:10)

barplot(ISD$DcalHistogram$BayesianNet,horiz=TRUE)

#into detailed step
imputeZero=T
FS = T
verbose = T
numberOfFolds =2
i = 1



loadCOADREAD = function() {
  survivalDataset <- read.csv(file = "Data/COADREAD.csv")
  survivalDataset$bcr_patient_barcode = NULL
  survivalDataset$bcr_patient_uuid = NULL
  survivalDataset$patient_id = NULL
  survivalDataset = survivalDataset[!is.na(survivalDataset$days_to_death) | !is.na(survivalDataset$days_to_last_followup),]
  survivalDataset$time = survivalDataset$days_to_death
  survivalDataset$time[is.na(survivalDataset$time)] = survivalDataset$days_to_last_followup[!is.na(survivalDataset$days_to_last_followup)]
  survivalDataset$days_to_death = NULL
  survivalDataset$days_to_last_followup = NULL
  survivalDataset$vital_status = replace(survivalDataset$vital_status, survivalDataset$vital_status=='alive', 0)
  survivalDataset$vital_status = replace(survivalDataset$vital_status, survivalDataset$vital_status=='dead', 1)
  survivalDataset$delta = survivalDataset$vital_status
  survivalDataset$vital_status = NULL
  for(c in colnames(survivalDataset)) {
    if(sum(is.na(survivalDataset[[c]]))>nrow(survivalDataset)*0.7){
      survivalDataset[[c]] = NULL
    }
    if(length(unique(survivalDataset[[c]]))==1){
      survivalDataset[[c]] = NULL
    }
    if(!is.numeric(survivalDataset[[c]]) & length(unique(survivalDataset[[c]]))>nrow(survivalDataset)*0.7) {
      survivalDataset[[c]] = NULL
    }
  }
  
  survivalDataset$bcr = NULL
  
  return(survivalDataset)
}
