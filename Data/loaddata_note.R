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

rm(list = ls())

source('analysisMaster.R')

extra = survivalDataset[c(2,4,6,13,14),]
survivalDataset = survivalDataset[survivalDataset$delta==1,]
survivalDataset = rbind(survivalDataset,extra)
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

#into detailed step
imputeZero=T
FS = T
verbose = T
numberOfFolds =2
i = 1

