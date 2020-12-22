# survivalDataset <- read.csv(file = "Data/All_Data_updated_may2011_CLEANED.csv")
# names(survivalDataset)[c(1,2)] = c("time", "delta")
# survivalDataset$delta = 1 - survivalDataset$delta 

survivalDataset <- read.csv(file = "Data/NACD_colorectal.csv")
names(survivalDataset)[c(1,2)] = c("time", "delta")
survivalDataset$delta = 1 - survivalDataset$delta 

survivalDataset <- read.csv(file = "Data/KIPAN.csv")

# for(k in 1:nrow(survivalDataset)) {
#   if(runif(1,0,1)<0.7 & survivalDataset[k,'delta']==1) {
#     censoreTime = runif(1,0.1,max(survivalDataset$time))
#     if(censoreTime<survivalDataset[k,'time']) {
#       survivalDataset[k,'delta'] = as.integer(0)
#       survivalDataset[k,'time'] = censoreTime
#     }
#   }
# }
print(sum(survivalDataset$delta==0)/nrow(survivalDataset))
# survivalDataset$STAGE_4 = NULL
# survivalDataset$STAGE_3 = NULL
# survivalDataset$STAGE_2 = NULL
# survivalDataset$STAGE_1 = NULL
# survivalDataset$PERFORMANCE_STATUS_4 = NULL
# survivalDataset$PERFORMANCE_STATUS_3 = NULL
# survivalDataset$PERFORMANCE_STATUS_2 = NULL
# survivalDataset$PERFORMANCE_STATUS_1 = NULL
# survivalDataset$PERFORMANCE_STATUS_0 = NULL
# survivalDataset$AGE65 = NULL

survivalDataset <- read.csv(file = "Data/covid_hospitalized_data.csv")
survivalDataset$delta = as.integer(as.logical(survivalDataset$event))
survivalDataset$event = NULL
survivalDataset$X = NULL

survivalDataset <- read.csv(file = "Data/LifeExpectancyData.csv")

source('Data/synthesizeData.R')
source('Data/synthesizeData2.R')
sythesize2(n=3000)
survivalDataset <- read.csv(file = "Data/syntheticData.csv")
survivalDataset$X = NULL
#lapply(survivalDataset, as.double)

source('analysisMaster.R')

ISD = analysisMaster(survivalDataset, CoxKP = F, MTLRModel=T,BayesianNetModel=F,KaplanMeier = F, FS = T, numberOfFolds = 5)

ISD = analysisMaster(survivalDataset, MTLRModel=T,BayesianNetModel=T,KaplanMeier = F, FS = F, numberOfFolds = 5)

ISD = analysisMaster(survivalDataset, CoxKP = T,CoxKPEN = F, KaplanMeier = F, RSFModel = F, AFTModel = F, MTLRModel = F, BayesianNetModel = F, FS=F, numberOfFolds = 5)

plotSurvivalCurves(ISD$survivalCurves$Bayes, 1:10)

barplot(ISD$DcalHistogram$BayesianNet,horiz=TRUE)


delta0curveX = delta0curve$time
delta0curveY = delta0curve[,2]
delta1curveX = delta1curve$time
delta1curveY = delta1curve[,2]

imputeZero=T
FS = T
verbose = T
numberOfFolds =2
i = 1


