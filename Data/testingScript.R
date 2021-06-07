survivalDataset <- read.csv(file = "Data/All_Data_updated_may2011_CLEANED.csv")
# names(survivalDataset)[c(1,2)] = c("time", "delta")
# survivalDataset$delta = 1 - survivalDataset$delta 

#survivalDataset <- read.csv(file = "Data/NACD_colorectal.csv")
names(survivalDataset)[c(1,2)] = c("time", "delta")
survivalDataset$delta = 1 - survivalDataset$delta 

survivalDataset <- read.csv(file = "Data/NSBCD.csv")
survivalDataset <- read.csv(file = "Data/EffectiveISD_data/GBM.csv")

survivalDataset <- read.csv(file = "Data/covid_cv_folds/discharge_wo_argentina.csv")
#survivalDataset2 <- read.csv(file = "Data/covid_cv_folds/asian_discharge_city.csv")

# tempEvent = survivalDataset[survivalDataset$event!=0&survivalDataset$event!=1,'time']
# survivalDataset[survivalDataset$event!=0&survivalDataset$event!=1,'time'] = survivalDataset[survivalDataset$event!=0&survivalDataset$event!=1,'event']
# survivalDataset[survivalDataset$event!=0&survivalDataset$event!=1,'event'] = tempEvent

survivalDataset$chronic_disease=NULL
survivalDataset$travel_hist_date=NULL
survivalDataset$travel_hist_location=NULL

survivalDataset$chronic_disease_binary=NULL
survivalDataset$latitude=NULL
survivalDataset$longitude=NULL

survivalDataset$countries = paste(as.character(as.integer(survivalDataset$latitude/3)),as.character(as.integer(survivalDataset$longitude/3)))

survivalDataset$GDP_per_capita_country=NULL
survivalDataset$GDP_total_country=NULL

survivalDataset$population_density_city=NULL
survivalDataset$population_density_country=NULL

survivalDataset$population_density_city=as.numeric(gsub(",", "", survivalDataset$population_density_city, fixed = TRUE))
survivalDataset$population_density_country=as.numeric(gsub(",", "", survivalDataset$population_density_country, fixed = TRUE))


fold3 = read.csv(file = "Data/covid_cv_folds/fold_0.csv")
fold2 = read.csv(file = "Data/covid_cv_folds/fold_1.csv")
fold1 = read.csv(file = "Data/covid_cv_folds/fold_2.csv")
fold4 = read.csv(file = "Data/covid_cv_folds/fold_3.csv")
fold5 = read.csv(file = "Data/covid_cv_folds/fold_4.csv")
n1=nrow(fold1);n2=nrow(fold2);n3=nrow(fold3);n4=nrow(fold4);n5=nrow(fold5)

survivalDataset = do.call("rbind", list(fold1,fold2,fold3,fold4,fold5))
survivalDataset$delta = survivalDataset$event
survivalDataset$event = NULL
foldIndex = list(c(1:n1),c(n1+1:n2),c(n1+n2+1:n3),c(n1+n2+n3+1:n4),c(n1+n2+n3+n4+1:n5))

survivalDataset[survivalDataset=="False"]=0
survivalDataset[survivalDataset=="True"]=1
survivalDataset$delta = as.integer(survivalDataset$delta)

survivalDataset$translon = cos(survivalDataset$latitude)*cos(survivalDataset$longitude)
survivalDataset$translat = cos(survivalDataset$latitude)*sin(survivalDataset$longitude)

new_df <- survivalDataset
new_df$city <- factor(new_df$city, exclude = NULL)
#new_df$countries <- factor(new_df$countries, exclude = NULL)
new_df <- model.matrix(~.-1, data = new_df[c("city")],
                       contrasts.arg = list(
                         city = contrasts(new_df$city, contrasts = FALSE)
                       ))
survivalDataset = cbind(survivalDataset,new_df)
survivalDataset$city = NULL
survivalDataset$countries = NULL
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

survivalDataset <- read.csv(file = "Data/LifeExpectancyData.csv")

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


