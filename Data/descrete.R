survivalDataset <- read.csv(file = "Data/All_Data_updated_may2011_CLEANED.csv")
names(survivalDataset)[c(1,2)] = c("time", "delta")
survivalDataset$delta = 1 - survivalDataset$delta 


survivalDataset$ALBUMIN = NULL
survivalDataset$CREATININE_SERUM = NULL
survivalDataset$HGB = NULL
survivalDataset$CALCIUM_SERUM = NULL
survivalDataset$WBC_COUNT = NULL
survivalDataset$PLATELET = NULL
survivalDataset$LYMPHOCYTES = NULL
survivalDataset$LDH_SERUM = NULL
survivalDataset$GRANULOCYTES = NULL
survivalDataset$AGE = NULL
survivalDataset$BMI = NULL

survivalDataset$STAGE_1 = NULL
survivalDataset$STAGE_2 = NULL
survivalDataset$STAGE_3 = NULL
survivalDataset$STAGE_4 = NULL
survivalDataset$PERFORMANCE_STATUS_0 = NULL
survivalDataset$PERFORMANCE_STATUS_1 = NULL
survivalDataset$PERFORMANCE_STATUS_2 = NULL
survivalDataset$PERFORMANCE_STATUS_3 = NULL
survivalDataset$PERFORMANCE_STATUS_4 = NULL

site = integer(nrow(survivalDataset))
for(i in 1:nrow(survivalDataset)){
  if(survivalDataset$SITE_BRUNCHUS_LUNG[i] == 1){
    site[i] = 1
  }else if(survivalDataset$SITE_COLORECTAL[i] == 1){
    site[i] = 2
  }else if(survivalDataset$SITE_HEAD_AND_NECK[i] == 1){
    site[i] = 3
  }else if(survivalDataset$SITE_ESOPHAGUS[i] == 1){
    site[i] = 4
  }else if(survivalDataset$SITE_PANCREAS[i] == 1){
    site[i] = 5
  }else if(survivalDataset$SITE_STOMACH[i] == 1){
    site[i] = 6
  }else if(survivalDataset$SITE_OTHER_DIGESTIVE[i] == 1){
    site[i] = 7
  }
}

survivalDataset$SITE = site
survivalDataset$SITE_BRUNCHUS_LUNG = NULL
survivalDataset$SITE_COLORECTAL = NULL
survivalDataset$SITE_HEAD_AND_NECK = NULL
survivalDataset$SITE_ESOPHAGUS = NULL
survivalDataset$SITE_PANCREAS = NULL
survivalDataset$SITE_STOMACH = NULL
survivalDataset$SITE_OTHER_DIGESTIVE = NULL

survivalDataset$BOX1_SCORE = NULL
survivalDataset$BOX2_SCORE = NULL
survivalDataset$BOX3_SCORE = NULL

survivalDataset$NO_APPETITE = NULL
survivalDataset$TASTE_FUNNY = NULL
survivalDataset$FEEL_FULL = NULL
survivalDataset$SSWALLOW = NULL
survivalDataset$SWALLOW = NULL
survivalDataset$DRY_MOUTH = NULL
survivalDataset$SMELL = NULL
survivalDataset$SORE_MOUTH = NULL
survivalDataset$PAIN = NULL
survivalDataset$OTHER = NULL
survivalDataset$DENTAL_PROBLEM = NULL
survivalDataset$VOMIT = NULL
survivalDataset$NAUSEA = NULL
survivalDataset$DIARRHEA = NULL
survivalDataset$CONSTIPATION = NULL
survivalDataset$AGE65 = NULL
survivalDataset$MISC = NULL
survivalDataset$WEIGHT_CHANGEPOINT = NULL

featureIdx = sample(3:ncol(survivalDataset),10)

survivalDataset = survivalDataset[,c(1,2,featureIdx)]

