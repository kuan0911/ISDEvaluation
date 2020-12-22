survivalDataset <- read.csv(file = "Data/covid_latestdata.csv", na.strings='')
survivalDataset <- read.csv(file = "Data/covid_clean_data_v3.csv", na.strings='')
survivalDataset <- read.csv(file = "Data/covid_hospitalized_data.csv", na.strings='')

survivalDataset$additional_information = NULL

survivalDataset = survivalDataset[!is.na(survivalDataset$sex),]
survivalDataset = survivalDataset[!is.na(survivalDataset$age),]
survivalDataset = survivalDataset[!is.na(survivalDataset$latitude),]
survivalDataset = survivalDataset[!is.na(survivalDataset$longitude),]
survivalDataset = survivalDataset[!is.na(survivalDataset$date_confirmation),]

survivalDataset = survivalDataset[!is.na(survivalDataset$date_admission_hospital),]

survivalDataset = survivalDataset[is.na(survivalDataset$outcome)|survivalDataset$outcome!='Death',]
survivalDataset = survivalDataset[is.na(survivalDataset$outcome)|survivalDataset$outcome!='death',]
survivalDataset = survivalDataset[is.na(survivalDataset$outcome)|survivalDataset$outcome!='died',]
survivalDataset = survivalDataset[is.na(survivalDataset$outcome)|survivalDataset$outcome!='Died',]
survivalDataset = survivalDataset[is.na(survivalDataset$outcome)|survivalDataset$outcome!='Dead',]
survivalDataset = survivalDataset[is.na(survivalDataset$outcome)|survivalDataset$outcome!='dead',]

survivalDataset = survivalDataset[is.na(survivalDataset$outcome)|survivalDataset$outcome=='discharge'|survivalDataset$outcome=='Discharged'|survivalDataset$outcome=='discharged',]


#survivalDataset = survivalDataset[,c('age','sex','event_time','censor')]

survivalDataset = survivalDataset[!grepl('-',survivalDataset$date_admission_hospital),]
survivalDataset = survivalDataset[!grepl('-',survivalDataset$date_onset_symptoms),]
survivalDataset = survivalDataset[!grepl('-',survivalDataset$date_confirmation),]
survivalDataset = survivalDataset[!grepl('-',survivalDataset$date_death_or_discharge),]
survivalDataset = survivalDataset[!grepl('-',survivalDataset$age),]




survivalDataset$date_admission_hospital = as.Date(survivalDataset$date_admission_hospital,tryFormats = c("%d.%m.%y"))
survivalDataset$date_onset_symptoms = as.Date(survivalDataset$date_onset_symptoms,tryFormats = c("%d.%m.%y"))
survivalDataset$date_confirmation = as.Date(survivalDataset$date_confirmation,tryFormats = c("%d.%m.%y"))
survivalDataset$date_death_or_discharge = as.Date(survivalDataset$date_death_or_discharge,tryFormats = c("%d.%m.%y"))

for(k in 1:nrow(survivalDataset)) {
  if(is.na(survivalDataset[k,'date_onset_symptoms'])&!is.na(survivalDataset[k,'date_admission_hospital'])&!is.na(survivalDataset[k,'date_confirmation'])) {
    survivalDataset[k,'date_onset_symptoms'] = survivalDataset[k,'date_admission_hospital']
  }
  if(is.na(survivalDataset[k,'date_onset_symptoms'])&is.na(survivalDataset[k,'date_admission_hospital'])&!is.na(survivalDataset[k,'date_confirmation'])) {
    survivalDataset[k,'date_onset_symptoms'] = survivalDataset[k,'date_confirmation']
  }
}
for(k in 1:nrow(survivalDataset)) {
  if(!is.na(survivalDataset[k,'outcome'])) {
    outcome = survivalDataset[k,'outcome']
    if(outcome=='Discharge'|outcome=='discharged'|outcome=='Discharged') {
      survivalDataset[k,'outcome'] = 'discharge'
    }
    if(outcome=='hospitalize'|outcome=='Hospitalize'|outcome=='hospitalized'|outcome=='Hospitalized') {
      survivalDataset[k,'outcome'] = 'discharge'
    }
    if(outcome=='recovered'|outcome=='recover'|outcome=='release'|outcome=='Released'|outcome=='Recovered'|outcome=='Recover'|outcome=='Release'|outcome=='Released') {
      survivalDataset[k,'outcome'] = 'discharge'
    }
    if(outcome=='Death'|outcome=='died'||outcome=='Died'|outcome=='dead'|outcome=='Dead') {
      survivalDataset[k,'outcome'] = 'death'
    }
  }
}
survivalDataset$event_time = rep(0,nrow(survivalDataset))
for(k in 1:nrow(survivalDataset)) {
  if(!is.na(survivalDataset[k,'date_death_or_discharge'])) {
    survivalDataset[k,'event_time'] = survivalDataset[k,'date_death_or_discharge'] - survivalDataset[k,'date_admission_hospital']
  }
}
maxTime = max(survivalDataset$event_time)
for(k in 1:nrow(survivalDataset)) {
  if(is.na(survivalDataset[k,'date_death_or_discharge'])) {
    survivalDataset[k,'event_time'] = maxTime
  }
}

survivalDataset = survivalDataset[survivalDataset$event_time>0,]

sum(!is.na(survivalDataset$outcome)&!is.na(survivalDataset$date_death_or_discharge))
sum(!is.na(survivalDataset$outcome))
sum(survivalDataset$outcome=='discharge',na.rm=TRUE)
sum(survivalDataset$outcome=='death',na.rm=TRUE)

survivalDataset = survivalDataset[!is.na(survivalDataset$outcome)&!is.na(survivalDataset$date_death_or_discharge),]
survivalDataset = survivalDataset[survivalDataset$country!='China',]
survivalDataset = survivalDataset[!is.na(survivalDataset$date_admission_hospital),]
survivalDataset = survivalDataset[is.na(survivalDataset$outcome),]



