survivalDataset <- read.csv(file = "Data/covid_cv_folds/discharge_wo_argentina.csv")

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