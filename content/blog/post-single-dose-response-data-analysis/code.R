#Loads libraries
library(drc)
library(ggplot2)
library(magrittr)
library(synergyfinder)
library(Cairo)



#### Synthetic curve ####
data_synthetic = data.frame(c(0.1,1,2,3,10, 100), c(0,20,40,60, 80, 100))
names(data_synthetic) <- c('concentration', 'inhibition')
plot(data_synthetic, log='x', xlab = 'Concentration (microM)', ylab= 'Inhibition')

Cairo::Cairo(
  16, #length
  11, #width
  file = paste("images/image_scatter_dose_response", ".png", sep = ""),
  type = "png", #tiff
  bg = "white", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)
plot(data_synthetic, log='x', xlab = 'Drug concentration (microM)', ylab= 'Inhibition') #p is your graph object 
dev.off()



#Loads toy data matrix
data = data.frame(read.csv('toy_drug_combination_data.csv'))
head(data)


#### Monotherapy analyses ####
##Drug row
data_row = data[data$conc_c == 0,]


#A general model fitting function for analysis of concentration/dose/time-effect/response data.
fitted_curve <- drm(formula = inhibition ~ conc_r,
             data = data_row,
             fct = LL.4()
             )
summary(fitted_curve)
plot(fitted_curve, xlab = 'Concentration (microM)', ylab= 'Inhibition')



#E0
#Effect in the absence of drug
#E0 = c
E0 <- fitted_curve$coefficients[2]

#Einf
#Maximum effect of a drug
#Einf = d
Einf <- fitted_curve$coefficients[3]

#EC50
#Half-maximal effective concentration is the concentration of a drug/antibody at which 50% of 
#its maximal response is induced. EC50 is normally measured in molar concentrations and is
#used as a measure of agonist drug potency - the lower the EC50 value, 
#the lower the concentration of drug required to elicit a 50% maximal response and the greater potency of the drug. 
#EC50 = e
EC50 <- fitted_curve$coefficients[4]


#Slope of DR curve
#A steeper curve signifies increased response rates upon increasing dose, indicating the drug has a greater 
#potency. A more gradual rising curve indicates that an increasing drug dose is required to elicit a response, 
#suggesting decreased potency. 
#Slope = hill coefficient = -b
H <-  -fitted_curve$coefficients[1]


#IC50
# Half-maximal inhibitory concentration is the concentration of an inhibitor (e.g. drug) required for 50% inhibition 
#of a biological target or biochemical function. 
#IC50 is normally measured in molar concentrations and is used as a measure of antagonist drug potency - 
#the lower the IC50 value, the more potent the drug/inhibitor is. 
#Careful here -- if the drug doesn't achieve 50% of inhibition, this value is extrapolated!

IC50 <- exp(log((((Einf - E0) / (50 - E0) ) - 1) * EC50**(-H))/ (-H))




##Drug col
data_col = data[data$conc_r == 0,]

#A general model fitting function for analysis of concentration/dose/time-effect/response data.
fitted_curve <- drm(formula = inhibition ~ conc_c,
                    data = data_col,
                    fct = LL.4()
)
summary(fitted_curve)
plot(fitted_curve, xlab = 'Concentration (microM)', ylab= 'Inhibition')


#E0
#Effect in the absence of drug
#E0 = c
E0 <- fitted_curve$coefficients[2]

#Einf
#Maximum effect of a drug
#Einf = d
Einf <- fitted_curve$coefficients[3]

#EC50
#Half-maximal effective concentration is the concentration of a drug/antibody at which 50% of 
#its maximal response is induced. EC50 is normally measured in molar concentrations and is
#used as a measure of agonist drug potency - the lower the EC50 value, 
#the lower the concentration of drug required to elicit a 50% maximal response and the greater potency of the drug. 
#EC50 = e
EC50 <- fitted_curve$coefficients[4]


#Slope of DR curve
#A steeper curve signifies increased response rates upon increasing dose, indicating the drug has a greater 
#potency. A more gradual rising curve indicates that an increasing drug dose is required to elicit a response, 
#suggesting decreased potency. 
#Slope = hill coefficient = -b
H <-  -fitted_curve$coefficients[1]


#IC50
# Half-maximal inhibitory concentration is the concentration of an inhibitor (e.g. drug) required for 50% inhibition 
#of a biological target or biochemical function. 
#IC50 is normally measured in molar concentrations and is used as a measure of antagonist drug potency - 
#the lower the IC50 value, the more potent the drug/inhibitor is. 
#Careful here -- if the drug doesn't achieve 50% of inhibition, this value is extrapolated!

IC50 <- exp(log((((Einf - E0) / (50 - E0) ) - 1) * EC50**(-H))/ -H)





#### Synergy analysis ####
data = data.frame(read.csv('toy_drug_combination_data.csv'))


names(data)[names(data) == "inhibition"] <- "response"

data <- ReshapeData(
  data,
  impute = TRUE,
  noise = TRUE, correction = "non", data.type = "inhibition"
)


ZIP <- CalculateSynergy(data, method = "ZIP", adjusted = TRUE)
synergy_score <- mean(ZIP$scores$'10'[ZIP$scores$'10' != 0])






