##---------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------
## 
## Created by: Ben Goulet-Scott; edited by Meghan Blumstein
## 
## created 9/13/2021 & updated 3/31/2023
##
## Description:
## calculate slope and R^2 for population random effects vs. distances
## Formerly 2021_9_13_population_effect_correlations.R
##
## Inputs: 
## From Step_1_Population_Effects.R; Step_1_Population_Random_Effects.csv
##
## Outputs: 
## Figures
##
## # references/links:
## https://www.learnbymarketing.com/tutorials/linear-regression-in-r/
## http://r-statistics.co/Linear-Regression.html
##
##---------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------

## Clears workspace
rm(list = ls())

##-----------------------------------------------------------------------------------------------
## USER INPUT 
##-----------------------------------------------------------------------------------------------

## Packages 
library(lme4)
require(plyr)

## Set Working Directory
wd <- #
setwd(wd)

## Load Data
pops <- read.csv("Step_1_Population_Random_Effects.csv")
#Runs models for results show in Figure 4 and Table S4
## Set Default R colors for easy plotting
r.colors <- palette()
r.colors[1:4] <- c("#F8766D", "#00BF7D", "#00B0F6", "#E76BF3")
palette(r.colors)

## Traits and Environmental Variables of Ultimate Interest
id_cols    <- c("label", "garden", "bed", "taxon", "population", "longitude", "latitude")
trait_cols <- c("herbivory_2019_ranef", "flowers_2019_ranef", "fruits_2019_ranef", "biomass_2019_ranef", "survival_2020_ranef")
env_cols   <- c("geographic_distance", "genetic_distance", "environmental_distance") 

##-----------------------------------------------------------------------------------------------
## CORRELATE POP RAN EFFECTS WITH ENV VARIABLES 
##-----------------------------------------------------------------------------------------------

## Function
correlate_vars <- function(data, variable_list, environment_list){
  
  ## For each Species
  species <- c("amoena", "pilosa", "deamii")
  
  out <- c()
  for(sp in species){
    
    ## Subset the Data
    sub <- data[data$taxon == sp,]
    
    for(i in variable_list){
      
      
      for(j in 1:length(environment_list)){
      
        ## Correlate
        model <- lm(sub[,i] ~ sub[,environment_list[j]], data = sub)
        
        ## Pull P values and slopes
        pval <- summary(model)$coefficients[2,4]
        m    <- coef(model)[2]
        
        ## Save information
        df <- data.frame(taxon = sp, env = environment_list[j], trait = i,pval = pval, m = m)
        
        ## Combine species
        out <- rbind(out, df)
      
      }
    }

  }

  
  ## Return formatted dataframe
  row.names(out) <- NULL
  return(out)
  
}

## Run Function
df <- correlate_vars(data = pops, variable_list = trait_cols, environment_list = env_cols)















