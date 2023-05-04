##---------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------
## 
## Created by: Ben Goulet-Scott; edited by Meghan Blumstein
## 
## created 9/2/2021 & updated 1/11/2023
##
## Description:
## Get population level random effects of traits (for examining local adaptation) 
##
## Inputs: 
## All measured garden trait data  and garden distance data using 
## Dataset_fitness_leaftraits_distance_lat_long.csv
##
##
## Outputs: 
##
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
traits <- read.csv("Dataset_fitness_leaftraits_distance_lat_long.csv")
##-----------------------------------------------------------------------------------------------
## FORMAT DATA
##-----------------------------------------------------------------------------------------------

## Merge distance with trait data
traits <- merge(traits, dist)

## Traits and Environmental Variables of Ultimate Interest
id_cols    <- c("label", "garden", "bed", "taxon", "population", "individual", "genotype", "longitude", "latitude")
trait_cols <- c("herbivory_2019", "flowers_2019", "fruits_2019", "biomass_2019", "survival_2020")
env_cols   <- c("geographic_distance", "genetic_distance", "environmental_distance") 

## Subset traits to just those of interest
traits     <- traits[,c(id_cols, trait_cols, env_cols)]

##-----------------------------------------------------------------------------------------------
## POPULATION RANDOM EFFECTS MODELS
##-----------------------------------------------------------------------------------------------

## Function
extract_randomEffects <- function(data = traits, variable = "herbivory_2019", family = "bionmial"){
  
  ## For each Species
  species <- c("amoena", "pilosa", "deamii")
  
  out <- c()
  for(sp in species){
    
  
    if(sp == "deamii"){
      
      ## Subset the Data
      sub <- data[data$taxon == sp,]
      
      ## Drop unused factors in subset
      sub <- droplevels(sub)
        
      ## deamii has replicates of genotypes within gardens
      model <- glmer(sub[,variable] ~ (1|garden/population/genotype), data = sub, family = family)
      
    }else{
      
      ## Subset the Data
      sub <- data[data$taxon == sp,]
      
      ## Drop unused factors in subset
      sub <- droplevels(sub)
      
      ## Run Model for each of the three species
      model <- glmer(sub[,variable] ~ (1|garden/population), data = sub, family = family)

    }
    
    
    ## Pull out Random Effects
    df <- as.data.frame(ranef(model)$'population:garden' + fixef(model))
    
    ## Format Dataframe
    df <- data.frame(pop_gar = rownames(df), taxon = sub$taxon[1], ranef = df$`(Intercept)`)
    
    ## save all the species data
    out <- rbind(out, df)
    
  }
  ## Rename the random effects with the variable name
  colnames(out)[3] <- paste0(variable, "_ranef")
  
  ## Return formatted random effects dataframe
  return(out)
  
}


## Run for each variable

  ## Herbivory 2019
  hist(traits[,trait_cols[1]], xlab = trait_cols[1])
  re1 <- extract_randomEffects(data = traits, variable = trait_cols[1], family = "binomial")

  ## Flowers 2019
  hist(traits[,trait_cols[2]], xlab = trait_cols[2])
  re2 <- extract_randomEffects(data = traits, variable = trait_cols[2], family = "poisson")
  
  ## Fruits 2019
  hist(traits[,trait_cols[3]], xlab = trait_cols[3])
  re3 <- extract_randomEffects(data = traits, variable = trait_cols[3], family = "poisson")
  
  ## Biomass
  ## Start by setting all 0 or negatives to 0
  ## And add 0.1 to all so you can log
  traits$biomass_2019[which(traits$biomass_2019 < 0)] <- 0
  traits$biomass_2019                                 <- log(traits$biomass_2019 + 1) ## Ben did 0.1, but might want to do 1
  hist(log( traits[,trait_cols[4]] ), xlab = trait_cols[4])
  re4 <- extract_randomEffects(data = traits, variable = trait_cols[4], family = "gaussian")

  ## Survival 2020
  hist(traits[,trait_cols[5]], xlab = trait_cols[5])
  re5 <- extract_randomEffects(data = traits, variable = trait_cols[5], family = "binomial")


##-----------------------------------------------------------------------------------------------
## FORMAT FINAL DATA FRAME
##-----------------------------------------------------------------------------------------------

## Get just info for garden, taxon, & populations
info <- ddply(traits, .(garden, taxon, population), summarize, 
              latitude = mean(latitude), 
              longitude = mean(longitude), 
              genetic_distance = mean(genetic_distance),
              geographic_distance = mean(geographic_distance), 
              environmental_distance = mean(environmental_distance)
              )
info$pop_gar <- paste0(info$population, ":", info$garden)
  
## Merge all the dataframes
output <- merge(info, re1)
output <- merge(output, re2)
output <- merge(output, re3)
output <- merge(output, re4)
output <- merge(output, re5)

## Write CSV
write.csv(output, "Step_1_Population_Random_Effects.csv", row.names = F)










