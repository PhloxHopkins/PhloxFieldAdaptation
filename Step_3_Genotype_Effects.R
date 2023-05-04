##---------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------
## 
## Created by: Ben Goulet-Scott; edited by Meghan Blumstein
## 
## created 6/23/2021 & updated 3/31/2023
##
## Description:
##  goal: to re-run "standard" models testing for interspecific local adaptation 
##
## Inputs: 
## 
##
## Outputs: 
## 
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
library(car)
library(multcomp)
library(ggplot2)

## Set Working Directory
wd <- "/Users/robinhopkins/Robin/BenGoulet/ch2_field/field_manuscript/Code_Data/meghan_code"
setwd(wd)

## Load Data
traits <- read.csv("Dataset_fitness_leaftraits_distance_lat_long.csv")
traits$taxonxgarden <- paste(traits$taxon, traits$garden, sep="_")
traits$taxonxgarden <- as.factor(traits$taxonxgarden)
## Set Default R colors for easy plotting
r.colors <- palette()
r.colors[1:4] <- c("#F8766D", "#00BF7D", "#00B0F6", "#E76BF3")
palette(r.colors)

## Traits and Environmental Variables of Ultimate Interest
id_cols    <- c("label", "garden", "bed", "taxon", "population", "longitude", "latitude","taxonxgraden")
trait_cols <- c("herbivory_2019", "flowers_2019", "fruits_2019", "biomass_2019", "survival_2020")
env_cols   <- c("geographic_distance", "genetic_distance", "environmental_distance") 

##-----------------------------------------------------------------------------------------------
## PREP DATA
##-----------------------------------------------------------------------------------------------

## Make biomass a log
traits$biomass_2019[which(traits$biomass_2019 < 0)] <- 0

traits$biomass_2019                                 <- log(traits$biomass_2019 + 1) ## Ben did 0.1, but might want to do 1


  ##-----------------------------------------------------------------------------------------------
  ##  EXTRACT FIXED EFFECTS AND CONTRASTS FROM POPULATION/GENOTYPE RANDOM EFFECTS MODELS
  ##-----------------------------------------------------------------------------------------------
  
  ## Function to run glm and extract fixed effects for each fitness trait with pop/geno random effects
  extract_fixedEffects <- function(data = traits, variable = "herbivory_2019", family = "bionmial"){
    
    ## Run Model
    model <- glmer(data[,variable] ~ taxon*garden + (1|population/genotype), data=data, family= family)
    rl <- Anova(model, type=3)
    return(rl)
    
  }    
  
  ##Function to extract contrasts for assessment of adaptive divergence for each fitness trait
  extract_Contrast<-function(data=traits, variable = "herbivory_2019", family="bionomial"){
    model <- glmer(data[,variable] ~ taxonxgarden + (1|population/genotype), data=data, family= family)
    l2 <- glht(model, linfct = mcp(taxonxgarden = "Tukey"))
    t1 <- summary(l2, test = adjusted(type = "none"))
    
    # put summary of contrasts into dataframe
    contrastDF <- data.frame(estimate = t1$test$coefficients, se = t1$test$sigma, TStat = t1$test$tstat, pval = t1$test$pvalues)
    
    # remove whitespace
    row.names(contrastDF) <- gsub(x = row.names(contrastDF), pattern = " ", replacement = "")
    gps <- t(as.data.frame(strsplit(x = row.names(contrastDF), split = '[x]|[_]|[-]')))
    row.names(gps) <- NULL
    
    cdf <- cbind(gps, contrastDF)
    cdf<-cdf[-c(2,4,6,8)]
    colnames(cdf)[1:4] <- c("species1", "garden1",  'species2', 'garden2')
    cdf$garden1<-sub("A","a",cdf$garden1)
    cdf$garden2<-sub("A","a",cdf$garden2)
    cdf$garden1<-sub("D","d",cdf$garden1)
    cdf$garden2<-sub("D","d",cdf$garden2)
    cdf$garden1<-sub("P","p",cdf$garden1)
    cdf$garden2<-sub("P","p",cdf$garden2)
    
    # subset contrast DF, so we get only the ones that interest us for adaptive divergence
    cdfslim <- cdf[as.character(cdf$species1) == as.character(cdf$species2) & as.character(cdf$species1) == as.character(cdf$garden1)|
                     as.character(cdf$species1) == as.character(cdf$species2) & as.character(cdf$species1) == as.character(cdf$garden2)|
                     as.character(cdf$garden1) == as.character(cdf$garden2) & as.character(cdf$species1) == as.character(cdf$garden1)|
                     as.character(cdf$garden1) == as.character(cdf$garden2) & as.character(cdf$species2) == as.character(cdf$garden1), ]
    # round p-values for table
    cdfslim$pval <- round(cdfslim$pval, digits = 5)
    contrasts<-cdfslim
    return(contrasts)
  }

  ## Run for each variable
  ## Returns a table with fixed effects 
  ## and a table of contrasts we care about
  ## 
#Herbivory 2019
FE_herbivory <- extract_fixedEffects(data = traits, variable = trait_cols[1], family = "binomial")
write.csv(FE_herbivory, file = "FE_Herbivory.csv")

con_herbivory<-extract_Contrast(data=traits, variable=trait_cols[1], family="binomial")
write.csv(con_herbivory, file = "contrast_Herbivory.csv")

## Flowers 2019
FE_flower <- extract_fixedEffects(data = traits, variable = trait_cols[2], family = "poisson")
write.csv(FE_flower, file="FE_Flower.csv")

con_flower<-extract_Contrast(data=traits, variable=trait_cols[2], family="poisson")
write.csv(con_flower, file = "contrast_flower.csv")

## Fruits 2019

FE_fruit <- extract_fixedEffects(data = traits, variable = trait_cols[3], family = "poisson")
write.csv(FE_fruit, file="FE_Fruit.csv")

con_fruit<-extract_Contrast(data=traits, variable=trait_cols[3], family="poisson")
write.csv(con_fruit, file = "contrast_fruit.csv")

## Biomass
## Start by setting all 0 or negatives to 0

FE_biomas <- extract_fixedEffects(data = traits, variable = trait_cols[4], family = "gaussian")
write.csv(FE_biomas, file="FE_biomass.csv")

con_biomass<-extract_Contrast(data=traits, variable=trait_cols[4], family="gaussian")
write.csv(con_biomass, file = "contrast_biomass.csv")

## Survival 2020
FE_survival <- extract_fixedEffects(data = traits, variable = trait_cols[5], family = "binomial")
write.csv(FE_survival, file="FE_survival.csv")

con_survival<-extract_Contrast(data=traits, variable=trait_cols[5], family="binomial")
write.csv(con_survival, file = "contrast_survival.csv")

## Graph estimates from contrasts
contSum <- read.csv("contrast_combined.csv")
ggplot(contSum, aes(x=LA.estimate, y=focal.species, color=pvalsig)) +
  geom_jitter(width=0, height=0.1, alpha=0.7, size=3.5) +
  scale_color_manual(values=c("#999999", "#000000"))+
  geom_vline(xintercept=0, linetype=2) +
  ggtitle("Estimates of local adaptation") +
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_blank())
  
contLA<-contSum[which(contSum$LvH=="L"),]


ggplot(contLA, aes(x=LA.estimate, y=focal.species, color=pvalsig)) +
  geom_jitter(width=0, height=0.1, alpha=0.7, size=3.5) +
  scale_color_manual(values=c("#999999", "#000000"))+
  geom_vline(xintercept=0, linetype=2) +
  ggtitle("Estimates of local adaptation") +
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_blank())
