## Comparison of leaf trait data collected from Field and Greenhouse plants
## All data from population 729 including both Amoena and Pilosa species
## Results show in Table S9 and Figure S2


## Packages 
library(lme4)
require(plyr)
library(car)
library(multcomp)
library(ggplot2)

## Set Working Directory
wd <- #
setwd(wd)
leaf <- read.csv("combined_leaftraits.csv")
leaf$taxonXsource<-paste(leaf$taxon, leaf$FvG, sep="_")
leaf$taxonXsource<-as.factor(leaf$taxonXsource)

lw<-lm(formula=length_width_ratio~taxon+ FvG +taxon*FvG, data=leaf)
Anova(lw, type=3)

w<-lm(formula=width~taxon+ FvG +taxon*FvG, data=leaf)
Anova(w, type=3)

area<-lm(formula=area~taxon+ FvG +taxon*FvG, data=leaf)
Anova(area, type=3)

SPA<-lm(formula=specific_leaf_area~taxon+ FvG +taxon*FvG, data=leaf)
Anova(SPA, type=3)

l<-lm(formula=length~taxon+ FvG +taxon*FvG, data=leaf)
Anova(l, type=3)

w2<-lm(formula=width~taxonXsource, data=leaf)
Anova(w2, type=3)
w2con <- glht(w2, linfct = mcp(taxonXsource = "Tukey"))
t1 <- summary(w2con, test = adjusted(type = "none"))

SPA2<-lm(formula=specific_leaf_area~taxonXsource, data=leaf)
Anova(SPA2, type=3)
SPA2con <- glht(SPA2, linfct = mcp(taxonXsource = "Tukey"))
t1 <- summary(SPA2con, test = adjusted(type = "none"))


greenhouse <- leaf[leaf$FvG == "G",]
field  <- leaf[leaf$FvG == "F",]

########### Figure S2
trts <- c("length", "width", "area", "length_width_ratio", "specific_leaf_area")
quartz(width = 10, height = 15)
par(mfrow = c(3,2), mar = c(3.5, 3, 0.1, 0.1), mgp = c(1.5, 0.5, 0))
for(i in 1:length(trts)){
  print(trts[i])
  mycols<-c("#991107","#F8766D","#003388","#619CFF")
  boxplot(leaf[,trts[i]] ~ interaction(leaf$FvG, leaf$taxon),
          names = rep(c("Field", "Greenhouse"),2),
          col = mycols, 
          #xlab = "Field vs. Greenhouse",
          xlab=NULL,
          ylab = trts[i])
  
}

