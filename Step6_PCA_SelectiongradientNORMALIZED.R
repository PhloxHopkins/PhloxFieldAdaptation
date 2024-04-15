# Leaf traits PCA and selection gradient analysis
# December 6, 2023 edited Dec. 2023 
# Goal to summarize leaf trait variation in PCA and reconstitute selection gradients

# load packages
library(tidyverse)
library(ggplot2)
library(reshape2)

wd<-"/Users/robinhopkins/Robin/BenGoulet/ch2_field/field_manuscript/revisions/NatureComm/Resubmit/PCA"
setwd(wd)

# load data
# File is named '.csv'
traits <- read.csv("pop_leaf_traits.csv")

# Perform PCA on leaf trait data from wild-collected plants grown in Weld Hill greenhouse
# to 'traits' data, add lw_ratio and sla
traits$lw_ratio <- traits$length/traits$width
traits$sla <- traits$area/traits$dry.mass

attach(traits)
X <- cbind(length, width, area, lw_ratio, chloro, sla)

# check correlations
cor(X)

# Principal component analysis
pca1 <- princomp(X, scores=TRUE, cor=TRUE)
summary(pca1)
# Proportion of Variance explained
# Comp.1 = 0.449726
# Comp.2 = 0.2780881
# Comp.3 = 0.1822384

# Loadings of first three principal components, making Matrix E
traitloadings<-pca1$loadings[,1:3]
print(traitloadings)

#read in dataset with values for each pc and fitness
fit_tra_dis_pca<-read.csv("Dataset_fitness_leaftraits_distance_lat_long.csv")
# normalize fitness measures
# https://www.statology.org/how-to-normalize-data-in-r/
# for all species/all gardens
normfit_tra_dis_pca <- as_tibble(fit_tra_dis_pca)
normfit_tra_dis_pca$flowers_2019_total <- scale(fit_tra_dis_pca$flowers_2019)
normfit_tra_dis_pca$fruits_2019_total <- scale(fit_tra_dis_pca$fruits_2019)
normfit_tra_dis_pca$biomass_minus_bag_plus_fruits_2019 <- scale(fit_tra_dis_pca$biomass_2019)

normfit_tra_dis_pca$biomass_2019_normal<-(fit_tra_dis_pca$biomass_2019/(mean(fit_tra_dis_pca$biomass_2019)))
normfit_tra_dis_pca$flowers_2019_normal<-(fit_tra_dis_pca$flowers_2019/(mean(fit_tra_dis_pca$flowers_2019)))
normfit_tra_dis_pca$fruits_2019_normal<-(fit_tra_dis_pca$fruits_2019/(mean(fit_tra_dis_pca$fruits_2019)))



#subset data for just pilosa garden
LBLamodeapil <- as_tibble(normfit_tra_dis_pca) %>%
  filter(garden=="3-Pilosa-LBL")
#subset data for just pilosa in the pilosa garden
LBLpil <- as_tibble(normfit_tra_dis_pca) %>%
  filter(garden=="3-Pilosa-LBL") %>%
  filter(taxon=="pilosa")
 

#Reconstitute selection gradients on traits using method from Chong et al 2018 Evol. Letters
######linear regression of fitness and pc scores to create A matrix with regression coefficient for each PC
pcval<-c("leaf_traits_pc1","leaf_traits_pc2", "leaf_traits_pc3")
fitnesstraits<-c("flowers_2019_normal", "fruits_2019_normal", "biomass_2019_normal")

#function to collect the regression coefficients into a vector A
coefA<-function(data,fit,pc){
    coA<-c()
    for(j in 1:length(pc)){
      
      f<-as.formula(paste(fit,paste(pcval[j]), sep="~"))
      m<-lm(formula=f,data=data)
      coef<-summary(m)$coefficients[2,1]
      coA[j]=coef
      
      }
return(coA)
}


#making beta vector of selection gradients
make_b_matrix<-function(fittrait){
#Generate A matrix for all species with #flowers as fitness trait
allcoA<-coefA(LBLamodeapil,fittrait,pcval)
#generate A matrix for pilosa species with #flowers as fitness trait
pilcoA<-coefA(LBLpil,fittrait,pcval)
#calculate beta vectors for all species and just pilosa for flower num. and make into matrix
b_all<-traitloadings%*%allcoA
b_pil<-traitloadings%*%pilcoA
b_matrix<-matrix(,nrow=6,ncol=2)
b_matrix[,1]=b_all
b_matrix[,2]=b_pil
colnames(b_matrix)<-c("bioallsp", "Pilosa")
rownames(b_matrix)<-c("length","width","area","lw_ratio","chloro","sla")
return(b_matrix)
#return(pilcoA)
}

flower<-make_b_matrix(fitnesstraits[1])
flower
fruit<-make_b_matrix(fitnesstraits[2])
fruit
biomass<-make_b_matrix(fitnesstraits[3])
biomass

make_heat_fit<-function(bmatrix){
matmelt<-melt(bmatrix)
colnames(matmelt)<-c("y","x","value")
heatplot<-ggplot(matmelt,aes(x=x,y=y,fill=value)) +
  geom_tile() + 
  scale_fill_gradient2(low="purple",
                       mid="white",
                       high="orange") +
  geom_text(aes(label=round(value,digits=2)),color="black", size=12, fontface="bold")
return(heatplot)
}

quartz(width = 5, height = 12)
make_heat_fit(biomass)

