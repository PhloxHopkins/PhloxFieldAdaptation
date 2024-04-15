# Leaf trait - fitness trait analyses (traits vs. fitness)
# Apr. 7, 2022
# goal: ask whether fitness (field) is associated with leaf traits (greenhouse) and whether relationships change based on within vs. between species variation

# load packages
library(tidyverse)
library(lme4)
require(plyr)
library(car)
library(multcomp)
library(gridExtra)
# load data
fit_tra_dis_pca <- read.csv("/Users/robinhopkins/Robin/BenGoulet/ch2_field/field_manuscript/Code_Data/Dataset_fitness_leaftraits_distance.csv")

# normalize fitness measures
# https://www.statology.org/how-to-normalize-data-in-r/
# for all species/all gardens
normfit_tra_dis_pca <- as_tibble(fit_tra_dis_pca)
normfit_tra_dis_pca$flowers_2019_total <- (fit_tra_dis_pca$flowers_2019/(mean(fit_tra_dis_pca$flowers_2019)))
normfit_tra_dis_pca$fruits_2019_total <- (fit_tra_dis_pca$fruits_2019/(mean(fit_tra_dis_pca$fruits_2019)))
normfit_tra_dis_pca$biomass_minus_bag_plus_fruits_2019 <- (fit_tra_dis_pca$biomass_2019/(mean(fit_tra_dis_pca$biomass_2019)))
# herbivory is binary, should not be scaled


#function to make relative fitness measures


# filter data
# 3-Pilosa-LBL
LBLamodeapil <- as_tibble(normfit_tra_dis_pca) %>%
  filter(garden=="3-Pilosa-LBL")
LBLpil <- as_tibble(normfit_tra_dis_pca) %>%
  filter(garden=="3-Pilosa-LBL") %>%
  filter(taxon=="pilosa")
LBLamo <- as_tibble(normfit_tra_dis_pca) %>%
  filter(garden=="3-Pilosa-LBL") %>%
  filter(taxon=="amoena")
LBLdea <- as_tibble(normfit_tra_dis_pca) %>%
  filter(garden=="3-Pilosa-LBL") %>%
  filter(taxon=="deamii")
############ models
#Determine % variation of leaf trait pcs that are explained by taxon identity. 
pc1<-lm(formula=leaf_traits_pc1~taxon, data=LBLamodeapil)
Anova(pc1)
#R^2 =0.72, Fvalue=417.34, p<0.001

pc2<-lm(formula=leaf_traits_pc2~taxon, data=LBLamodeapil)
Anova(pc2)
#R^2=0.005, Fvalue=0.8, p=0.451

################ models to determine: 
## 1) if fitness traits are explained by taxon, PC or interaction and (Results in Table S7)
## 2) simple linear models with all or just one species (Results in Table S8)


#test with traits: Flower length/width
fl1tax<-lm(formula=flowers_2019_total~area+ taxon +taxon*leaf_traits_pc1, data=LBLamodeapil)
Anova (fl1tax, type=3)
fl1all<-lm(formula=fruits_2019_total~area +length_width_ratio+chlorophyll+specific_leaf_area+ length+width, data=LBLamodeapil)
fl1p<-lm(formula=fruits_2019_total~area+length_width_ratio +chlorophyll+specific_leaf_area+ length+width, data=LBLpil)

summary(fl1all)
summary(fl1p)






#flower
#pc1
fl1tax<-lm(formula=flowers_2019_total~leaf_traits_pc1+ taxon +taxon*leaf_traits_pc1, data=LBLamodeapil)
Anova (fl1tax, type=3)
fl1all<-lm(formula=flowers_2019_total~leaf_traits_pc1, data=LBLamodeapil)
fl1p<-lm(formula=flowers_2019_total~leaf_traits_pc1, data=LBLpil)
fl1a<-lm(formula=flowers_2019_total~leaf_traits_pc1, data=LBLamo)
fl1d<-lm(formula=flowers_2019_total~leaf_traits_pc1, data=LBLdea)
summary(fl1all)
summary(fl1p)
summary(fl1a)
summary(fl1d)
#pc2
fl2tax<-lm(formula=flowers_2019_total~leaf_traits_pc2+ taxon +taxon*leaf_traits_pc2, data=LBLamodeapil)
Anova (fl2tax,type=3)
summary(fl2tax)

fl2all<-lm(formula=flowers_2019_total~leaf_traits_pc2, data=LBLamodeapil)
fl2p<-lm(formula=flowers_2019_total~leaf_traits_pc2, data=LBLpil)
fl2a<-lm(formula=flowers_2019_total~leaf_traits_pc2, data=LBLamo)
fl2d<-lm(formula=flowers_2019_total~leaf_traits_pc2, data=LBLdea)
summary(fl2p)
summary(fl2a)
summary(fl2d)
summary(fl2all)

#Fruit
#pc1
ft1tax<-lm(formula=fruits_2019_total~leaf_traits_pc1+ taxon +taxon*leaf_traits_pc1, data=LBLamodeapil)
Anova (ft1tax, type=3)
ft1all<-lm(formula=fruits_2019_total~leaf_traits_pc1, data=LBLamodeapil)
ft1p<-lm(formula=fruits_2019_total~leaf_traits_pc1, data=LBLpil)
ft1a<-lm(formula=fruits_2019_total~leaf_traits_pc1, data=LBLamo)
ft1d<-lm(formula=fruits_2019_total~leaf_traits_pc1, data=LBLdea)
summary(ft1p)
summary(ft1a)
summary(ft1d)
summary(ft1all)

#pc2
ft2tax<-lm(formula=fruits_2019_total~leaf_traits_pc2+ taxon +taxon*leaf_traits_pc2, data=LBLamodeapil)
Anova (ft2tax,type=3)
summary(ft2tax)

ft2all<-lm(formula=fruits_2019_total~leaf_traits_pc2, data=LBLamodeapil)
ft2p<-lm(formula=fruits_2019_total~leaf_traits_pc2, data=LBLpil)
ft2a<-lm(formula=fruits_2019_total~leaf_traits_pc2, data=LBLamo)
ft2d<-lm(formula=fruits_2019_total~leaf_traits_pc2, data=LBLdea)
summary(ft2p)
summary(ft2a)
summary(ft2d)
summary(ft2all)

###biomass
#pc1
Anova (bio1)
bio1tax<-lm(formula=biomass_minus_bag_plus_fruits_2019~leaf_traits_pc1+ taxon +taxon*leaf_traits_pc1, data=LBLamodeapil)
Anova (bio1tax, type=3)
bio1all<-lm(formula=biomass_minus_bag_plus_fruits_2019~leaf_traits_pc1, data=LBLamodeapil)
bio1p<-lm(formula=biomass_minus_bag_plus_fruits_2019~leaf_traits_pc1, data=LBLpil)
bio1a<-lm(formula=biomass_minus_bag_plus_fruits_2019~leaf_traits_pc1, data=LBLamo)
bio1d<-lm(formula=biomass_minus_bag_plus_fruits_2019~leaf_traits_pc1, data=LBLdea)
summary(bio1p)
summary(bio1a)
summary(bio1d)
summary(bio1all)
#pc2
bio2<-lm(formula=biomass_minus_bag_plus_fruits_2019~leaf_traits_pc2, data=LBLamodeapil)
Anova (bio2)
bio2tax<-lm(formula=biomass_minus_bag_plus_fruits_2019~leaf_traits_pc2+ taxon +taxon*leaf_traits_pc2, data=LBLamodeapil)
Anova (bio2tax, type=3)

bio2all<-lm(formula=biomass_minus_bag_plus_fruits_2019~leaf_traits_pc2, data=LBLamodeapil)
bio2p<-lm(formula=biomass_minus_bag_plus_fruits_2019~leaf_traits_pc2, data=LBLpil)
bio2a<-lm(formula=biomass_minus_bag_plus_fruits_2019~leaf_traits_pc2, data=LBLamo)
bio2d<-lm(formula=biomass_minus_bag_plus_fruits_2019~leaf_traits_pc2, data=LBLdea)
summary(bio2p)
summary(bio2a)
summary(bio2d)
summary(bio2all)

########################################### Plot

LBLamodeapil$biomass_2019g<-LBLamodeapil$biomass_2019/1000
LBLamodeapil$biomass_2019g<-round(LBLamodeapil$biomass_2019g, digits = 2)
Pilgar<-fitness[which(fitness$garden=="3-Pilosa-LBL"),]


graphtheme<-theme_bw()+theme(
  axis.title.x = element_text(size=22, face="bold"),
  axis.title.y = element_text(size=22, face="bold"),
  axis.text = element_text(size = 16),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.position="none")

graphtheme2<-theme_bw()+theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text = element_text(size = 16),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.position="none")
#fruits
pc1fruit<-ggplot(LBLamodeapil, aes(leaf_traits_pc1, fruits_2019,)) +
  xlab("Leaf Trait PC1") + ylab("Fruit #")+
  geom_point(size=1.5, alpha=0.25, aes(color=taxon))  +scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF"))+
  geom_smooth(method=lm, linewidth=1.25, se=FALSE, color="darkgray")+
  geom_smooth(method=lm, aes(color=taxon),linewidth=1, se=FALSE)


pc2fruit<-ggplot(LBLamodeapil, aes(leaf_traits_pc2, fruits_2019)) +
  xlab("Leaf Trait PC2") + ylab("Fruit #")+
  geom_point(size=1.5,alpha=0.25, aes(color=taxon))+scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF"))+
  geom_smooth(method=lm, linewidth=1.25, se=FALSE, color="darkgray")+
  geom_smooth(method=lm, aes(color=taxon),linewidth=1, se=FALSE)

#biomass
pc1biomass<-ggplot(LBLamodeapil, aes(leaf_traits_pc1, biomass_2019g,)) +
  xlab("Leaf Trait PC1") + ylab("Biomass")+
  geom_point(size=1.5,alpha=0.25, aes(color=taxon))+scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF"))+
  geom_smooth(method=lm, linewidth=1.25, se=FALSE, color="darkgray")+
  geom_smooth(method=lm, aes(color=taxon),linewidth=1, se=FALSE)

pc2biomass<-ggplot(LBLamodeapil, aes(leaf_traits_pc2, biomass_2019g)) +
  xlab("Leaf Trait PC2") + ylab("Biomass")+
  geom_point(size=1.5,alpha=0.25, aes(color=taxon))+scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF"))+
  geom_smooth(method=lm, linewidth=1.25, se=FALSE, color="darkgray")+
  geom_smooth(method=lm, aes(color=taxon),linewidth=1, se=FALSE)

#flowers
pc1flw<-ggplot(LBLamodeapil,alpha=0.25, aes(leaf_traits_pc1, flowers_2019,)) +
  xlab("Leaf Trait PC1") + ylab("Flower #")+
  geom_point(size=1.5,alpha=0.25, aes(color=taxon))+scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF"))+
  geom_smooth(method=lm, linewidth=1.25, se=FALSE, color="darkgray")+
  geom_smooth(method=lm, aes(color=taxon),linewidth=1, se=FALSE)

pc2flw<-ggplot(LBLamodeapil, aes(leaf_traits_pc2, flowers_2019)) +
  xlab("Leaf Trait PC2") + ylab("Flower #")+
  geom_point(size=1.5,alpha=0.25, aes(color=taxon))+scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF"))+
  geom_smooth(method=lm, linewidth=1.25, se=FALSE, color="darkgray")+
  geom_smooth(method=lm, aes(color=taxon),linewidth=1, se=FALSE)


#### Figure 5B,C,D,E,F,G #########

p1<-pc1fruit +graphtheme2
p2<-pc2fruit +graphtheme2
p3<-pc1flw+graphtheme2
p4<-pc2flw+graphtheme2
p5<-pc1biomass+graphtheme2
p6<-pc2biomass+graphtheme2

quartz(width = 10, height = 15) #this opens a window so you get a good size
grid.arrange(p1,p2,p3,p4,p5,p6, nrow=3)