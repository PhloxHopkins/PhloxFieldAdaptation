# PhloxFieldAdaptation
This code is for analyzing field collected data from the reciprocal transplant experiment with Phlox pilosa, P. amoena, and P. deamii from 2018-2020

The first two steps use the raw data from the field experiment to generate estimates for the random effect of source population. These estimates are then used to infer local adaptation across populations through correlations between population effect and distance. 

The third step analyzes the data to look for adaptive divergence between species. The models estimate the interaction between species and garden and use constrasts to test for an adventage of the local species

The four step looks at phenotypic predictors of fitness variation between and within species. It uses the first to principle components of leaf trait variation to ask if the same trait variation predicts fitness differences between species and fitness diferences within a specie across populations. 

The fifth step compares leaf trait variation measured in the greenhouse to trait variation measured in the natural field environment.

The final and sixth step uses a method described by Chong et al 2018 in Evolution Letters to reconstitute selection gradients from a PCA. This generates a matrix of eigenvectors and a vector of regression coeffecients of fitness traits onto the first three PC. The product of this matrix and vector is a vector of selection gradients for traits included in the PCA. 
