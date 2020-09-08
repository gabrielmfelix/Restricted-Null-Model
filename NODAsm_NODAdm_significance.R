##################################################################################################################
#                                                                                                                # 
#   SCRIPT FOR ESTIMATING P-VALUES FOR NESTEDNESS AT DIFFERENT MATRIX SCALES BY USING THE RESTRICTED NULL MODEL  #
#                                                                                                                # 
##################################################################################################################

##### Ecological Synthesis Lab (SintECO)
##### https://marcomellolab.wordpress.com
##### Authors: Gabriel Félix, Rafael Pinheiro & Marco Mello
##### E-mail: gabrielfelixmf@gmail.com
##### How to cite: Gabriel Félix, Rafael Pinheiro & Marco Mello. 2020. How to estimate P-values of nestedness at different matrix scales by using the restricted null model plus Monte Carlo procedures. Ecological Synthesis Lab at the University of São Paulo, Brazil.
##### Published on September 07th, 2020 (English version).
##### Run in R version 4.0.2 (2020-06-22) -- "Taking Off Again"

##### Disclaimer: You may use this script freely for non-comercial purposes at your own risk. We assume no responsibility or liability for the use of this software, convey no license or title under any patent, copyright, or mask work right to the product. We reserve the right to make changes in the software without notification. We also make no representation or warranty that such application will be suitable for the specified use without further testing or modification. If this script helps you produce any academic work (paper, book, chapter, dissertation etc.), please acknowledge the authors and cite the source.


#############################################################


####SUMMARY##################################################


#1. Preparation
#2. Procedure

#############################################################


####1. PREPARATION##########################################

#Set the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Delete all previous objects
rm(list= ls())


#Load the packages.
library(bipartite)
source("Rest_NullModel.R")
source("Posterior.Prob.R")

##############################################################


####2. P-value for one network (Dormann)####

#Create the object to be analyzed
data<-as.matrix(read.table("net1.txt", head=TRUE))
Modules <- rep(1:3,each = 10) # Defining an arbitrary partition for the rows and columns. 

#Calculate the desired metric (WNODA) for the original network (observed)
obs <- unlist(nest.smdm(x = data, constraints = rep(Modules,times = 2), weighted = T, decreasing = "abund"))

#Create randomized networks using a null model. 
#(It may take long, depending on your computer's processing power)
#The paramater "method" sets the null model and the parameter "N" sets the number of permutations
Pij <- PosteriorProb(M = data, R.partitions = Modules, C.partitions = Modules, Prior.Pij = "degreeprob", Conditional.level = "modules")
nulls <- RestNullModel(M = data, Pij.Prob = Pij, Numbernulls = 100, Print.null = T, allow.degeneration = F, return.nonrm.species = F, connectance = T,byarea = T, R.partitions = Modules, C.partitions = Modules)

#Calculate the same metric for all randomized networks
null <- sapply(nulls, function(x) nest.smdm(x = x, constraints = rep(Modules, times = 2), weighted = T, decreasing = "abund"))
WNODA.null <- unlist(null[3,])
WNODAsm.null <- unlist(null[8,])
WNODAdm.null <- unlist(null[9,])

#Plot the observed value against the distribution of randomized values
par(mfrow = c(1,3))
plot(density(WNODA.null), xlim=c(min(obs[3], min(WNODA.null)), max(obs[3], max(WNODA.null))), 
     main="Observed vs. randomized", xlab = "WNODA matrix")
abline(v=obs[3], col="red", lwd=2)    
plot(density(WNODAsm.null), xlim=c(min(obs[8], min(WNODAsm.null)), max(obs[8], max(WNODAsm.null))), 
     main="Observed vs. randomized", xlab = "WNODAsm matrix")
abline(v=obs[8], col="red", lwd=2)    
plot(density(WNODAdm.null), xlim=c(min(obs[9], min(WNODAdm.null)), max(obs[9], max(WNODAdm.null))), 
     main="Observed vs. randomized", xlab = "WNODAdm matrix")
abline(v=obs[9], col="red", lwd=2)    

#Estimate the P-value
praw.WNODA <- sum(WNODA.null>obs[3]) / length(WNODA.null)
ifelse(praw.WNODA > 0.5, 1- praw.WNODA, praw.WNODA)    # P-value
praw.WNODAsm <- sum(WNODAsm.null>obs[8]) / length(WNODAsm.null)
ifelse(praw.WNODAsm > 0.5, 1- praw.WNODAsm, praw.WNODAsm)    # P-value
praw.WNODAdm <- sum(WNODAdm.null>obs[9]) / length(WNODAdm.null)
ifelse(praw.WNODAdm > 0.5, 1- praw.WNODAdm, praw.WNODAdm)    # P-value

