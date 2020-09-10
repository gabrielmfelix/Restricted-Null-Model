#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Gabriel M. Felix, Rafael B. P. Pinheiro, and Marco A. R. Mello.

#### See README for further info: https://github.com/gabrielmfelix/Restricted-Null-Model/blob/master/README.md

#### How to test for a compound topology using the restricted null model

#### Follow the instructions in the sequence and will be able to run a compound topology test.

#### You can also replace the example network with your own network and test it. Just remember to keep the names consistent.


############### SUMMARY ##### 


#1. Preparing the data
#2. Modularity analysis
#3. Nestedness analysis
#4. Restricted null model analysis
#5. Plotting the network
#6. Source studies



############### 1. PREPARING THE DATA ############### 


#Set the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Delete all previous objects
rm(list= ls())

#Clear the console
cat("\014") 

#Load the required packages and functions.
library(bipartite)
source("RestNullModel.R")
source("PosteriorProb.R")

#Load the data
data<-as.matrix(read.table("net1.txt", head=TRUE))

#Visualize the data
data

#Check the data
dim(data)
class(data)

#Visualize the raw matrix
visweb(data)



############### 2. MODULARITY ANALYSIS ############### 


#Compute modularity
Mod <- bipartite::computeModules(data)

#Recover the partitions
Part <- bipartite::module2constraints(Mod)
row.Part <- Part[1:nrow(data)]
col.Part <- Part[(nrow(data)+1):(nrow(data)+ncol(data))]

#Test for the significance of modularity with a Monte Carlo procedure

#Generate randomized matrices
nulls <- nullmodel(data, N=9, method="r2d")

#Calculate the modularity of the randomized matrices
mod.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(mod.nulls, function(x) x@likelihood)

#Calculate the z-score of the randomized distribution
(z <- (Mod@likelihood - mean(like.nulls))/sd(like.nulls))

#Plot the observed modularity value against the distribution of randomized values
plot(density(like.nulls), xlim=c(min((Mod@likelihood), min(like.nulls)), max((Mod@likelihood), max(like.nulls))), 
     main="Observed vs. randomized")
abline(v=(Mod@likelihood), col="red", lwd=2)    

#Estimate the P-value
mean(like.nulls)
sd(like.nulls)
Mod@likelihood
praw <- sum(like.nulls>(Mod@likelihood)) / length(like.nulls)
ifelse(praw > 0.5, 1-praw, praw)



############### 3. NESTEDNESS ANALYSIS ############### 


#Calculate the desired nestedness metric (here WNODA) for the original network.
obs <- unlist(bipartite::nest.smdm(x = data, 
                                   constraints = Part, #Input the modular structured recovered from step 2
                                   weighted = T, #By considering the edge weights, you are choosing WNODA
                                   decreasing = "abund"))

#Check the scores
obs



############### 4. RESTRICTED NULL MODEL ANALYSIS ############### 


#Calculate constrained interaction probabilities considering the network's modular structure
Pij <- PosteriorProb(M = data, 
                     R.partitions = row.Part, C.partitions = col.Part, #Input the modular structured recovered from step 2
                     Prior.Pij = "degreeprob", #Choose the null model
                     Conditional.level = "modules") #Choose the kind of constraints

#Check what those probabilities look like
Pij

#Generate randomized networks with the null model of your choice, considering the interaction probabilities calculated before. 
nulls <- RestNullModel(M = data, 
                       Pij.Prob = Pij, #Recover the probabilities calculated in the previous command
                       Numbernulls = 9, #This step may take long, so start experimenting with low values
                       Print.null = T, allow.degeneration = F, #Choose whether you allow orphan rows and columns to be removed or not
                       return.nonrm.species = F, 
                       connectance = T, byarea = T, 
                       R.partitions = row.Part, C.partitions = col.Part)

#Calculate the same nestedness metric for all randomized networks
null <- sapply(nulls, function(x) bipartite::nest.smdm(x = x, constraints = Part, weighted = T, decreasing = "abund"))
WNODA.null <- unlist(null[3,])
WNODAsm.null <- unlist(null[8,])
WNODAdm.null <- unlist(null[9,])

#Plot the observed nestedness value against the distribution of randomized values
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

#Estimate the P-values

#Nestedness in th entire network
praw.WNODA <- sum(WNODA.null>obs[3]) / length(WNODA.null)
p.WNODA <- ifelse(praw.WNODA > 0.5, 1- praw.WNODA, praw.WNODA)    # P-value
p.WNODA

#Nestedness within the modules
praw.WNODAsm <- sum(WNODAsm.null>obs[8]) / length(WNODAsm.null)
p.WNODAsm <- ifelse(praw.WNODAsm > 0.5, 1- praw.WNODAsm, praw.WNODAsm)    # P-value
p.WNODAsm

#Nestedness between the modules
praw.WNODAdm <- sum(WNODAdm.null>obs[9]) / length(WNODAdm.null)
p.WNODAdm <- ifelse(praw.WNODAdm > 0.5, 1- praw.WNODAdm, praw.WNODAdm)    # P-value
p.WNODAdm


############### 5. PLOTTING THE NETWORK ############### 


par(mfrow = c(1,1))

#Sort the matrix in a way that facilitates visualizing the compound topology
data.comp <- bipartite::sortmatrix(matrix = data, topology = "compound", sort_by = "weights", row_partitions = row.Part, col_partitions = row.Part)

#Assign colors for the modules
modcol <- rainbow((length(unique(Part))), alpha=1)

#Plot the matrix
plotmatrix(data.comp$matrix, 
           row_partitions = data.comp$row_partitions, 
           col_partitions = data.comp$col_partitions, 
           border = T,
           binary = F,
           modules_colors = modcol,
           within_color = modcol, 
           between_color = "lightgrey")



############### 6. SOURCE STUDIES ############### 

#If you want to understand the background of those new analyses before using them, read the following studies. The first three paved the ground for the analysis of compound topologies developed later by our lab.

# 1. Lewinsohn, T. M., P. Inácio Prado, P. Jordano, J. Bascompte, and J. M. Olesen. 2006. Structure in plant-animal interaction assemblages. Oikos 113: 174–184. Available at: http://doi.wiley.com/10.1111/j.0030-1299.2006.14583.x.

# 2. Bezerra, E. L. S., I. C. Machado, and M. A. R. Mello. 2009. Pollination networks of oil-flowers: a tiny world within the smallest of all worlds. J. Anim. Ecol. 78: 1096–1101. Available at: http://www.ncbi.nlm.nih.gov/pubmed/19515098.

# 3. Flores, C. O., S. Valverde, and J. S. Weitz. 2013. Multi-scale structure and geographic drivers of cross-infection within marine bacteria and phages. ISME J. 7: 520–532. Available at: http://www.nature.com/doifinder/10.1038/ismej.2012.135 [Accessed June 9, 2016].

# 4. Pinheiro, R. B. P., G. M. F. Félix, A. V Chaves, G. A. Lacorte, F. R. Santos, É. M. Braga, and M. A. R. Mello. 2016. Trade-offs and resource breadth processes as drivers of performance and specificity in a host–parasite system: a new integrative hypothesis. Int. J. Parasitol. 46: 115–121. Available at: http://www.sciencedirect.com/science/article/pii/S0020751915002933.

# 5. Felix, G. M., R. B. P. Pinheiro, R. Poulin, B. R. Krasnov, and M. A. R. Mello. 2017. The compound topology of a continent-wide interaction network explained by an integrative hypothesis of specialization. bioRxiv 236687. Available at: https://doi.org/10.1101/236687.

# 6. Pinheiro, R. B. P. 2019. As topologias de redes de interações ecológicas e suas origens. PhD Thesis, Federal Univesity of Minas Gerais. URL: http://hdl.handle.net/1843/33333. 

# 7. Pinheiro, R. B. P., G. M. F. Felix, C. F. Dormann, and M. A. R. Mello. 2019. A new model explaining the origin of different topologies in interaction networks. Ecology 100: e02796. Available at: https://doi.org/10.1002/ecy.2796.

# 8. Mello, M. A. R., G. M. Felix, R. B. P. Pinheiro, R. L. Muylaert, C. Geiselman, S. E. Santana, M. Tschapka, N. Lotfi, F. A. Rodrigues, and R. D. Stevens. 2019. Insights into the assembly rules of a continent-wide multilayer network. Nat. Ecol. Evol. 3: 1525–1532. Available at: https://doi.org/10.1038/s41559-019-1002-3.

