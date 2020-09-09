setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list= ls())
cat("\014") 

library(bipartite)
source("RestNullModel.R")
source("PosteriorProb.R")
source("ModulesFromBipartite.R")


##### Teste com net1 #####

net1 <- as.matrix(read.delim("net1.txt", row.names = 1))
net1
class(net1)
visweb(net1)

net1.mod = computeModules(net1, method = "Beckett")
net1.mod@likelihood
printoutModuleInformation(net1.mod)
modulos.lpa.lista = listModuleInformation(net1.mod)
modules.lpa.dataframe = modules_from_bipartite(modulos.lpa.lista)
colnames(modules.lpa.dataframe$Rows_modules)=c("nodes", "modules")
colnames(modules.lpa.dataframe$Cols_modules)=c("nodes", "modules")
membership.lpa = rbind(modules.lpa.dataframe$Rows_modules, modules.lpa.dataframe$Cols_modules)
membership.lpa

rowspart = modules.lpa.dataframe$Rows_modules
colspart = modules.lpa.dataframe$Cols_modules

rowspart2 <- rowspart$modules
class(rowspart2)
colspart2 <- colspart$modules
class(colspart2)

dim(net1)
length(rowspart2)
length(colspart2)

net1.com <- sortmatrix(net1, topology = "compound",
                          row_partitions = rowspart2,
                          col_partitions = colspart2)
net1.com

modcol <- rainbow((length(unique(membership.lpa$modules))), alpha=1)
plotmatrix(net1.com, binary = F, border = T, 
           modules_colors = modcol,
           within_color = modcol)

mod <- module2constraints(net1.lpa)
nest.smdm(net1, constraints = mod, 
          weighted = T, 
          decreasing = "abund", 
          sort = T)

probabilities <- PosteriorProb(net1,
                              rowspart2,
                              colspart2,
                              Prior.Pij = "degreeprob.byarea",
                              Conditional.level = "modules")
probabilities


randomized <- RestNullModel(net1, probabilities, Numbernulls = 9, 
                         Print.null = T, allow.degeneration = T,
                         return.nonrm.species = T, connectance = T,
                         byarea = F, R.partitions = F, C.partitions = F)

class(randomized)
head(randomized)


nestedness.randomized <- sapply(randomized,
                                nest.smdm)

nulls <- nullmodel(data, N=9, method="r2d")
mod.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(mod.nulls, function(x) x@likelihood)
(z <- (mod@likelihood - mean(like.nulls))/sd(like.nulls))

#Plot the observed value against the distribution of randomized values
plot(density(like.nulls), xlim=c(min((mod@likelihood), min(like.nulls)), max((mod@likelihood), max(like.nulls))), 
     main="Observed vs. randomized")
abline(v=(mod@likelihood), col="red", lwd=2)    

#Estimate the P-value
mean(like.nulls)
sd(like.nulls)
mod@likelihood
praw <- sum(like.nulls>(mod@likelihood)) / length(like.nulls)
ifelse(praw > 0.5, 1-praw, praw)



##### Teste com matriz gerada automaticamente #####

set.seed(3)

matriz <- matrix(
  (sample(0:1000, 144, replace = T)),
  nrow = 12, ncol = 12)
matriz


R.partition <- c(sample(1:3, nrow(matriz), replace = T))
R.partition
C.partition <- c(sample(1:3, ncol(matriz), replace = T))
C.partition


Pij.Prob.equ <- PosteriorProb(matriz,
               R.partition,
               C.partition,
               Prior.Pij = "equiprobable",
               Conditional.level = "modules")
Pij.Prob.equ


Pij.Prob.deg <- PosteriorProb(matriz,
                             R.partition,
                             C.partition,
                             Prior.Pij = "degreeprob",
                             Conditional.level = "modules")
Pij.Prob.deg


Pij.Prob.are <- PosteriorProb(matriz,
                               R.partition,
                               C.partition,
                               Prior.Pij = "degreeprob.byarea",
                               Conditional.level = "modules")
Pij.Prob.are


modelos <- RestNullModel(matriz, Pij.Prob.are, Numbernulls = 100, 
                         Print.null = T, allow.degeneration = T,
                         return.nonrm.species = T, connectance = T,
                         byarea = F, R.partitions = F, C.partitions = F)



