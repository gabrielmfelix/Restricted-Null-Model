setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cat("\014")  

# Teste com net1

net1 <- read.delim("net1.txt", row.names = 1)
net1

library(bipartite)

net1 = read.delim("net1.txt", row.names=1)
net1

net1.lpa = computeModules(net1, method = "Beckett")
net1.lpa
printoutModuleInformation(net1.lpa)
modulos.lpa.lista = listModuleInformation(net1.lpa)
source("Extracting_modules_from_bipartite.R")
modules.lpa.dataframe = modules_from_bipartite(modulos.lpa.lista)
colnames(modules.lpa.dataframe$Rows_modules)=c("nodes", "modules")
colnames(modules.lpa.dataframe$Cols_modules)=c("nodes", "modules")
membership.lpa = rbind(modules.lpa.dataframe$Rows_modules, modules.lpa.dataframe$Cols_modules)
membership.lpa
write.table(membership.lpa, "partititions.txt", sep="\t", quote=F)

rowspart = read.delim("partitions.rows.txt", header = T)
colspart = read.delim("partitions.cols.txt", header = T)

rowspart2 <- rowspart$rows
colspart2 <- colspart$cols

dim(net1)
length(rowspart2)
length(colspart2)


net1.nes <- sortmatrix(net1, topology = "nested")
net1.nes

net1.com <- sortmatrix(net1, topology = "compound",
                          row_partitions = rowspart2,
                          col_partitions = colspart2)
net1.com

visweb(net1)

modcol <- c("purple", "orange", "yellow", "blue")

plotmatrix(net1.nes, binary = F, base_color = "blue")

plotmatrix(net1.com, binary = F, border = T, 
           modules_colors = modcol,
           within_color = modcol)


mod <- module2constraints(net1.lpa)

nest.smdm(net1, constraints = mod, 
          weighted = T, 
          decreasing = "abund", 
          sort = T)



# Teste com matriz gerada automaticamente 

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



