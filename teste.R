setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list= ls())
cat("\014")  


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
                               Prior.Pij = "degreeprob.byarea", #não está dando certo com esta opção. Cheque, por favor.
                               Conditional.level = "modules")
Pij.Prob.are
