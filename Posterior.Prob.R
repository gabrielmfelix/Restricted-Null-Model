#### Computing probability of interactions based on modular structure
#### Felix, G. M., Pinheiro, R. B. P., Poulin, R., Krasnov, B. R., & de Mello, M. A. R. (2017). The compound topology of a continent-wide interaction network explained by an integrative hypothesis of specialization. bioRxiv, 236687.

### M -> The interaction matrix for which posterior probability will be computed

### R.partition and C.partition -> rows and columns partitions 

### Prior.Pij -> Method to be used when computed the "a priori" probability of interaction among species i and j. Can be defined as: 
# (i) "equiprobable", probability of interaction identical to all species  
# (ii) "degreeprob" probability of interaction proportional to overall species degrees
# (iii) "degreeprob.byarea", probability of interaction proportional to species degrees in each matrix area

### Conditional.level -> The level to which conditional probability of interaction among species i and j will be conditionated. Can be defined as: 
# (i)   "matrix": conditional probabilities identical in all matrix areas
# (ii)  "modules": conditional probabilities differing between areas within and outside modules
# (iii) "areas": a different conditional probability in each matrix area

Posterior.prob <- function(M, R.partitions, C.partitions, Prior.Pij, Conditional.level){
  
  # Test of assumptions
  if (!is.matrix(M)){stop("M is not a matrix")}
  if (0 %in% rowSums(M) | 0 %in% colSums(M)) {stop("M is degenerated. There are rows and/or columns withouth interactions in the matrix. Remove then before continue")}
  if (!is.numeric(R.partitions) | !is.numeric(C.partitions)) {stop("Partitions are not numeric")}
  if (length(R.partitions) != nrow(M) | length(C.partitions) != ncol(M)) {stop("Partitions and matrix dimensions have different sizes")}
  if (!(Conditional.level %in% c("matrix","modules","areas"))) {stop("Conditional.level should be 'matrix','modules' or 'areas'")}
  if (Prior.Pij != "degreeprob" & Prior.Pij != "equiprobable" & Prior.Pij != "degreeprob.byarea") {stop("Pij.probs should be 'equiprobable' or 'degreeprob' or 'degreeprob.byarea")}
  
  # M dimensions
  r <- dim(M)[1] # Number of rows
  c <- dim(M)[2] # Number of columns
  array()
  
  # Constructing a array with r rows, c columns and 3 slices. This array represents the modular structure
  # The first slice informs if a given matrix cell M(rc) is within (1) or outside (0) a module
  # The second slice informs to which module the species in the row (r) of a given cell M(rc) belongs
  # The third slice informs to which module the species in the column (c) of a given cell M(rc) belongs 
  
  Matrix.mod <- array(0, dim = c(r, c, 3))
  for (rr in 1:r){
    for (cc in 1:c){
      Matrix.mod[rr,cc,1] <- ifelse(R.partitions[rr] == C.partitions[cc], 1,0)
      Matrix.mod[rr,cc,2] <- R.partitions[rr]
      Matrix.mod[rr,cc,3] <- C.partitions[cc]
    }
  }
  
  # Defining a priori Pij probabilities
  if (Prior.Pij == "equiprobable"){
    Pi <- rep(1 / r, times = r) 
    Pj <- rep(1 / c, times = c)
    Prior.Pij.species <- tcrossprod(Pi, Pj)
  }else if (Prior.Pij == "degreeprob"){
    Pi <- rowSums(M) / sum(rowSums(M)) 
    Pj <- colSums(M) / sum(colSums(M))
    Prior.Pij.species <- tcrossprod(Pi, Pj)
  }else if(Prior.Pij == "degreeprob.byarea"){
    Prior.Pij.species <- M
    RMod <- sort(unique(R.partitions))
    CMod <- sort(unique(C.partitions))
    for (rr in RMod){
      for (cc in CMod){
        M.rr.cc <- M[R.partitions == rr,C.partitions == cc]
        Pi.rr.cc <- rowSums(M.rr.cc) / sum(rowSums(M.rr.cc)) 
        Pj.rr.cc <- colSums(M.rr.cc) / sum(colSums(M.rr.cc))
        Prior.Pij.species[R.partitions == rr, C.partitions == cc] <- tcrossprod(Pi.rr.cc, Pj.rr.cc)
      }
    }
  }

  # Defining a conditional probabilities by area based on species degrees and conectance by area
  
  if (Conditional.level == "matrix"){
    Post.Pij <- Prior.Pij.species
  }else {
    Prior.Pij.area <- matrix(NA,r,c)
    Cond.Pij.area <- matrix(NA,r,c)
    if (Conditional.level == "modules"){
      WMod.prior <- sum(Prior.Pij.species[Matrix.mod[,,1] == 1])
      OMod.prior <- sum(Prior.Pij.species[Matrix.mod[,,1] == 0])
      Prior.Pij.area[Matrix.mod[,,1] == 1] <- WMod.prior
      Prior.Pij.area[Matrix.mod[,,1] == 0] <- OMod.prior
    
      WMod.cond <- sum(M[Matrix.mod[,,1] == 1]) / sum(M)
      OMod.cond <- sum(M[Matrix.mod[,,1] == 0]) / sum(M)
      Cond.Pij.area[Matrix.mod[,,1] == 1] <- WMod.cond
      Cond.Pij.area[Matrix.mod[,,1] == 0] <- OMod.cond
    }else if (Conditional.level == "areas"){
      RMod <- sort(unique(R.partitions))
      CMod <- sort(unique(C.partitions))
      for (rr in RMod){
        for (cc in CMod){
          WArea.prior <- sum(Prior.Pij.species[Matrix.mod[,,2] == rr & Matrix.mod[,,3] == cc])
          Prior.Pij.area[Matrix.mod[,,2] == rr & Matrix.mod[,,3] == cc] <- WArea.prior
        
          WArea.cond <- sum(M[Matrix.mod[,,2] == rr & Matrix.mod[,,3] == cc]) / sum(M)
          Cond.Pij.area[Matrix.mod[,,2] == rr & Matrix.mod[,,3] == cc] <- WArea.cond
        }
      } 
    }
  
    # Adjusting the prior Pij prob by conditional probabilities 
  
    Post.Pij <- Prior.Pij.species * (Cond.Pij.area / Prior.Pij.area)
  }
  
  return(Post.Pij = Post.Pij)
}


