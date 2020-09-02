#### Vaznull algorithm of bipartite modified to run the restricted null model
#### Felix, G. M., Pinheiro, R. B. P., Poulin, R., Krasnov, B. R., & de Mello, M. A. R. (2017). The compound topology of a continent-wide interaction network explained by an integrative hypothesis of specialization. bioRxiv, 236687.

### Arguments

# M: Matrix. A matrix. Interaction matrix to be randomized

# Pij.Prob: A matrix. Matrix of probabilities, with the same dimensions of M, computed by function "Posterior.Prob"

# Numbernulls: Interger. Number of null matrices to be produced

# Print.null: Logical. If simulation progress should be printed. Default is FALSE

# allow.degeneration: Logical. If null matrices are allowed to degenerate. Default is FALSE

# return.nonrm.species: Logical. If the index of non-removed rows and columns should be returned in the output. Default is T

# connectance: Logical. If connectance of the null matrices should be either exactly (TRUE) or aproximately (FALSE) the same as the original matrix. Default is T

# byarea: Logical. If interactions should be drawn independently in each matrix area. Default is F

# R.partitions: Vector of Intergers. Partition of rows. Only applied if byarea = T

# C.partitions: Vector of Intergers. Partition of columns. Only applied if byarea = T

rest_null <- function(M, Pij.Prob, Numbernulls, Print.null = F, allow.degeneration = F, 
                              return.nonrm.species = T, connectance = T, byarea = F, R.partitions = NULL, C.partitions = NULL){
  
  ### Test of assumptions
  
  if (!is.matrix(M)){stop("M is not a matrix")}
  if (0 %in% rowSums(M) | 0 %in% colSums(M)) {stop("M is degenerated")}
  
  if (!is.matrix(Pij.Prob)){stop("Pij is not a matrix")}
  if (T %in% c(Pij.Prob < 0)){stop("Pij must contain only numbers >= 0")}
  
  if (nrow(M) != nrow(Pij.Prob) | ncol(M) != ncol(Pij.Prob)){stop("Dimensions of M and Pij.Prob should be identical")}
  
  if (byarea == T){
    if(is.null(C.partitions) | is.null(R.partitions)){stop("Partitions missing")}
    if (length(unique(c(length(R.partitions),nrow(M),nrow(Pij.Prob)))) != 1){stop("The number of elements of R.partition should be the same as the number of rows of M and Pij.prob")}
    if (length(unique(c(length(C.partitions),ncol(M),ncol(Pij.Prob)))) != 1){stop("The number of elements of C.partition should be the same as the number of column of M and Pij.prob")}
    if(!identical(unique(R.partitions), unique(C.partitions))){stop("The number and labels of modules in R.partition and C.partition must be the same")}
  }  
  
  if (Numbernulls <= 0 | !is.numeric(Numbernulls)) {stop("Numbernulls should be a number > 0")}
  if (!is.logical(connectance)){stop("connectance should be logical (T or F)")}
  if (!is.logical(allow.degeneration)){stop("allow.degeneration should be logical (T or F)")}
  if (!is.logical(return.nonrm.species)){stop("return.nonrm.species should be logical (T or F)")}
  if (!is.logical(byarea)){stop("byarea should be logical (T or F)")}
  
  ### M dimensions
  r <- dim(M)[1] # Number of rows
  c <- dim(M)[2] # Number of collums
  
  ### Constructing a array with r rows, c columns and 2 slices. This array represents the matrix area structure
  if (byarea == T){
    Matrix.area <- array(0, dim = c(r, c, 2))
    for (rr in 1:r){
      for (cc in 1:c){
        Matrix.area[rr,cc,1] <- R.partitions[rr]
        Matrix.area[rr,cc,2] <- C.partitions[cc]
      }
    }
    
  }else if (byarea == F){
    ## Assigning all rows and columns to the same partition in order to run the code bellow
    Matrix.area <- array(1, dim = c(r, c, 2))
    R.partitions <- rep(1, nrow(M))
    C.partitions <- rep(1, ncol(M))
    
  }
  
  ### Null model simulation
  
  NullMatrices <- list() # list where the null matrices will be storage
  length(NullMatrices) <- Numbernulls #assigning the number of null matrices to be saved in NullMatrices 
  
  ## Drawing interaction in each null matrix 
  for (nn in 1:Numbernulls){
      
    R.part <- sort(unique(as.vector(Matrix.area[,,1])))
    C.part <- sort(unique(as.vector(Matrix.area[,,2])))
    finalmat <- matrix(NA, r, c)
    
    for (R.p in R.part){
      for (C.p in C.part){
          
        M.a <- as.matrix(M[R.partitions == R.p, C.partitions == C.p])
        Pij.a <- Pij.Prob[R.partitions == R.p, C.partitions == C.p]
          
        r.a <- dim(M.a)[1]
        c.a <- dim(M.a)[2]
        
        P.a <- P1.a <- Pij.a
        finalmat.a <- matrix(0, r.a, c.a)
        
        if(allow.degeneration == F & R.p == C.p){
          
          ## Ensuring that the dimensions of the null matrix will be the same of the original matrix
          
          D.int.finalmat.a <- 0 # The number of rows + columns occupied of the null matrix 
          while (D.int.finalmat.a < sum(dim(M.a))) { # While the dimensions of the null matrix was smaller then the original matrix, keep going
            sel <- sample(1:length(M.a), 1, prob = P.a) # Sample an cell of M.a with probability P.a
            selc <- floor((sel - 1)/(dim(M.a)[1])) + 1 # Recovering column and 
            selr <- ((sel - 1)%%dim(M.a)[1]) + 1 # row of the cell sampled
            if (sum(finalmat.a[, selc]) == 0 | sum(finalmat.a[selr,]) == 0) { # Checking if row or column of the sampled cell is empty
              finalmat.a[sel] <- 1 
              P.a[sel] <- 0
            }
          D.int.finalmat.a <- sum(rowSums(finalmat.a) > 0) + sum(colSums(finalmat.a) > 0) # Setting the new number of dimensions occupied
          }
          # When the number of occupied dimensions of the null matrix was the same as the original matrix, continue
        }
        
        conn.remain <- sum(M.a > 0) - sum(finalmat.a > 0) # The number of cells remaining to be occupied to mantain the original connectance
          
          if (conn.remain > 0) {
              if(connectance == T){
                if (length(which(finalmat.a == 0)) == 1) {
                  add <- which(finalmat.a == 0)
                } else {
                add <- sample(which(finalmat.a == 0), conn.remain, 
                              prob = P1.a[finalmat.a == 0], replace = F)
                }
              }else {
                add <- sample(1:length(finalmat.a), conn.remain, 
                              prob = P1.a, replace = T)
              }
            for (add1 in add){
              finalmat.a[add1] <- finalmat.a[add1] + 1
            }
          }
          
          ### Checking if there are still interactions to be drawn. If applicable, draw.
          int.remain <- (sum(M.a) - sum(finalmat.a))
          if (int.remain > 0) { 
            add <- sample(which(finalmat.a > 0), int.remain, prob = P1.a[which(finalmat.a >0)], replace = T)
            finalmat.a[as.numeric(names(table(add)))] <- finalmat.a[as.numeric(names(table(add)))] + (table(add))
          }
      
        finalmat[R.partitions == R.p, C.partitions == C.p] <- finalmat.a
      }
    }
    
    # Saving outputs
    R2keep <- which(rowSums(finalmat) != 0)
    C2keep <- which(colSums(finalmat) != 0)
    finalmat2 <- finalmat[R2keep,C2keep]
    if (return.nonrm.species == T){
      NullMatrices[[nn]] = list(NullMatrix = finalmat2, RKeep = R2keep, CKeep = C2keep)
    }else if(return.nonrm.species == F){
      NullMatrices[[nn]] = finalmat2
    }
    if (Print.null == T){print(nn)}
  }
  return(NullMatrices = NullMatrices)
}




