# Restricted-Null-Model
The functions provided in this repository are the same used in Felix et al (2017) DOI:https://doi.org/10.1101/236687
These functions can be used to simulate random interaction matrices which conserve both the modular structure and the marginal totals of a given real interaction matrix. The restricted null model is particularly useful when testing the compound topology hypothesis.
Posterior.Prob computes pairwise probabilities of interaction among species. This probabilities will be used for Rest_NUllModel when drawning interactions on the null matrices.

## Posterior.Prob
Computing probability of interactions based on modular structure

M -> The interaction matrix for which posterior probability will be computed

R.partition and C.partition -> rows and columns partitions 

Prior.Pij -> Method to be used when computed the "a priori" probability of interaction among species i and j. Can be defined as: 
(i) "equiprobable", probability of interaction identical to all species  
(ii) "degreeprob" probability of interaction proportional to overall species degrees
(iii) "degreeprob.byarea", probability of interaction proportional to species degrees in each matrix area

Conditional.level -> The level to which conditional probability of interaction among species i and j will be conditionated. Can be defined as: 
(i)   "matrix": conditional probabilities identical in all matrix areas
(ii)  "modules": conditional probabilities differing between areas within and outside modules
(iii) "areas": a different conditional probability in each matrix area

## RestNullModel
Vaznull algorithm of bipartite modified to run the restricted null model

M: Matrix. A matrix. Interaction matrix to be randomized

Pij.Prob: A matrix. Matrix of probabilities, with the same dimensions of M, computed by function "Posterior.Prob"

Numbernulls: Interger. Number of null matrices to be produced

Print.null: Logical. If simulation progress should be printed. Default is FALSE

allow.degeneration: Logical. If null matrices are allowed to degenerate. Default is FALSE

return.nonrm.species: Logical. If the index of non-removed rows and columns should be returned in the output. Default is T

connectance: Logical. If connectance of the null matrices should be either exactly (TRUE) or aproximately (FALSE) the same as the original matrix. Default is T

byarea: Logical. If interactions should be drawn independently in each matrix area. Default is F

R.partitions: Vector of Intergers. Partition of rows. Only applied if byarea = T

C.partitions: Vector of Intergers. Partition of columns. Only applied if byarea = T
