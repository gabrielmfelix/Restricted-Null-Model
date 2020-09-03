# Restricted-Null-Model

Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com.

Authors: Gabriel M. Felix, Rafael B. P. Pinheiro, and Marco A. R. Mello.

E-mail: gabriel.felixf@hotmail.com.  

Published on September 3rd, 2020 (English version).

Run in R version 4.0.2 (2020-06-22) -- "Taking Off Again"

Disclaimer: You may use this script freely for commercial or non-commercial purposes at your own risk. We assume no responsibility or liability for the use of this software, convey no license or title under any patent, copyright, or mask work right to the product. We reserve the right to make changes in the software without notification. We also make no representation or warranty that such application will be suitable for the specified use without further testing or modification. If this script helps you produce any academic work (paper, book, chapter, dissertation etc.), please acknowledge the authors and cite the source.


## Functionality and origin

R code provided in this repository can be used to simulate null interaction matrices which conserve both the modular structure and the marginal totals of a given  interaction matrix.

This is the **restricted null model** used in [Felix et al 2017](https://doi.org/10.1101/236687), [Pinheiro et al 2019](https://doi.org/10.1002/ecy.2796), and [Mello et al 2019](https://doi.org/10.1038/s41559-019-1002-3). It was derived from the [vaznull model](https://doi.org/10.1111/j.0030-1299.2007.15828.x).

The restricted null model is particularly useful for testing for a compound topology, i.e., a modular structure with internally nested modules. Our model allows comparing observed and expected values of nestedness between species of the same module (NODFsm), and between species of different modules (NODFdm). 

Functions to compute NODFsm and NODFdm have already been implemented in the [bipartite package for R](https://cran.r-project.org/web/packages/bipartite/index.html).

"Posterior.Prob" computes pairwise probabilities of interaction among species. These pairwise probabilities will be used by Rest_NUllModel when drawning interactions on the null matrices.


## Posterior.Prob

Computing probabilities of interaction based on a modular structure.

### Arguments
M -> The interaction matrix for which posterior probability will be computed

R.partition and C.partition -> rows and columns partitions 

Prior.Pij -> Method to be used when computed the "a priori" probability of interaction among species i and j. Can be defined as: 

1. "equiprobable", probability of interaction identical to all species  

2. "degreeprob" probability of interaction proportional to overall species degrees

3. "degreeprob.byarea", probability of interaction proportional to species degrees in each matrix area

Conditional.level -> The level to which conditional probability of interaction among species i and j will be conditionated. Can be defined as: 

1. "matrix": conditional probabilities identical in all matrix areas

2. "modules": conditional probabilities differing between areas within and outside modules

3. "areas": a different conditional probability in each matrix area


## RestNullModel

Restricted null model derived from the vaznull model.

### Arguments
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
