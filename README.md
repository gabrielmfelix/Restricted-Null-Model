# Restricted-Null-Model

[Ecological Synthesis Lab](https://marcomellolab.wordpress.com) (SintECO)

Authors: Gabriel M. Felix, Rafael B. P. Pinheiro, and Marco A. R. Mello.

E-mail: gabriel.felixf@hotmail.com.  

Published originaly on September 3rd, 2020 (English version).

Run in R version 4.0.2 (2020-06-22) -- "Taking Off Again"

Disclaimer: You may use this script freely for commercial or non-commercial purposes at your own risk. We assume no responsibility or liability for the use of this software, convey no license or title under any patent, copyright, or mask work right to the product. We reserve the right to make changes in the software without notification. We also make no representation or warranty that such application will be suitable for the specified use without further testing or modification. If this script helps you produce any academic work (paper, book, chapter, dissertation etc.), please acknowledge the authors and cite the source.


## List of files

See further info in the respective sections.

1. CompoundTest.R -> commented script with the steps required to test for a compound topology in an interaction matrix.

2. CompoundTest.pdf -> tutorial in PDF format with the steps required to test for a compound topology in an interaction matrix.

3. CompoundTest.Rmd -> tutorial in RMD format with the steps required to test for a compound topology in an interaction matrix.

4. PosteriorProb.R -> script for calculating interaction probabilities to be used with "RestNullModel.R".

5. RestNullModel.R -> script for generating randomized matrices based on the restricted null model. 

6. net1.txt -> example network with a typical compound topology.


## Functionality and origin

R code provided in this repository can be used to generate randomized matrices that conserve both the modular structure and the marginal totals of an original matrix.

This is the **restricted null model** used in [Felix et al 2017](https://doi.org/10.1101/236687), [Pinheiro 2019](http://hdl.handle.net/1843/33333), [Pinheiro et al 2019](https://doi.org/10.1002/ecy.2796), [Mello et al 2019](https://doi.org/10.1038/s41559-019-1002-3), and [Queiroz et al 2020](https://github.com/marmello77/queiroz-et-al-2020). It was derived from the [vaznull model](https://doi.org/10.1111/j.0030-1299.2007.15828.x). The synthesis presented in these new functions and models, reported in a series of studies, was based on the ideas first proposed by [Lewinsohn et al. 2006](http://doi.wiley.com/10.1111/j.0030-1299.2006.14583.x), [Mello et al. 2009](http://doi.wiley.com/10.1111/j.1365-2656.2009.01567.x), [Flores et al. 2013](http://www.nature.com/doifinder/10.1038/ismej.2012.135), and [Pinheiro et al 2016](https://doi.org/10.1016/j.ijpara.2015.10.002).

Our restricted null model was designed for testing for a compound topology, i.e., a modular network structure with internally nested modules. It allows comparing observed and expected values of nestedness between species of the same module (NODFsm), and between species of different modules (NODFdm). 

The function *nest.smdm* for computing NODFsm and NODFdm has already been implemented in the [package bipartite for R](https://cran.r-project.org/web/packages/bipartite/index.html), as well as the functions *sortmatrix* and *plotmatrix* for drawing matrices in a way that helps visualizing a compound topology.

In this repo, we integrated all those functions aiming to make the analysis of compound topologies easier.

In addition to the functions already implemented in the package bipartite for R, this repo contains 2 new functions and 1 integrative script and tutorial to be used sequentially.


## Instructions

1. If you are fully familiar with R and the publications mentioned here, run the functions separately and experiment with them;

2. If you would like to see how these functions work together, run the commented script "CompoundTest.R";

3. Alternatively, run the tutorial provided in "CompoundTest.Rmd";

4. If you are not so familiar with R, read the tutorial provided in "CompoundTest.pdf".


## (1) PosteriorProb

Computes pairwise probabilities of interaction among species for a matrix with a modular structure. 

### Arguments

1. M -> matrix. The original matrix for which posterior probabilities will be calculated.

2. R.partitions -> vector of integers. It must containd info on row partitions.

3. C.partitions -> vector of integers. It must containd info on column partitions.

4. Prior.Pij -> method for computing "a priori" probabilities of interaction between species i and j. Can be defined as: 

    a. "equiprobable": probability of interaction identical to all species.  

    b. "degreeprob": probability of interaction proportional to overall species degrees.

    c. "degreeprob.byarea": probability of interaction proportional to species degrees in each matrix area - see "areas" in "Conditional.level" for a definition of matrix areas.

5. Conditional.level -> level to which conditional probability of interaction among species i and j will be conditioned. Can be defined as: 

    a. "matrix": conditional probabilities identical in all matrix areas.

    b. "modules": conditional probabilities differing between areas within and outside modules.

    c. "areas": a different set of conditional probabilities for each matrix area. A matrix area is a submatrix M[AB] of M formed by all rows of module A and all columns of module B. If A = B, then M[AB] is a module area, otherwise M[AB] is the area between two modules. Therefore, when Conditional.level =  "areas", each area have their own conditional probabilities of interaction.  


## (2) RestNullModel

Restricted null model derived from the vaznull model. Uses the pairwise probabilities generated by PosteriorProb to draw interactions for the null matrices.

### Arguments

1. M: Matrix -> matrix. The original matrix to be randomized.

2. Pij.Prob -> matrix. Matrix of probabilities with the same dimensions of M, computed by PosteriorProb.

3. Numbernulls -> integer. Number of null matrices to be produced.

4. Print.null -> logical. If simulation progress should be printed. Default is FALSE.

5. allow.degeneration -> logical. If null matrices are allowed to degenerate. Default is FALSE. If TRUE interactions are drawn without assure that all rowns and columns must have at least one interaction each.

6. return.nonrm.species -> logical. If the index of non-removed rows and columns should be returned in the output. Default is TRUE.

7. connectance -> logical. If connectance of the null matrices should be either exactly (TRUE) or aproximately (FALSE) the same as the original matrix. Default is TRUE.

8. byarea- > logical. If interactions should be drawn independently for each matrix area (i.e., in each submatrix M[ij] of M formed by all rows of module i and all columns of module j). Default is FALSE

9. R.partitions -> vector of integers. Partition of rows. Used only if byarea = TRUE.

10. C.partitions -> vector of integers. Partition of columns. Used only if byarea = TRUE.


## (3) CompoundTest 

This script, provided in 2 formats, can be used to test for a compound topology in an interaction network. 

If you are fully confortable using R, run the original script in R format.

If you would prefer a more didatic example, run the script in RMD format.

It integrates the two new functions presented here and other functions that have already been implemented in the package bipartite for R.

The net1.txt file contains an example network with a compound topology, which can be used in a test drive. If you want to test your own network, replace this file with another one formatted in the same way. And don't forget to keep object names consistent.

Follow the instructions given in the script to run a compound topology test.


## Acknowledgements

We thank all colleagues who, together with us, wrote the papers used to lay the ground for this synthesis. Special thanks go to Renata Muylaert, Pavel Dodonov, Alexandre Palaoro, Danilo Muniz, and the StackOverflow community, who helped us improve our scripts and coding skills in general. Last, but no least, we thank Carsten Dormann for incorporating many of our codes into the package bipartite for R.


## Source studies

If you want to understand the background of those new functions before using them, read the following studies. The first three paved the ground for the analysis of compound topologies developed later by our lab.

1. Lewinsohn, T. M., P. Inácio Prado, P. Jordano, J. Bascompte, and J. M. Olesen. 2006. Structure in plant-animal interaction assemblages. Oikos 113: 174–184. Available at: http://doi.wiley.com/10.1111/j.0030-1299.2006.14583.x.

2. Bezerra, E. L. S., I. C. Machado, and M. A. R. Mello. 2009. Pollination networks of oil-flowers: a tiny world within the smallest of all worlds. J. Anim. Ecol. 78: 1096–1101. Available at: http://www.ncbi.nlm.nih.gov/pubmed/19515098.

3. Flores, C. O., S. Valverde, and J. S. Weitz. 2013. Multi-scale structure and geographic drivers of cross-infection within marine bacteria and phages. ISME J. 7: 520–532. Available at: http://www.nature.com/doifinder/10.1038/ismej.2012.135.

4. Pinheiro, R. B. P., G. M. F. Félix, A. V Chaves, G. A. Lacorte, F. R. Santos, É. M. Braga, and M. A. R. Mello. 2016. Trade-offs and resource breadth processes as drivers of performance and specificity in a host–parasite system: a new integrative hypothesis. Int. J. Parasitol. 46: 115–121. Available at: http://www.sciencedirect.com/science/article/pii/S0020751915002933.

5. Felix, G. M., R. B. P. Pinheiro, R. Poulin, B. R. Krasnov, and M. A. R. Mello. 2017. The compound topology of a continent-wide interaction network explained by an integrative hypothesis of specialization. bioRxiv 236687. Available at: https://doi.org/10.1101/236687.

6. Pinheiro, R. B. P. 2019. As topologias de redes de interações ecológicas e suas origens. PhD Thesis, Federal Univesity of Minas Gerais. URL: http://hdl.handle.net/1843/33333. 

7. Pinheiro, R. B. P., G. M. F. Felix, C. F. Dormann, and M. A. R. Mello. 2019. A new model explaining the origin of different topologies in interaction networks. Ecology 100: e02796. Available at: https://doi.org/10.1002/ecy.2796.

8. Mello, M. A. R., G. M. Felix, R. B. P. Pinheiro, R. L. Muylaert, C. Geiselman, S. E. Santana, M. Tschapka, N. Lotfi, F. A. Rodrigues, and R. D. Stevens. 2019. Insights into the assembly rules of a continent-wide multilayer network. Nat. Ecol. Evol. 3: 1525–1532. Available at: https://doi.org/10.1038/s41559-019-1002-3.

