# Relatedness & BLUPs
## Supplementary Material 

This repository contains supplementary material for the "Efficient Prediction of Dominance Effects in Large Breeding Populations: Unveiling the Power of Mixed-Model Methodology" research paper. This supplementary material provides a collection of additional R functions to enhance the research pipeline. The provided functions facilitate breeding program simulations and genetic evaluations for quantitative traits.

## Features

### Breeding Program Simulations

*	Simulate breeding programs for quantitative traits based on the infinitesimal model (Fisher
1918; Bulmer 1980).
*	Incorporate both additive and dominance effects to model the genetic architecture of traits.
*	Utilize stochastic simulation techniques to mimic realistic breeding population dynamics.

### Genetic Relatedness Matrices
Generate the dominance- and family-relatedness matrices using two approaches:
*	Jacquard's nine condensed coefficients of identity (Jacquard, 1974) which quantify accuratelly the genetic relationships among individuals.
*	Coancestry coefficients based on the approach by Cockerham (1954).

### Predictions
Get the best linear unbiased estimate (BLUE) and best linear unbiased predictor (BLUP) by solving the mixed model equations under an animal model (Henderson, 1950).

## Usage
To use the functions, follow these steps:
* Install R programming language (https://www.r-project.org/) and ensure it is properly set up.
* Clone or download this repository.
* Load the required functions into your R environment using the appropriate import statements.
Refer to the documentation and examples provided within each function to understand the input requirements and usage guidelines.

## References
Bulmer, M. G. (1980). The mathematical theory of quantitative genetics. Oxford Univ. Press, Oxford, UK.   

Cockerham, C.C. (1954). An extension of the concept of partitioning hereditary variance for analysis of covariances among relatives when epistasis is present. Genetics, 39, 859–882.

Fisher, R. A. (1918). The correlation between relatives on the supposition of Mendelian inheritance. Trans. R. Soc. Edin. 52, 399-433.

Henderson, C.R. (1950). Estimation of Genetic Parameters. Annals of Mathematical Statistics, 21, 309-310.

Jacquard, A. (1974). The Genetic Structure of Populations. Springer-Verlag, New York.

Please ensure to appropriately cite the relevant publications and authors mentioned in the description.
Feel free to explore the functions and modify them to suit your specific research needs.
If you have any questions or encounter any issues, please open an issue on this repository.
Note: You can provide additional information such as installation instructions, examples, and a license file within the repository to make it more comprehensive and user-friendly.








