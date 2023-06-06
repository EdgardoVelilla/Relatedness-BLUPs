
#' rmvChol  function
#'                                                                                                                            
#' Author: Edgardo Velilla P.                                                                                                 
#' email{edgardo.velilla@cmpc.cl}                                                                                             
#' Created: 06-May-2022                                                                                                       
#' Modified: 31-Mar-2023    
#' 
#' General description:                                                                                                       
#' Generates n samples from a d-dimensional multivariate normal (MVN) 
#' distribution ~ Nd(mu, C).                                                                                                            
#'                                                                                                                            
#' Algorithm to generate multivariate normal samples                                                                          
#'                                                                                                                            
#' 1. Generate an nxd matrix Z such that Zij ~ N(0,1) .                                                          
#' 2. Compute the factorization Gu = Q'Q = LL', where Gu = Var(u)= sigma2_u*G 
#' and u ~ MVN(0, sigma2_u*G) (Gentler 2017) .                                                                                                                                                                                                                                       
#' 3. Apply the transformation X = ZQ + Jmu'. Each row of X is a random variate 
#' from the Nd ~ (mu, Gu) distribution (Searle 1982)).             
#'                                                                                                                            
#' Algorithm to generate multivariate normal samples whose distribution 
#' "exactly matching" to a pre-specified parametrization (Ripley 1987).
#'                                                                                                                            
#' 1. Proceede as before (step 1)                                                                                             
#' 2. Then, subtracting the sample mean of Z, such as Z* = Z - mean(Z).                                                        
#' 3. Computing the Cholesky decomposition of Z*. Thas is,  chol(Z*) = L*, 
#'    where chol() is the Cholesky decomposition.                                                                         
#' 4. Calculating z0 = [(L*) - 1]Z* (z0 should have sample mean 0 and identity 
#'    sample covariance).                                
#' 5. Calculating X = [(L*) - 1]Z* + Jmu'. X should be matched with the desired
#'    parametrization.                                 
#'                                                                                                                            
#' References:
#' 
#' Chen, Y., Davis, T.A., Hager, W.W., & Rajamanickam, S. (2008). Algorithm 887
#' : CHOLMOD, Supernodal Sparse Cholesky Factorization and Update/Downdate. ACM 
#' Trans. Math. Softw., 35, 22:1-22:14.                                                                                                                   
#'                                                                                                                            
#' Gentle J. (2017). Matrix algebra: theory, computations and applications in 
#' statistics, 2nd edn. Springer, New York.        
#'                                                                                                                            
#' B. D. Ripley (1987) Sthochastic Simulation. Wiley.                                                                         
#'                                                                                                                            
#' Searle, S. R. (1982). Matrix Algebra Useful for Statistics. John Wiley, New
#' York, NY.  
#'                                     
#'                                                                                                                                                                                                                                                       
#' Arguments                                                                                                                  
#'                                                                                                                            
#' @param n                                                                                                                   
#'                                                                                                                            
#' Integer number that defines the sample size.                                                                               
#'                                                                                                                            
#' @param mu                                                                                                                  
#'                                                                                                                            
#' Numerical vector for the expected mean with length equal to the d-dimension
#' of G. The default is a vector of zeros.        
#'                                                                                                                            
#' @param C                                                                                                                   
#'                                                                                                                            
#' A symmetric positive definite matrix of order dxd.                                                                         
#'                                                                                                                            
#' @param adj                                                                                                                 
#'                                                                                                                            
#' Logical indicating if the sampling distribution should be exactly matching 
#' with the desired parametrization. The default is @adj = FALSE.                                                                                                                
#'                                                                                                                            
#' @return                                                                                                                    
#'                                                                                                                            
#' A list with the following components:                                                                                      
#'                                                                                                                            
#' MVN.dt                                                                                                                     
#'                                                                                                                            
#' A data.table with "wide format", i.e. dim(MVN.dt)= (nxd).                                                                   
#'                                                                                                                            
#' MVN.col                                                                                                                    
#'                                                                                                                            
#' A data.table with "long format", i.e. dim(MVN.col)= (n*dx2)                                                                
#'                                                                                                                                                                                                                                              
#'                                                                                                                            
#' A forestry @example                                                                                                         
#'                                                                                                                            
#' Simulation a bivariate distribution between DBH (cm) and Total Height (m) 
#' for E. globulus at 3.8 years old                 
#'                                                                                                                            
#' Assumptions:                                                                                                               
#' var.DBH <- 2.34 # variance assumed for DBH (cm)                                                                            
#' var.height <- 2.62 # variance assumed for Total Height (m)                                                                 
#' var <- c(var.DBH, var.height)                                                                                              
#' mu <- c(10.3,117.7) # mean vector assumed for DBH and Total Height                                                         
#' cor_DBH.TH <- 0.78 # correlation assumed between DBH and Height                                                            
#'                                                                                                                            
#' Hence, the correlation matrix between DBH & Total Height is                                                                
#'  C <- Matrix(c( # DBH      Height                                                                                              
#'                   1.0,     0.7841,                                                                                             
#'                   0.7841,  1.0), 2, 2, sparse = TRUE)	                                                                      
#'                                                                                                                            
#' Rescaling the correlation matrix (C) by pre- and post-multiplying by a 
#' diagonal matrix that contains the standard deviations for DBH & Height  
#'                                                                                               
#' D1 <- .symDiagonal(C@Dim[1], sqrt(variances)) # diagonal matrix for standard
#' deviations of DBH and Total Height        
#'     
#' G <- D1%*%C%*%D1 # G is the covariance matrix for DBH and Total Height                                                     
#'                                                                                                                            
#' > G                                                                                                                        
#' [1,] 2.340000 1.941467                                                                                                     
#' [2,] 1.941467 2.620000                                                                                                     
#'                                                                                                                            
#' set.seed(123)                                                                                                              
#' DBH.Height_sim <- rmvChol(100, mu, G, adj=FALSE) # simulation of 100 trees                                                 
#' Check variance-covariance matrix for this particular realization                                                           
#'                                                                                                                            
#' > var(DBH.Heightsim)                                                                                                       
#'            DBH   Height                                                                                                    
#' DBH    1.949765 1.550507                                                                                                   
#' Height 1.550507 2.174348                                                                                                   
#'                                                                                                                            
# a little bit different from the varianza-covarianza matrix (G) assumed...?,
#' i.e.  
#'                                           
# > G                                                                                                                         
# [1,] 2.340000 1.941467                                                                                                      
# [2,] 1.941467 2.620000                                                                                                      
#'                                                                                                                            
#' Now try with option adj = TRUE                                                                                             
#' set.seed(123)                                                                                                              
#' DBH.Height_sim <- rmvChol(100, mu, G, adj = TRUE) # simulation of 100 trees                                                
#' check G matrix again...                                                                                                    
#'                                                                                                                            
#' > var(DBH.Height_sim)                                                                                                      
#'            DBH   Height                                                                                                    
#' DBH    2.340000 1.941467                                                                                                   
#' Height 1.941467 2.62000                                                                                                    

rmvChol <- function(n,
                    mu = NULL,
                    C,
                    adj = FALSE) {

  library(Matrix, quietly = TRUE)
  library(data.table, quietly = TRUE)
  source("matrixMethods.r")

  if (inherits(C, "matrix")) {
    C <- Matrix(C, sparse = TRUE, doDiag = FALSE)
  } # -> sparse symm. "dsCMatrix"
  if (!is.null(mu)) {
    d <- length(mu)
    if (any(is.na(match(dim(C), d)))) {
      stop("non-conformable arguments...!")
    }
  } else {
    mu <- rep(0L, dim(C)[1])
    d <- length(mu)
  }
  if (adj & n < dim(C)[1] + 1L) {
    stop(paste0(
      "the sample size must be >= ",
      dim(C)[1] + 1L, " to adjust the empirical distribution"
    ))
  }
  if (is.null(C@Dimnames[[1]])) {
    dimnames(C) <- list(
      seq(1L, nrow(C)),
      seq(1L, ncol(C))
    )
  }

  # create sparse supernodal (LL') Cholesky factorization (Chen et al. 2008)
  CHM <- Cholesky(C, perm = FALSE, super = TRUE)
  Lp <- t(as(CHM, "Matrix"))
  Z <- Matrix(rnorm(n * d, mean = 0, sd = 1),
    nrow = n, 
    ncol = d, 
    sparse = TRUE
  )
  MVN <- Z %*% Lp + Matrix(mu, n, d, byrow = TRUE)
  MVN@Dimnames[[2]] <- C@Dimnames[[1]]

  if (adj) {
    # subtracting the sample mean of z
    Z <- scale(Z, TRUE, FALSE)
    G <- as(cov(Z), "sparseMatrix") # ->  class 'dgCMatrix')
    # Cholesky decomposition works on symmetric matrices...!
    G <- forceSymmetric(new(
      "dgCMatrix",
      i = G@i,
      p = G@p,
      Dim = G@Dim,
      x = G@x,
      Dimnames = list(
        G@Dimnames[[1]],
        G@Dimnames[[2]]
      )
    ))
    G <- drop0(G,
      tol = 1e-15,
      is.Csparse = NA
    )

    CHM0 <- Cholesky(G, perm = FALSE, super = TRUE)
    Lp0 <- t(as(CHM0, "Matrix")) # (dxd)
    # Z0 should have sample mean 0 and identity sample covariance
    Z0 <- Z %*% solve(Lp0)
    # MVN must have a sample with the desired mean & covariance
    MVN <- Z0 %*% Lp + Matrix(mu, n, d, byrow = TRUE)
    MVN@Dimnames[[2]] <- C@Dimnames[[1]]
    MVN <- as.data.table(as.matrix(MVN))
    MVN.col <- copy(MVN)[, ii := seq(1, nrow(MVN))]
    MVN.col <- data.table::melt(MVN.col,
      id.vars = c("ii"),
      measures.vars = c(colnames(MVN.col))
    )
    MVN.col[, c("ii") := NULL]
    setnames(MVN.col, c("level", "value"))
    return(list(
      # wide format of data.table MVN
      MVN.dt = MVN,
      # long format of data.table MVN
      MVN.col = MVN.col
    ))
  } else {
    MVN <- as.data.table(as.matrix(MVN))
    MVN.col <- copy(MVN)[, ii := seq(1, nrow(MVN))]
    MVN.col <- data.table::melt(MVN.col,
      id.vars = c("ii"),
      measures.vars = c(colnames(MVN.col))
    )
    MVN.col[, c("ii") := NULL]
    setnames(MVN.col, c("level", "value"))
    return(list(
      MVN.dt = MVN,
      MVN.col = MVN.col
    ))
  }
}
