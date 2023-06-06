
#' Condition number for inversion 
#' 
#' General description 
#' Computes the condition number of a matrix using the specified norm
#'
#' Deails
#' "Mathematically", if the condition number is less than infinite, the matrix is invertible. "Numerically", 
#' there are roundoff errors which occur. A large condition number also indicates that a small change in the 
#' coefficient matrix A can lead to larger changes on its inverse. If a matrix is singular, then its condition  
#' number is infinite. A finite large condition number means that the matrix is close to being singular 
#' (Anderson et al. 1994).
#'
#' References
#' 
#' Anderson. E. and ten others (1999) LAPACK Users' Guide. Third Edition. SIAM.
#' Available on-line at https://www.netlib.org/lapack/lug/lapack_lug.html.
#'
#' https://en.wikipedia.org/wiki/Condition_number
#'
#' https://www.wikiwand.com/en/Singular-value_decomposition
#'
#' https://www.mathworks.com/help/matlab/ref/cond.html#mw_9f5f2cf0-2731-4226-8916-05802ee6768a
#'
#'
#' Arguments
#'
#' @A
#'
#' A square numeric matrix
#'
#' @norm 
#'
#' character string ("1", "2", or "Inf"), specifying the type of matrix norm to be computed
#' 
#' "1"
#' specifies the one norm, (maximum absolute column sum)
#'
#' "2"
#' specifies the “spectral norm” or 2-norm, which is the largest singular value (svd) of A.
#' 
#' "Inf"
#' specifies the infinity norm (maximum absolute row sum)
#'
#' @Value
#'
#' Returns a numeric scalar equal to the condition number
#'
#' @example taken from MATLAB (see 4th reference)
#'
#' create a 2 by 2 matrix
#' A <- matrix(c(4.1,9.7,2.8,6.6), byrow=FALSE, nrow=2)
#' > A
#'     [,1] [,2]
#' [1,]  4.1  2.8
#' [2,]  9.7  6.6
#'
#' calculate condition number of A
#' cnd <- condNumber(A, norm = "2")
#' > cnd
#' [1] 1622.999
#'
#' Since the condition number of A is much larger than 1, the matrix is sensitive to the inverse 
#' calculation. To see this, calculate the inverse of A
#' 
#' > solve(A)
#'     [,1] [,2]
#'[1,]  -66   28
#'[2,]   97  -41
#'
#' Now, make a "small change" in the second row of A and calculate the inverse again.
#' 
#' v <- c(1.2, 3.4)
#' A2 <- rbind(A[1, ], v)
#'
#' > A2
#'      [,1]  [,2]
#' [1,] 4.100 2.800
#' [2,] 9.671 6.608
#'
#' calculate the inverse of A2
#'
#' > solve(A2)
#'.          [,1]      [,2]
#' [1,]  472.0000 -200.0000
#' [2,] -690.7857  292.8571
#'
#' draw your own conclusions...

condNumber <- function(A, norm = "2") {

   library(Matrix, quietly = TRUE)

   if (any(is.na(A))) {
      stop("Matrix 'A' contains missing values.")
   }
   if (ncol(A) != nrow(A)) {
      stop("Matrix 'A' must be square.")
   }
   if (inherits(A, "matrix") && !inherits(A, "sparseMatrix")) {
      As <- as(A, "sparseMatrix")
   } else if (inherits(A, "sparseMatrix")) {
      As <- A
   } else {
      stop("Argument 'A' must be a matrix or a sparseMatrix.")
   }

   norms <- c("1", "2", "Inf")
   if (!is.character(norm) || !(norm %in% norms)) {
      stop("Argument 'norm' must be one of '1', '2', or 'Inf' ")
   }

   ## vector of singular values (sval) from SVD decomposition of matrix A
   sval <- svd(A)[["d"]]

   ## definition of the condition number (cnd) based on the type of norm
   if (norm == "1") {
      cnd <- max(colSums(abs(A))) / min(colSums(abs(A)))
   } else if (norm == "2") {
      cnd <- max(sval) / min(sval)
   } else if (norm == "Inf") {
      cnd <- max(abs(sval)) / min(abs(sval))
   }

   return(cnd)
}