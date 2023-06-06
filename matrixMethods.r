requireNamespace("MatrixModels", quietly = TRUE)
library(Matrix, quietly = TRUE)
library(data.table, quietly = TRUE)
library(pedigree, quietly = TRUE)

#' function solve.CHMperm()
#'
#' Solve the system Ax = b with A sparse positive definite matrix via Cholesky
#' factorization.
#'
#' Details
#'
#' A is a sparse symmetric positive definite matrix of order (qxq). If the
#' goal is to get the solutions for the system, Ax = b. Then, b represents the
#' right-hand side (RHS) of the Mixed Model Equations (MME). Conversely, if the
#' purpose is to create the inverse of the matrix A, b represents an identity
#' matrix of order (qxq), and its value is ommit (b = NULL) in the argument
#' (see usage).
#'
#'                      Sparse Cholesky factorization
#'
#' If A is sparse and positive definite, it is usually factored as
#'
#'                              A = PLL'P'
#' where
#' L is a lower triangular with positive diagonal elements. Thus, L is called
#' the Cholesky factor of A, P a permutation matrix, L' and P' are transpose
#' matrices of L and P.
#'
#' Cholesky() function implemented in the 'Matrix' package uses the approximate
#' minimum degree algorithm to produce large blocks of zeros in the matrix
#' (Davis 1996, 2006).
#'
#' Interpretation: we permute the rows and columns of A and Cholesky factor
#' so that
#'                              P'AP = LL'
#' where
#' LL' is a sparse supernodal Cholesky factorization (Chen et al. 2008).
#'
#'              Solving sparse positive definite equations
#'
#' Solve Ax = b with A is a sparse positive definite matrix via factorization
#' A = PLL'P' (Boyd and Vandenberghe 2018; Vandenberghe 2021).
#'
#'
#'  Algorithm
#'  1. compute sparse Cholesky supernodal factorization (Chen et al. 2008):
#'     A = PLL'P'
#'  2. Permute right-hand side (RHS): c: = P'b
#'  3. solve Ly = c by forward substitution
#'  4. solve L'z = y by backward substitution
#'  5. permute solution: x: = Pz
#'
#' References:
#'
#' Chen, Y., Davis, T.A., Hager, W.W., & Rajamanickam, S. (2008). Algorithm 887
#' : CHOLMOD, Supernodal Sparse Cholesky Factorization and Update/Downdate. ACM
#' Trans. Math. Softw., 35, 22:1-22:14.
#'
#' Timothy A. Davis (2006). Direct Methods for Sparse Linear Systems, SIAM
#' Series “Fundamentals of Algorithms”. p. 21-22.
#'
#' Tim Davis (1996). An approximate minimal degree ordering algorithm, SIAM
#' J. Matrix Analysis and Applications, 17, 4, 886–905.
#'
#' Lecture notes from ECE133A (2021) - Applied Numerical Computing, Prof. L.
#' Vandenberghe, UCLA. https://www.seas.ucla.edu/~vandenbe/ee133a.html.
#' 7. Linear Equations; 13. Cholesky factorization
#'
#' S. Boyd and L. Vandenberghe (2018). Introduction to Applied Linear Algebra
#' – Vectors, Matrices, and Least Squares. p. 207 - 210.
#' https://web.stanford.edu/~boyd/vmls/
#'
#' https://www.rdocumentation.org/packages/Matrix/versions/1.5-4/topics/Cholesky
#'
#'
#' Arguments
#'
#' @A
#'
#' A sparse symmetric positive definite matrix
#'
#' @b
#'
#' The right-hand side (RHS) of the Mixed Model Equations (MME) or b = NULL.
#'
#' @condNumber
#'
#' Logical indicating if the condition number for matrix invertion will be
#' returned. The default is TRUE.
#'
#' @tol
#'
#' The tolerance for detecting linear dependencies in the columns of A.
#' The default is .Machine[['double.eps']]
#'
#' @Usage:
#'
#' solve.CHMperm(A, b, condNumber = FALSE, tol = NULL)
#' return the vector of solutions for the system Ax = b
#'
#' solve.CHMperm(A, condNumber = FALSE, tol = NULL)
#' return the inverse of the A matrix
#'
#' @example
#'
#' Example 12.1 (pag 206) from "Linear Models for the Prediction of Animal
#' Breeding Values" by Raphael A. Mrode, 3rd Edition (2014).
#'
#' library(data.table)
#' mum <- c(0,0,0,0,2,4,5,5,8,8,8,8)# mothers
#' dad <- c(0,0,0,0,1,3,6,0,3,3,6,6) # fathers
#' ID <- c(seq(1:length(dad))) # individuals
#' pedigree <- data.table(ID, mum, dad)
#' F <- makeF(pedigree)[["F"]]
#' Finv <- solve.CHMperm(F) # inverse of the F matrix

solve.CHMperm <- function(C,
                          b = NULL,
                          cond.number = FALSE,
                          tol = NULL) {
  library(Matrix, quietly = TRUE)
  library(data.table, quietly = TRUE)

  source("condNumber.r")

  if (inherits(C, "matrix")) {
    C <- Matrix(C, sparse = TRUE)
  } else if (!inherits(C, "dgCMatrix") ||
    inherits(C, "dgTMatrix") ||
    inherits(C, "dgRMatrix")) {
    C <- as(as(C, "sparseMatrix"), "dgCMatrix")
  }

  # Cholesky factorization works only for symmetric matrices!
  if (!isSymmetric(C)) {
    C <- forceSymmetric(new(
      "dgCMatrix",
      i = C@i,
      p = C@p,
      Dim = C@Dim,
      x = C@x,
      Dimnames = list(
        C@Dimnames[[1]],
        C@Dimnames[[2]]
      )
    ))
    C <- drop0(C,
      tol = 1e-15,
      is.Csparse = NA
    )
  }

  if (!is.null(tol)) {
    if ((!is.numeric(tol) || length(tol) != 1 || tol < 0)) {
      stop("Argument 'tol' must be a non-negative numeric scalar")
    }
  } else {
    tol <- .Machine[["double.eps"]]
  }

  # compute and check the condition number (cnd) for matrix invertion
  if (cond.number) {
    cnd <- condNumber(C, norm = "2")
    if (cnd > 4.5e15) {
      warning("The Matrix is ill-conditioned")
    }
  }

  # build a sparse supernodal (LL' = P'AP) Cholesky factorization of C matrix
  # (Chen et al. 2008)
  CHMp <- Cholesky(C,
    perm = TRUE,
    super = TRUE
  )
  L <- as(CHMp, "Matrix")
  P <- as(CHMp, "pMatrix")
  if (is.null(b)) {
    Ic <- .symDiagonal(dim(C)[1L])
    c <- solve(CHMp, Ic, system = "P", tol = tol)
  } else {
    c <- solve(CHMp, b, system = "P", tol = tol)
  }
  y <- solve(L, c, system = "L", tol = tol)
  z <- solve(t(L), y, system = "Lt", tol = tol)
  if (!is.null(b)) { # get solutions for the system Ax = b
    z0 <- as.vector(solve(P, z, system = "P", tol = tol))
    z0[]
  } else { # return the inverse of the A matrix
    z0 <- as(
      solve(P, z, system = "P", tol = tol),
      "symmetricMatrix"
    )
    z0@Dimnames[[1L]] <- z0@Dimnames[[2L]] <-
      C@Dimnames[[1L]]
    z0[]
  }
}


#' loglik() function
#'
#' function to compute the REML log-likelihood for strategies I, II or III
#' under the gamma parametrization
#'
#' References:
#'
#' Gilmour, A., Thompson, R. and Cullis, B. (1995). Average information REML:
#' An efficient algorithm for variance parameter estimation in linear mixed
#' models. Biometrics 51 1440–1450.
#'
#' Harville, D. (1977). Maximum likelihood approaches to variance component
#' estimation and to related problems. Journal of the American Statistical
#' Association 72 320–340.

loglik <- function(pedigree,
                   h2, # narrow-sense heritability
                   H2, # broad-sense heritability
                   X, # incidence matrix for fixed effects
                   y, # a vector of phenotypes
                   Z, # incidence matrix for additive/dominance/mendel effects
                   F = NULL, # relatedness matrix between full-sib families
                   D = NULL, # dominance relatedness matrix
                   U = NULL, # incidence matrix for sca effects
                   eig.tol = 1e-6, # tolerance for positive of eigenvalues
                   C) { # the coefficient matrix for the mixed model equations

  source("makeKinship.nonadj.r")
  source("makeRelatedness.r")

  C <- forceSymmetric(new(
    "dgCMatrix",
    i = C@i,
    p = C@p,
    Dim = C@Dim,
    x = C@x,
    Dimnames = list(
      C@Dimnames[[1]],
      C@Dimnames[[2]]
    )
  ))
  C <- drop0(C,
    tol = 1e-15,
    is.Csparse = NA
  )

  if (!is.null(F)) {
    q <- dim(F)[2]
  } else {
    q <- 0
  }
  p <- dim(X)[2]
  m <- dim(pedigree)[1]
  n <- dim(Z)[1]

  # check if pedigree has the inbreeding coefficient ("f")
  # if not, create it
  if (("f" %in% names(pedigree))) {
    f <- as.vector(pedigree[, f])
  } else {
    pedigree[, f := calcInbreeding(pedigree[, c(1L:3L)])]
    f <- as.vector(pedigree[, f])
  }

  # create the additive reledtness matrix
  A <- makeKinship.nonadj(
    pedigree = pedigree,
    familyPed = FALSE
  )[["A"]]

  # setting sigmas
  sigma2.add <- h2
  sigma2.dom <- H2 - h2
  sigma2.sca <- sigma2.dom / 4
  sigma2.e <- 1 - H2

  # setting gammas
  gamma.add <- sigma2.add / sigma2.e
  gamma.sca <- sigma2.sca / sigma2.e
  gamma.dom <- sigma2.dom / sigma2.e

  if (!is.null(D) & !is.null(F) & (dim(C)[1] == 2 * m + p + q |
    dim(C)[1] == 2 * m + p + q + 1)) { # strategy S3
    # y ~ Xb + Za + Zd + e 
    # C is the coefficient matrix for additive (a) and dominance (d) effects
    # : C is of order (2m + p + q)x(2m + p + q) # augmented by the MMEs from sca
    cat("Strategy S3 ...\n")

    # retrieve the MMEs for mu & additive effects
    if (dim(C)[1] == 2 * m + p + q) {
      C0 <- C[1:(p + m), 1:(p + 2 * m)] # (p + m)x(p + m) matrix
      # retrieve the MMEs for mu, b & dditive effects
    } else {
      C0 <- C[1:(p + 1 + m), 1:(p + 1 + 2 * m)] # (p + 1 + m)x(p + 1 + 2m) matrix
    }
    # retrieve the MMEs for dominance effects
    # create "directly" the inverse of the dominance matrix (D)
    Dinv <- solve.CHMperm(D)
    alpha2 <- sigma2.e / sigma2.dom
    if (dim(C)[1] == 2 * m + p + q) { # y ~ 1mu + Za + Zd + e
      C01 <- cbind(
        crossprod(Z, X), # mx1
        crossprod(Z), # mxm
        crossprod(Z) + alpha2 * Dinv # mxm
      )
      W <- cbind(X, Z, Z)
    } else { # inbreeding coefficient (f) is included as a covariable, i.e.,
      # y ~ 1mu + Zfb + Za + Zd + e
      Zf <- Z %*% f
      C01 <- cbind( # mx(p + 1 + 2m)
        crossprod(Z, X), # mx1
        crossprod(Z, Zf), # mx1
        crossprod(Z), # mxm
        crossprod(Z) + alpha2 * Dinv # mxm
      )
      W <- cbind(X, Zf, Z, Z)
    }
    # recovering the original coefficient matrix (C)
    C <- rbind(C0, C01)
    strategy <- "S3"
  } else if (!is.null(D) & is.null(F) & (dim(C)[1] == 2 * m + p |
    dim(C)[1] == 2 * m + p + 1)) { # strategy S1

    if (dim(C)[1] == 2 * m + p) { # y ~ Xb + Za + Zd + e
      W <- cbind(X, Z, Z)
    } else { # y ~ 1mu + Zfb + Za + Zd + e
      Zf <- Z %*% f
      W <- cbind(X, Zf, Z, Z)
    }
    cat("Strategy S1 ...\n")
    strategy <- "S1"
  } else if (is.null(D) & !is.null(F) & !is.null(U) &
    (dim(C)[1] == p + 2 * m + q |
      dim(C)[1] == p + 2 * m + q + 1)) { # strategy S2
    if (dim(C)[1] == p + 2 * m + q) { # y ~ 1mu + Za + Usca + Zm + e
      W <- cbind(X, Z, U, Z)
    } else { # y ~ 1mu + Zfb + Za + Usca + Zm + e
      Zf <- Z %*% f
      W <- cbind(X, Zf, Z, U, Z)
    }
    cat("Strategy S2 ...\n")
    strategy <- "S2"
  }

  if (strategy == "S1" | strategy == "S3") {
    # create G-submatrices
    G.add.gamm <- gamma.add * A
    G.dom.gamm <- gamma.dom * D
    # create a null matrix
    zzeros2 <- zeros.matrix(m, m)
    # make the G matrix
    G.gamm <- rbind(
      cbind(G.add.gamm, zzeros2),
      cbind(t(zzeros2), G.dom.gamm)
    )
  } else if (strategy == "S2") {
    # cat("Strategy S2 ...\n")
    strategy <- "S2"
    # create the relatedness (M) matrix
    M <- Mend.adj(pedigree)

    # create G-submatrices
    G.add.gamm <- gamma.add * A
    G.sca.gamm <- gamma.sca * F
    G.mend.gamm <- gamma.dom * M

    # create null matrices to build the G matrix
    q <- dim(F)[2]
    zzeros1 <- zeros.matrix(m, q)
    zzeros3 <- zeros.matrix(q, m)
    zzeros2 <- zeros.matrix(m, m)

    # make the G matrix
    G.gamm <- rbind(
      cbind(G.add.gamm, zzeros1, zzeros2),
      cbind(t(zzeros1), G.sca.gamm, zzeros3),
      cbind(t(zzeros2), t(zzeros3), G.mend.gamm)
    )
  }

  # create the inverse of the coefficient matrix C
  Cinv <- solve.CHMperm(C)

  # the covariance matrix, R (error terms) it assumed as an identity matrix
  # i.e., R = sigma2.e / sigma2.e * I = I
  Rinv <- R.gamm <- Diagonal(n)

  # An alternative (but more efficient) expression to compute the matrix P
  # is given by (Gilmour et al., 1995, p. 1441)
  P <- Rinv - Rinv %*% W %*% Cinv %*% t(W) %*% Rinv

  # the REML log-likelihood (Harville 1977)
  # loglik = -0.5(log|G| + log|R| + log|C| + ypPy)

  # Cholesky decomposition is used to compute the log determinant
  # more efficiently and accurately...!

  # Note that
  # 0.5*log(det(G)) = sum(log(diag(chol.G)))
  # 0.5*log(det(R)) = sum(log(diag(chol.R)))
  # 0.5*log(det(C)) = sum(log(diag(chol.C)))

  # Cholesky decomposition of G, R and C
  chol.G <- chol(G.gamm)
  chol.R <- chol(R.gamm)
  chol.C <- chol(C)

  # check if "P" matrix is positive definite, if not use check.PD() function to
  # approximate it to a nearest positive definite matrix (Higham 2002)
  P <- check.PD(P, eig.tol = 1e-6)

  # Cholesky decomposition can also be used to compute yp*P*y, i.e.,
  # yp*P*y = yp*(Rp*R)*y = zp*z
  # where yp = t(y); z = Rp*y; Rp = t(R); R = chol(P); zp = t(z)
  chol.P <- chol(P)
  z <- chol.P %*% y
  ypPy <- as.vector(crossprod(z))

  loglik <- -sum(log(diag(chol.G))) -
    sum(log(diag(chol.R))) -
    sum(log(diag(chol.C))) -
    0.5 * ypPy

  loglik[]
}


#' check.PD() function
#'
#' function to check if a matrix is positive definite
#' if not, approximate it to a nearest positive definite matrix
#'
#' References:
#'
#' Higham, Nick (2002) Computing the nearest correlation matrix - a problem
#' from finance; IMA Journal of Numerical Analysis 22, 329–343.

check.PD <- function(C, # the matrix to be checked
                     eig.tol = NULL) { # specify the tolerance to define which
  # eigenvalues are "numerically" zero
  # compute the eigenvalues of the C matrix
  eigen <- eigen(C, symmetric = TRUE)[["values"]]
  if (!is.null(eig.tol)) {
    if ((!is.numeric(eig.tol) || length(eig.tol) != 1 || eig.tol < 0)) {
      stop("Argument 'eig.tol' must be a non-negative numeric scalar")
    }
  } else { # compute eig.tol from eigenvalues of the C matrix
    eig.tol <- length(eigen) * eps(eigen) # for details, see eps() function
  }
  # checking if some eigenvalue is less than eig.tol
  PD <- ifelse(any(eigen < eig.tol), "TRUE", "FALSE")
  if (PD) { # if so, compute the nearest positive definite matrix (Higham 2002)
    C.pd <- as(
      nearPD(C,
        corr = FALSE, # C is a covariance matrix
        keepDiag = FALSE, # do not keep the original diagonal values
        eig.tol = eig.tol,
        conv.tol = 1e-5, # convergence tolerance for Higham's (2002) algorithm
        maxit = 300 # default maximum number of iterations
      )[['mat']],
      "sparseMatrix"
    )
  } else {
    C.pd <- C
  }
  C.pd[]
}


# function that returns the machine epsilon of x object
# x can be a numeric vector, matrix or array

eps <- function(x = 1.0) {
  stopifnot(is.numeric(x))
  x <- max(abs(x))
  if (x < .Machine[['double.xmin']]) {  
    ep <- .Machine[['double.xmin']]
  } else {
    ep <- 2^floor(log2(x)) * .Machine[['double.eps']]
  }
  ep[]
}

# function to create the M relatedness or its inverse (ginv = TRUE)
# adjusted (adj = TRUE) or not adjusted by inbreeding (adj = FALSE)

makeM <- function(pedigree,
                  ginv = FALSE, # ginv = TRUE returns the inverse of M
                  adj = FALSE) { # adj = TRUE returns the M matrix adjusted
  # by inbreeding

  if (!is.data.table(pedigree)) {
    if (is.data.frame(pedigree) || is.matrix(pedigree)) {
      pedigree <- as.data.table(pedigree)
    } else {
      stop("nothing to do...")
    }
  }
  setnames(pedigree,
    old = c(1:3),
    new = c("TreeID", "mum", "dad")
  )
  # m0 is the number of phantom parents
  m0 <- length(pedigree[mum == 0L | dad == 0L, TreeID])
  # s is the number of progenies
  s <- dim(pedigree)[1L] - m0
  M <- .symDiagonal(m0 + s)
  M@Dimnames[[1L]] <- M@Dimnames[[2L]] <-
    as.character(pedigree[[1L]])

  if (adj) { # diagonal elements of M are adjusted by inbreeding
    if (!("f" %in% names(pedigree))) {
      f <- as.vector(calcInbreeding(pedigree[, c(1L:3L)]))
    } else {
      f <- as.vector(pedigree[, f])
    }
    if (!ginv) { # return the M matrix adjusted by inbreeding
      M[cbind(
        seq(m0 + 1L, m0 + s),
        seq(m0 + 1L, m0 + s)
      )] <- 3 / 4 * (1 - f[seq(m0 + 1L, m0 + s)])
      return(M)
    } else { # return the inverse of the M matrix adjusted by inbreeding
      M[cbind(
        seq(m0 + 1L, m0 + s),
        seq(m0 + 1L, m0 + s)
      )] <- 1 / (3 / 4 * (1 - f[seq(m0 + 1L, m0 + s)]))
      return(M)
    }
  } else { # M is not adjusted by inbreeding
    if (!ginv) { # return the M matrix
      M[cbind(
        seq(m0 + 1L, m0 + s),
        seq(m0 + 1L, m0 + s)
      )] <- 3 / 4
      return(M)
    } else { # return the inverse of the M matrix
      M[cbind(
        seq(m0 + 1L, m0 + s),
        seq(m0 + 1L, m0 + s)
      )] <- 4 / 3
      return(M)
    }
  }
}


# function to create the family relatedness (F) matrix assuming that
# families are unrelated (ide)

makeF.ide <- function(pedigree) {
  library(data.table)
  if (!is.data.table(pedigree)) {
    if (is.data.frame(pedigree) || is.matrix(pedigree)) {
      pedigree <- as.data.table(pedigree)
    } else {
      stop("nothing to do...")
    }
  }
  ped3 <- pedigree[, c(1L:3L)]
  setnames(
      ped3,
      c("TreeID", "mum", "dad")
  )
  cross <- as.vector(unique(makeFam(ped3)[, cross]))
  # label rows/cols of the F matrix
  cross.label <- cross[!is.na(cross)]
  # set the dimension of the F (ide) matrix
  q <- length(cross.label)
  F.ide <- .symDiagonal(q)
  F.ide@Dimnames[[1L]] <- F.ide@Dimnames[[2L]] <-
    cross.label
  F.ide[]
}


makeF.ide2 <- function(pedigree) {
  library(data.table)
  if (!is.data.table(pedigree)) {
    if (is.data.frame(pedigree) || is.matrix(pedigree)) {
      pedigree <- as.data.table(pedigree)
    } else {
      stop("nothing to do...")
    }
  }
  ped3 <- pedigree[, c(1L:3L)]
  setnames(
      ped3,
      c("TreeID", "mum", "dad")
  )
  # m0 is the number of phantom parents
  m0 <- length(pedigree[mum == 0L | dad == 0L, TreeID])
  # s is the number of progenies
  s <- dim(pedigree)[1L] - m0
  #cross <- as.vector(unique(makeFam(ped3)[, cross]))
  # label rows/cols of the F matrix
  #cross.label <- cross[!is.na(cross)]
  # set the dimension of the F (ide) matrix
  #q <- length(cross.label)
  F.ide <- .symDiagonal(s)
  #F.ide@Dimnames[[1L]] <- F.ide@Dimnames[[2L]] <-
   # cross.label
  F.ide[]
}


# function to create the dominance relatedness matrix (D.ide) assuming that
# full-sib families are unrelated

makeD.ide <- function(pedigree) {
  library(data.table)
  if (!is.data.table(pedigree)) {
    if (is.data.frame(pedigree) || is.matrix(pedigree)) {
      pedigree <- as.data.table(pedigree)
    } else {
      stop("nothing to do...")
    }
  }
  m <- dim(pedigree)[1]
  D.ide <- .symDiagonal(m)
  D.ide@Dimnames[[1L]] <- D.ide@Dimnames[[2L]] <-
    as.character(pedigree[, TreeID])
  D.ide[]
}


# function to create the Z.dom matrix assuming that
# full-sib families are unrelated

Zdom.ide <- function(pedigree) {
  if (!is.data.table(pedigree)) {
    if (is.data.frame(pedigree) || is.matrix(pedigree)) {
      pedigree <- as.data.table(pedigree)
    } else {
      stop("nothing to do...")
    }
  }
  ped3 <- pedigree[, c(1L:3L)]
  setnames(
      ped3,
      c("TreeID", "mum", "dad")
  )
  cross <- as.vector(unique(makeFam(ped3)[, cross]))
  # label rows/cols of the F matrix
  cross.label <- cross[!is.na(cross)]
  # set the dimension of the F (ide) matrix
  q <- length(cross.label)
  # add the gen column to the ped3 data.table
  ped3[, gen := gen.add(ped3)[, gen]]
  n <- dim(ped3[gen > 0L, ])[1]
  m <- dim(pedigree)[1L]
  zeros_base <- zeros.matrix(
    n.rows = m - q,
    n.cols = q
  )
  Iq <- .symDiagonal(q)
  Zdom_ide <- rbind(zeros_base, Iq)
  Zdom_ide@Dimnames[[1L]] <-
    as.character(ped3[, TreeID])
  Zdom_ide@Dimnames[[2L]] <- cross.label
  return(Zdom_ide)
}

Zdom.ide2 <- function(pedigree) {
  if (!is.data.table(pedigree)) {
    if (is.data.frame(pedigree) || is.matrix(pedigree)) {
      pedigree <- as.data.table(pedigree)
    } else {
      stop("nothing to do...")
    }
  }
  
  setnames(
      pedigree, c(1:3),
      c("TreeID", "mum", "dad")
  )
  
  # m is the number of individuals in the pedigree
  m <- dim(pedigree)[1L]
  # m0 is the number of phantom parents
  m0 <- length(pedigree[mum == 0L | dad == 0L, TreeID])
  # s is the number of progenies 
  s <- m - m0
  Iq <- .symDiagonal(s)

  zeros_base <- zeros.matrix(
    n.rows = m - s,
    n.cols = s
  )
  Zdom_ide <- rbind(zeros_base, Iq)
  #Zdom_ide@Dimnames[[1L]] <-
   # as.character(pedigree[, TreeID])
  #Zdom_ide@Dimnames[[2L]] <- cross.label
  return(Zdom_ide)
}


# function to build the incidence matrix (Z.dom) relating an individual
# in the pedigree to its family (if it proceeds)

Z.dom <- function(pedigree) {
  if (!is.data.table(pedigree)) {
    if (is.data.frame(pedigree) || is.matrix(pedigree)) {
      pedigree <- as.data.table(pedigree)
    } else {
      stop("nothing to do...")
    }
  }
  datafam <- pedigree[, c(1L:3L)]
  setnames(
      datafam,
      c("TreeID", "mum", "dad")
  )
  datafam[, cross := makeFam(datafam)[, cross]]
  form <- formula(~ cross - 1)
  termsf <- terms(form,
    keep.order = TRUE
  )
  mf <- model.frame(
    termsf,
    data = datafam,
    na.action = na.pass
  )
  Zdom <- MatrixModels::model.Matrix(
    form, mf,
    sparse = TRUE
  )
  Zdom <- as(Zdom, "sparseMatrix")
  Zdom@Dimnames[[1L]] <- c(paste0(
    "TreeID",
    pedigree[[1L]]
  ))
  fam <- as.vector(unique(datafam[, cross]))
  Zdom@Dimnames[[2L]] <- fam[!is.na(fam)]
  return(list(
    Zdom = Zdom,
    datafam = datafam
  ))
}


# function to build the incidence matrix (Z.fam) for full-sib families

Z.fam <- function(trial) {
  if (!is.data.table(trial)) {
    if (is.data.frame(trial) || is.matrix(trial)) {
      trial <- as.data.table(trial)
    } else {
      stop("nothing to do...")
    }
  }
  trial.tmp <- copy(trial)[, c(1L:3L)]
  setnames(
      trial.tmp,
      c("TreeID", "mum", "dad")
  )
  # check if trial has a column named "cross"
  # if not, create it
  if (!("cross" %in% names(trial.tmp))) {
    trial.tmp[, cross := makeFam(trial.tmp)[, cross]]
  }
  # create the incidence matrix for full-sib families
  form <- formula(~ cross - 1)
  termsf <- terms(form,
    keep.order = TRUE
  )
  mf <- model.frame(
    termsf,
    data = trial.tmp,
    na.action = na.pass
  )
  Zfam <- MatrixModels::model.Matrix(
    form, mf,
    sparse = TRUE
  )
  Zfam <- as(Zfam, "dgCMatrix")
  Zfam@Dimnames[[1L]] <- c(paste0(
    "TreeID",
    trial.tmp[[1L]]
  ))
  fam <- as.vector(unique(trial.tmp[, cross]))
  Zfam@Dimnames[[2L]] <- fam[!is.na(fam)]
  Zfam[]
}


# function to build the design matrix (Z) under the assumption that each
# individual has one record

Z.mat <- function(pedigree) {
  if (!is.data.table(pedigree)) {
    if (is.data.frame(pedigree) || is.matrix(pedigree)) {
      pedigree <- as.data.table(pedigree)
    } else {
      stop("nothing to do...")
    }
  }
  if (!any(names(pedigree) == "gen")) {
    pedigree <- gen.add(pedigree)
  }
  ped.trial <- pedigree[gen > 0L, c(1L:3L)]
  setnames(
      ped.trial,
      c("TreeID", "mum", "dad")
  )
  ped.trial[, TreeID := factor(TreeID,
    levels = unique(TreeID)
  )]
  form <- formula(~ TreeID - 1)
  termsf <- terms(form, keep.order = TRUE)
  mf <- model.frame(
    termsf,
    data = ped.trial,
    na.action = na.pass
  )
  Zdata <- MatrixModels::model.Matrix(
    form, mf,
    sparse = TRUE
  )
  m <- dim(pedigree)[1L]
  n <- dim(ped.trial)[1L]
  Zbase <- zeros.matrix(
    n.rows = n,
    n.cols = m - n
  )
  Z <- cbind(Zbase, Zdata)
  Z@Dimnames[[2L]] <-
    c(paste0(
      "TreeID",
      pedigree[[1L]]
    ))
  Z@Dimnames[[1L]] <-
    c(paste0(
      "TreeID",
      ped.trial[[1L]]
    ))
  return(Z)
}


# function to built the inverse of the numerator relationship matrix
# based on Quaas (1976)

# Reference:
# Quaas, R.L. (1976) Computing the diagonal elements of a large numerator
# relationship matrix. Biometrics 32, 949–953.

Ainverse <- function(pedigree) {
  if (!is.data.table(pedigree)) {
    if (is.data.frame(pedigree) || is.matrix(pedigree)) {
      pedigree <- as.data.table(pedigree)
    } else {
      stop("nothing to do...")
    }
  }
  pedigree::makeAinv(pedigree[, c(1:3)])
  output <- fread("Ainv.txt")
  setnames(output,
    old = colnames(output),
    new = c("x", "y", "ai")
  )
  m <- dim(pedigree)[1L]

  # build vectors index of unique values of x and y
  ux <- as.vector(unique(output[, x]))
  idx <- base::order(ux)
  index.x <- data.table(
    x = ux,
    idx = idx
  )
  output <- index.x[output,
    on = .(x = x)
  ]
  uy <- as.vector(unique(output[, y]))
  idy <- base::order(uy)
  index.y <- data.table(
    y = uy,
    idy = idy
  )

  output <- index.y[output,
    on = .(y = y)
  ]
  cols <- c("idx", "idy")
  output[, (cols) := lapply(.SD, as.integer),
    .SDcols = cols
  ]

  Ainv <- (with(
    output,
    Matrix::sparseMatrix(
      i = idx,
      j = idy,
      x = ai,
      dims = c(m, m),
      dimnames = list(
        pedigree[[1L]],
        pedigree[[1L]]
      ),
      triangular = FALSE,
      check = TRUE
    )
  ))

  Ainv.s <- forceSymmetric(Ainv,
    uplo = "L"
  )
  Ainv.s <- drop0(Ainv.s,
    tol = 1e-15,
    is.Csparse = NA
  )
  return(Ainv.s)
}


# function to make a matrix of ones in the compressed sparse column format
# ("dgCMatrix")

ones.matrix <- function(n.rows, n.cols) {
  m <- Matrix::Matrix(
    nrow = n.rows,
    ncol = n.cols,
    data = 1L,
    sparse = TRUE
  )
  return(m)
}


# function to make a matrix of zeros in the compressed sparse column format
# ("dgCMatrix")

zeros.matrix <- function(n.rows, n.cols) {
  m <- Matrix::Matrix(
    nrow = n.rows,
    ncol = n.cols,
    data = 0L,
    sparse = TRUE
  )
  return(m)
}


# function to add a "generation" field to the pedigree

gen.add <- function(pedigree) {
  if (!is.data.table(pedigree)) {
    if (is.data.frame(pedigree) || is.matrix(pedigree)) {
      pedigree <- as.data.table(pedigree)
    } else {
      stop("nothing to do...")
    }
  }
  ped3 <- pedigree[, c(1L:3L)]
  setnames(
      ped3,
      c("TreeID", "mum", "dad")
  )
  ped3[, gen := 0L]
  for (i in 0:nrow(ped3)) {
    ped3[
      mum %in% TreeID[gen == i] | dad %in% TreeID[gen == i],
      gen := gen + 1L
    ]
    if (.Last.updated == 0L) break
  }
  ped3[]
}


# function to add a "cross" field to trial/pedigree

makeFam <- function(pedigree) {
  ped.cross <- copy(pedigree)
  setnames(
    ped.cross,
    c(1:3), c("TreeID", "mum", "dad")
  )
  ped.cross[
    ,
    cross := ifelse((mum == 0L | dad == 0L),
      NA,
      paste(pmin(mum, dad),
        pmax(mum, dad),
        sep = "x"
      )
    )
  ]
  ped.cross[, cross := factor(cross,
    levels = unique(cross)
  )]
  ped.cross[]
}
