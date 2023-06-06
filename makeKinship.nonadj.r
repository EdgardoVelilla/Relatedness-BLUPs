#' makeKinship function
#'
#' Author: Edgardo Velilla P.
#' email{edgardo.velilla@cmpc.cl}
#' Created: 12-Mar-2021
#' Modified: 05-May-2022
#' License: GPLv3
#' References:
#'
#' Cockerham, C.C.(1954). An extension of the concept of partitioning
#' hereditary variance for analysis of covariances among relatives when
#' epistasis is present. Genetics 39, 859–882.
#'
#' Henderson, C. R. (1985). Best linear unbiased prediction of nonadditive
#' genetic merits. Journal of Animal Science, 60, 111–117.
#'
#' Hoeschele, I. and VanRaden, P.M. (1991). Rapid inverse of dominance
#' relationship matrices for noninbred populations by including sire by dam
#' subclass effects. Journal of Dairy Science 74, 557–569.
#'
#' General description:
#' Generates several relatedness matrices in sparse and/or triplet format from
#' a pedigree frame. The critical assumption is that the algorithm for
#' computing the relatedness matrix between full-sib families (or dominance
#' effects) is valid only for a non-inbred populations.
#'
#' Arguments
#'
#' @param  pedigree
#'
#' A data.table/dataframe/matrix where the first three columns correspond to
#' the identifiers for the individual, mother and father, respectively. The row
#' giving the pedigree of an individual must appears before any row where that
#' individual appears as parent. Founders use 0 (zero) in the parental columns.
#'
#' @param init.ghost
#'
#' Integer number indicating the initial ID code for ghost fathers to match
#' pedigree for selection from 1st open-pollination population. The default is
#' @init.ghost = NULL.
#'
#' @par familyPed
#'
#' Logical indicating if a family pedigree (but keeping all individuals across
#' the pedigree, which represents the critical link for each family) should be
#' created to build the relateddness matrix, F, between full-sib failies. The
#' default is @familyPed = FALSE.
#'
#' @param Finv
#'
#' Logical indicating if the inverse of the dominance relatedness matrix among
#' subclass effects (F), -as described in Hoeschele & VanRaden (1991)- should
#' be returned. This matrix is returned in sparse triplet format according to
#' ASReml-R requirements. The default is @Finv = FALSE.
#'
#' @param fi
#'
#' Logical indicating if the individual inbreeding coefficient (fi) should be
#' returned. The default is @fi = FALSE.
#'
#' @param NRM
#'
#' Logical indicating if the numerator relationship matrix (A) should be
#' returned. The default is @NRM = FALSE.
#'
#' @param Dinv
#'
#' Logical indicating if the inverse of dominance relatedness matrix (D) should
#' be returned. This matrix is returned in sparse triplet format according to
#' ASReml-R requirements. The default is the @Dinv = FALSE.
#'
#'
#' @return
#'
#' A potential list with the following components:
#'
#' The numerator relationship matrix (A) in sparse matrix format.
#'
#' The relatedness matrix among parent subclass effects (F) in sparse matrix
#' format. This matrix is equivalent to covariance between full-sib families
#' due to dominance as described in Hoeschele & VanRaden (1991).
#'
#' The inverse of the F matrix (Finv.as) in three column sparse coordinate form
#' matrix in row major order, with attribute "INVERSE" set to TRUE and
#' attribute rowNames.
#'
#' The relatedness matrix corresponding to Mendelian Sampling Term (Dw) for
#' the non-additive effect in sparse matrix format.
#'
#' The inverse of the Dw matrix in sparse matrix format (Dwinv).
#'
#' The inverse of the Dw matrix in three column sparse coordinate form matrix
#' in row major order (Dwinv.as), with attribute "INVERSE" set to TRUE and
#' attribute rowNames.
#'
#' The dominance relatedness matrix (D) between individuals in sparse matrix
#' format.
#'
#' The inverse of dominance relatedness matrix (Dinv) in sparse matrix format.
#'
#' The inverse of the dominance relatedness matrix in three column sparse
#' coordinate form matrix in row major order (Dinv.as), with attribute
#' "INVERSE" set to TRUE and attribute rowNames.
#'
#' A numeric vector (f) containing the inbreeding coefficient for each
#' individual.
#'
#' @example
#' Example 12.1 (pag 206) from "Linear Models for the Prediction of Animal
#' Breeding Values" by Raphael A. Mrode, 3rd Edition (2014)
#'
#' # pedigre information
#' sire <- c(0,0,0,0,1,3,6,0,3,3,6,6)
#' dam <- c(0,0,0,0,2,4,5,5,8,8,8,8)
#' ID <- c(seq(1:12)) # individuals
#' pedm9 <- data.table(ID=ID, mum=dam, dad=sire)
#' matm9 <- makeKinship.nonadj(pedm9, fi=TRUE, NRM=TRUE, Dinv=TRUE)
#' > names(matm9)
#' [1] "A"  "F" "Finv" "Finv.as"  "Dw"  "Dwinv" "Dwinv.as"  "D"  "Dinv"
#' "Dinv.as"  "f"

makeKinship.nonadj <- function(pedigree,
                               familyPed = FALSE,
                               Finv = FALSE,
                               fi = FALSE,
                               NRM = FALSE,
                               Dinv = FALSE,
                               tol = NULL) {

  library(data.table, quietly = TRUE)
  requireNamespace("MatrixModels", quietly = TRUE)
  library(Matrix, quietly = TRUE)

  source("forASReml.r")
  source("matrixMethods.r")
  source("pedFam.r")
  source("condNumber.r")

  if (is.data.frame(pedigree)) {
    pedigree <- as.data.table(pedigree)
  } else if (is.matrix(pedigree)) {
    pedigree <- as.data.table(pedigree)
  } else if (!is.data.table(pedigree)) {
    stop("nothing to do...")
  }

  pedi <- copy(pedigree)
  pedi <- pedi[, c(1:3)]
  setnames(pedi, c("TreeID", "mum", "dad"))

  if (familyPed) {
    ped <- pedFam(pedi)
  } else {
    ped <- pedi
  }
  TreeID <- ped[, TreeID]
  mum <- match(ped[, mum],
    ped[, TreeID],
    nomatch = NA
  )
  mum[is.na(mum)] <- 0L
  dad <- match(ped[, dad],
    ped[, TreeID],
    nomatch = NA
  )
  dad[is.na(dad)] <- 0L
  m <- dim(ped)[1]
  TreeID2 <- 1:m
  # family = mum x dad=dad x mum (no maternal or reciprocal effects...!)
  parents <- ped[, .(mum, dad)]
  parents[, family := ifelse((mum == 0L | dad == 0L), NA,
    paste(pmin(mum, dad), pmax(mum, dad), sep = "x")
  )]
  parents[, family := factor(family,
    levels = unique(family)
  )]
  datafam <- parents[, IID := TreeID2]
  IDfam <- datafam[!duplicated(family), ][!is.na(family), IID]
  if (familyPed) {
    # cat('starting to build F ...\n')
  } else {
    #  cat('starting to build D ...\n')
  }
  # cat("be patient ...\n")
  A <- Dbtemp <- Db <- f <- Dw <- D <- array(0, dim = c(m, m))
  rownames(A) <- rownames(D) <- rownames(Dbtemp) <- 
  rownames(Dw) <- ped[[1]]
  colnames(A) <- colnames(D) <- colnames(Dbtemp) <- 
  colnames(Dw) <- ped[[1]]
  rownames(Db) <- colnames(Db) <- TreeID2
  for (i in 1:m) {
    if (dad[i] != 0 && mum[i] != 0) { # Both parents are known
      for (j in 1:(i - 1)) {
        A[j, i] <- A[i, j] <- 0.5 * (A[j, dad[i]] + A[j, mum[i]])
        if (length(A[dad[i], dad[j]]) > 0) {
          Dss <- A[dad[i], dad[j]]
        } else {
          Dss <- 0L
        }
        if (length(A[mum[i], mum[j]]) > 0) {
          Ddd <- A[mum[i], mum[j]]
        } else {
          Ddd <- 0L
        }
        if (length(A[dad[i], mum[j]]) > 0) {
          Dsd <- A[dad[i], mum[j]]
        } else {
          Dsd <- 0L
        }
        if (length(A[mum[i], dad[j]]) > 0) {
          Dds <- A[mum[i], dad[j]]
        } else {
          Dds <- 0L
        }
        # (Cockerham 1954; Henderson 1985)
        Dbtemp[j, i] <- Dbtemp[i, j] <- D[i, j] <- D[j, i] <-
          0.25 * (Dss * Ddd + Dsd * Dds)
        # (Dss*Ddd + Dsd*Dds)=F
        # ==> var[f(S,D),f(K,L)]=(Dss*Ddd + Dsd*Dds)*var(f)=F*var(f)
        # (Hoeschele & VanRaden 1991)
        Db[j, i] <- Db[i, j] <- (Dss * Ddd + Dsd * Dds)
      }
      f[i, i] <- 0.5 * A[dad[i], mum[i]]
      A[TreeID2[i], TreeID2[i]] <- 1 + f[i, i]
      D[i, i] <- Db[i, i] <-
        (A[dad[i], dad[i]] * A[mum[i], mum[i]] + (A[dad[i], mum[i]])^2)
      Dw[i, i] <- (D[i, i] - 0.25 * Db[i, i])
    }
    if (dad[i] != 0 && mum[i] == 0) { # Only dad is known
      for (j in 1:(i - 1)) {
        A[j, i] <- A[i, j] <- 0.5 * (A[j, dad[i]])
        Db[j, i] <- Db[i, j] <- D[j, i] <- D[i, j] <- 0L
      }
      A[TreeID2[i], TreeID2[i]] <- Db[i, i] <- Dw[i, i] <- D[i, i] <- 1L
      f[i, i] <- 0L
    }
    if (dad[i] == 0 && mum[i] != 0) { # Only mum is known
      for (j in 1:(i - 1)) {
        A[j, i] <- A[i, j] <- 0.5 * (A[j, mum[i]])
        Db[j, i] <- Db[i, j] <- D[j, i] <- D[i, j] <- 0L
      }
      A[TreeID2[i], TreeID2[i]] <- Db[i, i] <- Dw[i, i] <- D[i, i] <- 1L
      f[i, i] <- 0L
    }
    if (dad[i] == 0 && mum[i] == 0) { # Both parents are unknown
      A[TreeID2[i], TreeID2[i]] <- Db[i, i] <- Dw[i, i] <- D[i, i] <- 1L
      f[i, i] <- 0L
    }
  }
  Dws <- as(Dw, "sparseMatrix")
  As <- as(A, "sparseMatrix")
  # the family relatedness matrix of order (qxq)
  F <- Db[IDfam, IDfam]
  Fs <- as(F, "sparseMatrix")
  if (!isSymmetric(Fs)) {
    Fs <- forceSymmetric(new(
      "dgCMatrix",
      i = Fs@i,
      p = Fs@p,
      Dim = Fs@Dim,
      x = Fs@x,
      Dimnames = list(
        Fs@Dimnames[[1]],
        Fs@Dimnames[[2]]
      )
    ))
    Fs <- drop0(Fs,
      tol = 1e-15,
      is.Csparse = NA
    )
  }
  fam <- as.vector(unique(datafam[, family]))
  Fs@Dimnames[[1]] <- Fs@Dimnames[[2]] <- fam[!is.na(fam)]
  q <- dim(Fs)[1]
  if (!fi) {
    f <- NULL
  } else {
    f <- as.vector(diag(f))
  }
  if (!NRM) A <- NULL
  if (familyPed) {
    M <- Mend.m(pedigree)

    ## compute the condition number for matrix inversion
    cnd <- condNumber(Fs, norm = "2")
    if (cnd > 4.5e15) {
      warning("Matrix 'F' is ill-conditioned")
    }

    ## defines which eigenvalues are numerically zero
    if (!is.null(tol)) {
      if ((!is.numeric(tol) || length(tol) != 1 || tol < 0)) {
        stop("Argument 'tol' must be a non-negative numeric scalar")
      }
    } else {
      eigen <- eigen(Fs, symmetric = TRUE)[["values"]]
      tol <- length(eigen) * eps(eigen)
    }
    ## make a masking using the positive relative eigenvalues
    eigen[abs(eigen) < tol] <- 0

    ## examing whether the matrix F is positive definite, negative definite,
    # and so on.
    if (!any(eigen > 0) & any(eigen == 0) & any(eigen < 0)) {
      stop("Matrix F seems negative semi-definite")
    } else if (any(eigen > 0) & any(eigen == 0) & !any(eigen < 0)) {
      warning("Matrix F seems positive semi-definite")
      status <- "positive semi-definite"
    } else if (!any(eigen < 0) & !any(eigen == 0)) {
      status <- "positive-definite"
    } else if (all(eigen < 0)) {
      status <- "negative definite"
    } else {
      status <- "indefinite"
    }
    if (any(abs(eigen) < tol)) {
      warning("Matrix 'F' is singular or nearly singular")
    }
    if (Finv) {
      F.inv <- solve.CHMperm(F,
        CheckPD = TRUE
      )
      Finv.as <- forASReml(
        F.inv,
        ginv = TRUE,
        colnames = c(
          "row",
          "column", "Finv"
        )
      )
      Dw.inv <- solve.CHMperm(M,
        CheckPD = TRUE
      )
      Dwinv.as <- forASReml(
        Dw.inv,
        ginv = TRUE,
        colnames = c(
          "row",
          "column", "Dwinv"
        )
      )
      cat("job done...! \n")
      return(list(
        A = As,
        F = Fs,
        Dw = M,
        Finv = F.inv,
        Finv.as = Finv.as,
        Dwinv = Dw.inv,
        Dwinv.as = Dwinv.as,
        f = f
      ))
    } else {
      cat("job done...! \n")
      return(list(
        A = As,
        F = Fs,
        cond = cnd,
        status = status,
        Dw = M,
        f = f
      ))
    }
  } else {
    Ds <- as(D, "sparseMatrix")
    if (Dinv) {
      cat("starting to invert D ...\n")
      cat("be patient ...\n")
      Dw.inv <- solve.CHMperm(Dws)
      Dwinv.as <- forASReml(
        Dw.inv,
        ginv = TRUE,
        colnames = c(
          "row",
          "column", "Dwinv"
        )
      )
      D.inv <- solve.CHMperm(Ds)
      Dinv.as <- forASReml(
        D.inv,
        ginv = TRUE,
        colnames = c(
          "row",
          "column", "Dinv"
        )
      )

      F.inv <- solve.CHMperm(
        Fs,
        CheckPD = TRUE
      )
      Finv.as <- forASReml(
        F.inv,
        ginv = TRUE,
        colnames = c(
          "row",
          "column", "Finv"
        )
      )
      return(list(
        A = As,
        F = Fs,
        Finv = F.inv,
        Finv.as = Finv.as,
        Dw = Dws,
        Dwinv = Dw.inv,
        Dwinv.as = Dwinv.as,
        D = Ds,
        Dinv = D.inv,
        Dinv.as = Dinv.as,
        f = f
      ))
    } else {
      return(list(
        A = As,
        F = Fs,
        Dw = Dws,
        D = Ds,
        f = f
      ))
    }
  }
}
