prediction.S1 <- function(pedigree,
                       ped.jacquard,
                       trial,
                       D = NULL,
                       Ainv = NULL,
                       init.ghost = NULL,
                       type = c("ide", "ks", "jq"),
                       ID = FALSE,
                       h2,
                       H2,
                       n.mod,
                       center = FALSE,
                       Loglik= FALSE) {

  library(Matrix, quietly = TRUE)
  library(data.table, quietly = TRUE)

  source("matrixMethods.r")
  source("makeKinship.nonadj.r")
  source("makeRelatedness.r")

  # check if pedigrees and trial objects are a data.table
  # if not, convert it.
  # Also, rename the first three columns as c("TreeID", "mum", "dad")
  convertAndRename <- function(object) {
      if (!is.data.table(object)) {
          if (is.data.frame(object) || is.matrix(object)) {
              object <- as.data.table(object)
          } else {
              stop("nothing to do...")
          }
      }
      setnames(object,
          old = c(1:3),
          new = c("TreeID", "mum", "dad")
      )
      object[, TreeID := as.integer(TreeID)]
      return(object)
  }

  pedigree <- convertAndRename(pedigree)
  ped.jacquard <- convertAndRename(ped.jacquard)
  trial <- convertAndRename(trial)

  # check if pedigree has a column named "gen"
  if (!("gen" %in% names(pedigree))) {
      pedigree[, gen := gen.add(pedigree)[, gen]] # if not, create it
      gen <- max(pedigree[, gen])
  } else {
      gen <- max(pedigree[, gen])
  }
  # check if trial has a column named "gen"
  trial.tmp <- copy(trial)
  if (!("gen" %in% names(trial.tmp))) {  # if not, recovery from pedigree
      trial.tmp[, gen := pedigree[, .(TreeID, gen)][trial.tmp[, .(TreeID)],
          on = .(TreeID = TreeID)
      ][, gen]]
  }

  # If the "cross" column does not exist, add it...!
  if (!("cross" %in% names(trial.tmp))) {
      trial.tmp[, cross := makeFam(trial.tmp)[, cross]]
  }

  # check if the trial contains the vector of phenotypic values ("pheno")
  if (!("pheno" %in% names(trial.tmp))) {
      stop("trial must have a field with phenotypic values")
  }
  # check if quantitative trait ("pheno") must be scaled
  if (center) { # if so, do it!
      trial.tmp[, pheno := scale(pheno,
          center = TRUE, scale = TRUE
      ),
      by = .(gen)
      ]
  }

  # create or check the inverse of A matrix
  if (is.null(Ainv)) {
      cat("starting to build Ainv ...\n")
      Ainv <- Ainverse(pedigree)
      cat("done...! \n")
  } else {
      if (!inherits(Ainv, "dsCMatrix")) {
          if (inherits(Ainv, "dgCMatrix") || inherits(
              Ainv, "dgTMatrix"
          ) || inherits(Ainv, "dgRMatrix")) {
              Ainv <- as(Ainv, "symmetricMatrix")
          } else if (inherits(Ainv, "matrix")) {
              Ainv <- as(
                  as(Ainv, "sparseMatrix"),
                  "symmetricMatrix"
              )
          } else {
              stop(
                  substitute(Ainv),
                  " must be a symmetric matrix from the Matrix class"
              )
          }
      }
  }

  # check if the D matrix is provided
  # if not, create it
  if (!is.null(D)) { # D matrix is not NULL
      if (inherits(D, "matrix")) { # check if D is a matrix
          if (!isSymmetric(D)) { # check if D is symmetric
              stop("matrix F must be symmetric!")
          } else { # D is a symmetric matrix
              D <- Matrix(D, sparse = TRUE) # create a sparse matrix
          }
      } else if (inherits(D, "dgCMatrix")) { # check if D is a dgCMatrix
          D <- as(D, "symmetricMatrix") # create a symmetric matrix
      } else if (!inherits(D, "dsCMatrix")) { # check if D is a dsCMatrix
          stop("D must be a symmetric matrix")
      }
  } else if (is.null(D) & type == "ide") { 
      # full-sib families are treated as unrelated
      # so the dominance relatedness matrix D is created as an identity matrix
      Dinv <- D <- D.ide(
          pedigree = pedigree
      )
      cat("ide ok! ...\n")
  } else if (is.null(D) & type == "ks") { # Cockerham (1954)
        # create the dominance matrix "D" from kinship coefficients
      D <- makeKinship.nonadj(
          pedigree = pedigree,
          familyPed = FALSE
      )[["D"]]
        # create the inverse of relatedness matrix (Dinv)
      Dinv <- solve.CHMperm(D)
      cat("ks ok! ...\n")
  } else if (is.null(D) & type == "jq") {
        # create the dominance matrix "D" from Jacquard's identity coefficients
      D <- makeRelatedness(
          pedigree = ped.jacquard
      )[["Dr"]]
      # create the inverse of D matrix (Dinv)
      Dinv <- solve.CHMperm(D)
      cat("Jacquard ok! ...\n")
  } else {
      stop(paste0("type=", type, " is not an option"))
  }

  if (ID) { # check if pedigree has the inbreeding coefficient ("f")
      if (!("f" %in% names(pedigree))) { # if not, create it
          pedigree[, f := calcInbreeding(pedigree[, c(1L:3L)])]
          f <- as.vector(pedigree[, f])
      } else {
          f <- as.vector(pedigree[, f])
      }
      if (!("f" %in% names(trial.tmp))) {
          trial.tmp[, f := pedigree[
              ,
              .(TreeID, f)
          ][trial.tmp[, .(TreeID)],
              on = .(TreeID = TreeID)
          ][, f]]
      }
  } else {
      f <- NULL
  }
  # ensure that the coefficient matrix (C) is non-singular
  if (!ID || (ID && all(f == 0))) { 
      model.base <- TRUE
  } else {
      model.base <- FALSE
  }

  # incidence matrices
  n <- dim(trial)[1]
  X <- ones.matrix(n, 1)
  Z <- Z.mat(pedigree)
  p <- dim(X)[2]
  m <- dim(pedigree)[1]
  y <- as.vector(trial.tmp[, pheno])

  # setting sigmas
  sigma2.add <- h2
  sigma2.dom <- H2 - h2
  sigma2.e <- 1 - H2

  alpha1 <- sigma2.e / sigma2.add
  alpha2 <- sigma2.e / sigma2.dom

  # cross-products
  XpX <- crossprod(X)
  XpZ <- crossprod(X, Z)
  Xpy <- crossprod(X, y)

  ZpX <- crossprod(Z, X)
  ZpZ <- crossprod(Z)
  Zpy <- crossprod(Z, y)

  if (model.base) {
    cat('starting to build & solve MME ...\n') 

# the Mixed Model
# y ~ 1mu + Za + Zd + e

# The mixed model equations (MME) solved are:

# | XpX    XpZ                  XpZ               | mu  | Xpy |
# | ZpX    ZpZ + Ainv*alpha1    ZpZ               | a   | Zpy |
# | ZpX    ZpZ                  ZpZ + Dinv*alpha2 | d   | Zpy |

    # create the rows of the coefficient matrix (C)
    c1 <- cbind(XpX, XpZ, XpZ)
    c2 <- cbind(ZpX, ZpZ + Ainv * alpha1, ZpZ)
    c3 <- cbind(ZpX, ZpZ, ZpZ + Dinv * alpha2)

    # the coefficient matrixsol
    C <- rbind(c1, c2, c3) # (p + 2m) x (p + 2m)
    RHS <- rbind(Xpy, Zpy, Zpy)

    # solving the MME
    sol <- solve.CHMperm(C = C, b = RHS)

    # the overall mean
    mu <- sol[1]
    # retrieve additive & non-additive Mendelian effects
    BLUP_add <- as.data.table(sol[seq(
        p + 1, m + p
    )])
    BLUP_dom <- as.data.table(sol[seq(
        p + m + 1, 2 * m + p 
    )])
    # merge additive & non-additive BLUPs
    BLUP <- data.table(BLUP_add, BLUP_dom)
    setnames(BLUP, c("BLUP_add", "BLUP_dom"))
    
    BLUP[, TreeID := as.double(pedigree[[1]])]
    trial[, TreeID := as.double(TreeID)]
    BLUP[, cross := trial.tmp[, .(TreeID, cross)][BLUP[, .(TreeID)],
        on = .(TreeID = TreeID)
    ][, cross]]
    BLUP[, BLUP_dom := ifelse(!is.na(cross), BLUP_dom, NA)]
    
    # computing the log-likelihood (suggested only for small datasets...!)
    if (Loglik) {
        loglik <- loglik(
            pedigree = pedigree,
            h2 = h2,
            H2 = H2,
            y = y,
            X = X,
            Z = Z,
            D = D,
            C = C
        )
    } else {
        loglik <- NULL
    }
    setcolorder(
        BLUP,
         c(3, 1, 4, 2)
    )
    lab <- names(BLUP)
        lab.m <- paste0(".m", n.mod)
        new.names <- c(
            lab[1],
            paste0(lab[2], lab.m), 
            lab[3], paste0(lab[4], lab.m) 
    )
    setnames(BLUP, new.names)
    cat('job done...! \n')
    if (gen < 3L) {
       return(list(
           BLUP = BLUP,
           mu = mu,
           loglik = loglik
       ))
    } else {
        return(list(
           BLUP = BLUP,
           mu = mu,
           F = F,
           loglik = loglik
        ))
    }
  } else {

      Zf <- Z %*% f
      XpZf <- crossprod(X, Zf)
      fpZpZf <- crossprod(Zf)
      ZpZf <- crossprod(Z, Zf)

      fpZpX <- crossprod(Zf, X)
      fpZpZ <- crossprod(Zf, Z)
      fpZpy <- crossprod(Zf, y)

# the Mixed Model
# y ~ 1mu + Zfb + Za + Zd + e

#  The mixed model equations (MME) solved are:

# | XpX     XpZf     XpZ                 XpZ                | mu | Xpy   
# | fpZpX   fpZpZf   fpZpZ               fpZpZ              | b  | fpZpy 
# | ZpX     ZpZf     ZpZ + Ainv*alpha1   ZpZ                | a  | Zpy   
# | ZpX     ZpZf     ZpZ                 ZpZ + Dinv*alpha2  | d  | Zpy   

      # create the rows of the coefficient matrix (C)
      c1 <- cbind(XpX, XpZf, XpZ, XpZ)
      c2 <- cbind(fpZpX, fpZpZf, fpZpZ, fpZpZ)
      c3 <- cbind(ZpX, ZpZf, ZpZ + Ainv * alpha1, ZpZ)
      c4 <- cbind(ZpX, ZpZf, ZpZ, ZpZ + Dinv * alpha2)

      # the coefficient matrix
      C <- rbind(c1, c2, c3, c4) # (2m+p+1)x(2m+p+1)
      RHS <- rbind(Xpy, fpZpy, Zpy, Zpy)

      # solving the MME
      sol <- solve.CHMperm(C = C, b = RHS)

      # retrieve the overall mean
      mu <- sol[1]
      # retrieve the coefficient for inbreeding depression ("id")
      id <- sol[2] 
      # estimate the effect of "id"
      dom0.id <- Zf * id 
      zeros.id <- zeros.matrix(m - n, 1)
      BLUE_id <- as.data.table(as.matrix(
          rbind(zeros.id, dom0.id)
      ))
      # retrieve breeding values
      BLUP_add <- as.data.table(sol[seq(
          p + 1 + 1, p + 1 + m
      )])
      # retrieve dominance deviations
      BLUP_dom <- as.data.table(sol[seq(
          p + 1 + m + 1, p + 1 + 2 * m 
      )])

      # merge additive & non-additive effects
      BLUP <- data.table(BLUP_add, BLUP_dom, BLUE_id)
      setnames(BLUP, c("BLUP_add", "BLUP_dom", "id"))

      BLUP[, TreeID := as.double(pedigree[[1]])]
      trial[, TreeID := as.double(TreeID)]
      BLUP[, cross := trial.tmp[, .(TreeID, cross)][BLUP[, .(TreeID)],
          on = .(TreeID = TreeID)
      ][, cross]]
      BLUP[, BLUP_dom := ifelse(!is.na(cross), BLUP_dom, NA)]

      # adjusting dominance effects by "id"
      BLUP[, BLUP_dom.adj := BLUP_dom + id]

      if (Loglik) {
          loglik <- loglik(
              pedigree = pedigree,
              h2 = h2,
              H2 = H2,
              X = X,
              y = y,
              Z = Z,
              D = D,
              C = C
          )
      } else {
          loglik <- NULL
      }
    setcolorder(
        BLUP,
         c(4, 1, 5, 2, 3, 6)
    )
      lab <- names(BLUP)
      lab.m <- paste0(".m", n.mod)

      new.names <- c(
          lab[1], 
          paste0(lab[2], lab.m), 
          lab[3],
          paste0(lab[4:6], lab.m)
      )
      setnames(BLUP, new.names)
      cat('job done...! \n')
      if (gen < 3L) {
          return(list(
              BLUP = BLUP,
              mu = mu,
              loglik = loglik
          ))
      } else {
          return(list(
              BLUP = BLUP,
              mu = mu,
              loglik = loglik,
              id = id,
              F = F
          ))
      }  

  }
}

 
 


