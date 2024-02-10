prediction.S3 <- function(pedigree,
                          ped.jacquard,
                          trial,
                          Fj = NULL,
                          Ainv = NULL,
                          init.ghost = NULL,
                          type = c("ide", "ks", "jq"),
                          ID = FALSE,
                          h2,
                          H2,
                          n.mod,
                          center = FALSE,
                          Loglik = FALSE) {
    library(Matrix, quietly = TRUE)
    library(data.table, quietly = TRUE)

    source("matrixMethods.r")
    source("makeKinship.nonadj.r")
    source("makeF.r")
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
    if (!("gen" %in% names(trial.tmp))) { # if not, recovery from pedigree
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
    # check if the F matrix is provided
    # if not, create it
    if (!is.null(Fj)) { # F matrix is not NULL
        if (inherits(Fj, "matrix")) { # check if Fj is a matrix
            if (!isSymmetric(Fj)) { # check if Fj is symmetric
                stop("matrix F must be symmetric!")
            } else { # Fj is a symmetric matrix
                F <- Matrix(Fj, sparse = TRUE) # create a sparse matrix
            }
        } else if (inherits(Fj, "dgCMatrix")) { # check if Fj is a dgCMatrix
            F <- as(Fj, "symmetricMatrix") # create a symmetric matrix
        } else if (!inherits(Fj, "dsCMatrix")) { # check if Fj is a dsCMatrix
            stop("F must be a symmetric matrix")
        }
    } else if (is.null(Fj) & type == "jq") {
        if (!is.null(init.ghost)) {
            F.out <- makeF( # F.out is a list
                pedigree = ped.jacquard,
                init.ghost = init.ghost
            )
        } else {
            F.out <- makeF(
                pedigree = ped.jacquard
            )
        }
        F <- F.out[["F"]] # the F relatedness matrix from Jacquard (1974)
        cond <- F.out[["cond"]] # condition number for matrix inversion
        if (cond > 4.5e15) {
            warning("Matrix 'F' is ill-conditioned")
        }
        Finv <- solve.CHMperm(F) # create the inverse of F matrix
        P11 <- makeM( # create the inverse of M matrix (adjusted)
            pedigree = pedigree,
            ginv = TRUE,
            adj = TRUE
        ) 
        cat("Jacquard ok! ...\n")
    } else if (is.null(Fj) & type == "ide") {
         Finv <- F <- makeF.ide( # create the F matrix as an identity matrix
             pedigree = pedigree
         )
         P11 <- makeM(
             pedigree = pedigree, # create the inverse of M matrix
             ginv = TRUE        # assuming that full-sib families are unrelated
         )
         cat("ide ok! ...\n")
    } else if (is.null(Fj) & type == "ks") {
        famped.ks <- makeKinship.nonadj( # famped.ks is a list
            pedigree = pedigree,
            familyPed = TRUE
        )
        # the F relatedness matrix from kinship coefficients (Cokerham 1954)
        F <- famped.ks[["F"]]
        # compute and check the condition number of the F matrix
        cnd <- famped.ks[["cond"]]
        if (cnd > 4.5e15) {
            warning("Matrix 'F' is ill-conditioned")
        } # create the inverse of F matrix
        Finv <- solve.CHMperm(F)
        P11 <- makeM( # create the inverse of the M matrix non-adjusted 
            pedigree = pedigree,
            ginv = TRUE,
        )
        cat("ks ok! ...\n")
    } else if (is.null(Fj) & type == "") {
        stop(paste0("type=", type, " is not an option"))
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
    # building the incidence matrices
    n <- dim(trial)[1]
    X <- ones.matrix(n, 1)
    Z <- Z.mat(pedigree)
    Zdom <- Z.dom(pedigree)[["Zdom"]]
    # Construction of block matrices
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    #'                                                                           '#
    #'   Dblock= | D11 D12  |     Dblock.inv= | P11 P12 |  = | Dinv11 Dinv12 |   '#
    #'           | D21 D22  |                 | P21 P22 |    | Dinv21 D22inv |   '#
    #'                                                                           '#
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    # create the he inverse of block matrices by using the properties
    # of the  inverse of a partioned matrix (Searle 1971)
    P12 <- -P11 %*% Zdom
    P21 <- t(P12)
    P22 <- t(Zdom) %*% P11 %*% Zdom + 4 * Finv # (8x28)x(28x28)x(28x8)
    y <- as.vector(trial.tmp[, pheno])
    # number of full-sib families (q) and number of fixed effects (p)
    q <- dim(Zdom)[2]
    p <- dim(X)[2]
    # number of individuals (m) in the pedigree
    m <- dim(pedigree)[1]
    # setting sigmas
    sigma2.add <- h2
    sigma2.dom <- H2 - h2
    sigma2.sca <- sigma2.dom / 4
    sigma2.e <- 1 - H2
    alpha1 <- sigma2.e / sigma2.add
    alpha2 <- sigma2.e / sigma2.dom
    # build null matrices for setting the MME
    zeros1 <- zeros.matrix(p, q)
    zeros2 <- zeros.matrix(m, q)
    zeros3 <- zeros.matrix(q, 1)
    zeros1p <- t(zeros1)
    zeros2p <- t(zeros2)
    zeros3p <- t(zeros3)
    # cross-products
    XpX <- crossprod(X)
    XpZ <- crossprod(X, Z)
    Xpy <- crossprod(X, y)
    ZpX <- crossprod(Z, X)
    ZpZ <- crossprod(Z)
    Zpy <- crossprod(Z, y)
    if (model.base) {
        cat("starting to build & solve MME ...\n")
        # the Mixed Model
        # y ~ 1mu + Za + Zd + e
        # The mixed model equations (MME) solved are:

        # | XpX   XpZ                  XpZ                     0       | mu  | Xpy  |
        # | ZpX   ZpZ  + Ainv*alpha1   ZpZ                     0       | a   | Zpy  |
        # | ZpX   ZpZ                  ZpZ  + P11*alpha2  P12*alpha2   | d   | Zpy  |
        # |  0     0                   P21*alpha2         P22*alpha2   | sca | 0    |

        # create the rows of the coefficient matrix (C)
        c1 <- cbind(XpX, XpZ, XpZ, zeros1) # 1x(2m+q+p)
        c2 <- cbind(ZpX, ZpZ + Ainv * alpha1, ZpZ, zeros2) # (m, p+2m+q)
        c3 <- cbind(ZpX, ZpZ, ZpZ + P11 * alpha2, P12 * alpha2) # (m, p+2m+q)
        c4 <- cbind(zeros1p, zeros2p, P21 * alpha2, P22 * alpha2) # (q, p+2m+q)
        # the coefficient matrix
        C <- rbind(c1, c2, c3, c4) # (2m+p+q)x(2m+p+q)
        RHS <- rbind(Xpy, Zpy, Zpy, zeros3) # (p+2m+q)x1
        # solving the MME
        sol <- solve.CHMperm(C = C, b = RHS)
        # retrieve the overall mean
        mu <- sol[1]
        # retrieve breeding values
        BLUP_add <- as.data.table(sol[seq(
            p + 1, m + p
        )])
        # retrieve dominance deviations
        BLUP_dom <- as.data.table(sol[seq(
            m + p + 1, 2 * m + p
        )])
        BLUP <- data.table(BLUP_add, BLUP_dom)
        setnames(BLUP, c("BLUP_add", "BLUP_dom"))
        BLUP[, TreeID := as.double(pedigree[[1]])]
        # retrieve sca effects
        BLUP_sca <- as.data.table(sol[seq(
            p + 2 * m + 1, p + 2 * m + q
        )])
        setnames(BLUP_sca, "BLUP_sca")
        BLUP_sca[, cross := as.vector(rownames(F))]
        # adding cross field to BLUP
        trial[, TreeID := as.double(TreeID)]
        BLUP[, cross := trial.tmp[, .(TreeID, cross)][BLUP[, .(TreeID)],
            on = .(TreeID = TreeID)
        ][, cross]]
        # look up the sca effects for each cross
        BLUP[, BLUP_sca := BLUP_sca[, .(cross, BLUP_sca)][BLUP[, .(cross)],
            on = .(cross = cross)
        ][, BLUP_sca]]
        BLUP[, BLUP_dom := ifelse(!is.na(BLUP_sca), BLUP_dom, NA)]

        # compute the REML log-likelihood (suggested only for small datasets...!)
        if (Loglik) {
            # recovering the dominance relatedness matrix (D11) by using
            # the properties of the inverse of a partioned matrix (Searle 1971)
            D11 <- as(
                solve(P11 - P12 %*% solve(P22) %*% P21),
                "sparseMatrix"
            )
            loglik <- loglik(
                pedigree = pedigree,
                h2 = h2,
                H2 = H2,
                X = X,
                y = y,
                Z = Z,
                F = F,
                D = D11,
                C = C
            )
        } else {
            loglik <- NULL
        }
        setcolorder(
            BLUP,
            c(3, 1, 4, 5, 2)
        )
        lab <- names(BLUP)
        lab.m <- paste0(".m", n.mod)
        new.names <- c(
            lab[1],
            paste0(lab[2], lab.m), 
            lab[3], paste0(lab[4:5], lab.m) 
        )
        setnames(BLUP, new.names)
        cat("job done...! \n")
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

        zeros4 <- zeros.matrix(1, q)
        zeros4p <- zeros.matrix(q, 1)
        # the Mixed Model
        # y ~ 1mu + Zfb + Za + Zd + e
        # The mixed model equations (MME) solved are:
      
        # | XpX   XpZf    XpZ               XpZ                  0       | mu  | Xpy  |
        # | fpZpX fpZpZf  fpZpZ             fpZpZ                0       | b   | fpZpy|
        # | ZpX   ZpZf    ZpZ + Ainv*alpha1 ZpZ                  0       | a   | Zpy  |
        # | ZpX   ZpZf    ZpZ               ZpZ + P11*alpha2  P12*alpha2 | d   | Zpy  |
        # |  0     0      0                 P21*alpha2        P22*alpha2 | sca | 0    |

        # create the rows of the coefficient matrix (C)
        c1 <- cbind(XpX, XpZf, XpZ, XpZ, zeros1) # px(p+1+2m+q)
        c2 <- cbind(fpZpX, fpZpZf, fpZpZ, fpZpZ, zeros4) # (1, p+1+2m+q)
        c3 <- cbind(ZpX, ZpZf, ZpZ + Ainv * alpha1, ZpZ, zeros2) # (m, p+1+2m+q)
        c4 <- cbind(ZpX, ZpZf, ZpZ, ZpZ + P11 * alpha2, P12 * alpha2) # (m, p+1+2m+q)
        c5 <- cbind(zeros1p, zeros4p, zeros2p, P21 * alpha2, P22 * alpha2) # (q, p+1+2m+q)

        # the coefficient matrix
        C <- rbind(c1, c2, c3, c4, c5) # (2m+p+q+1)x(2m+p+q+1)
        RHS <- rbind(Xpy, fpZpy, Zpy, Zpy, zeros4p)
        # solving the MME
        sol <- solve.CHMperm(C = C, b = RHS)
        # the overall mean
        mu <- sol[1]
        # retrieve the coefficient for inbreeding depression ("id")
        id <- sol[2]
        # estimate the effect of "id"
        dom0.id <- Zf * id
        zeros.id <- zeros.matrix(m - n, 1)
        BLUE_id <- as.data.table(as.matrix(
            rbind(zeros.id, dom0.id)
        ))
        # retrieve breeding values (bv)
        BLUP_add <- as.data.table(sol[seq(
            p + 1 + 1, p + 1 + m
        )])
        # retrieve dominance deviations (dd)
        BLUP_dom <- as.data.table(sol[seq(
            p + 1 + m + 1, p + 1 + 2 * m
        )])
        # merge bv, dd and id
        BLUP <- data.table(BLUP_add, BLUP_dom, BLUE_id)
        setnames(BLUP, c("BLUP_add", "BLUP_dom", "id"))
        # retrieve sca effects
        BLUP_sca <- data.table(BLUP_sca = sol[seq(
            p + 1 + 2 * m + 1, p + 1 + 2 * m + q
        )])
        BLUP_sca[, cross := as.vector(rownames(F))]
        # merge all effects
        BLUP[, TreeID := as.double(pedigree[[1]])]
        trial[, TreeID := as.double(TreeID)]
        BLUP[, cross := trial.tmp[, .(TreeID, cross)][BLUP[, .(TreeID)],
            on = .(TreeID = TreeID)
        ][, cross]]
        BLUP[, BLUP_sca := BLUP_sca[, .(cross, BLUP_sca)][BLUP[, .(cross)],
            on = .(cross = cross)
        ][, BLUP_sca]]
        BLUP[, BLUP_dom := ifelse(!is.na(BLUP_sca), BLUP_dom, NA)]
        # adjusting dominance effects by "id"
        BLUP[, BLUP_dom.adj := BLUP_dom + id]
        BLUP[, BLUP_sca.adj := BLUP_sca + id]
        if (Loglik) {
            D11 <- as(
                solve(P11 - P12 %*% solve(P22) %*% P21),
                "sparseMatrix"
            )
            loglik <- loglik(
                pedigree = pedigree,
                h2 = h2,
                H2 = H2,
                X = X,
                y = y,
                Z = Z,
                F = F,
                D = D11,
                C = C
            )
        } else {
            loglik <- NULL
        }
        setcolorder(
            BLUP,
            c(4, 1, 5, 6, 3, 8, 2, 7)
        )
        lab <- names(BLUP)
        lab.m <- paste0(".m", n.mod)
        new.names <- c(
            lab[1],
            paste0(lab[2], lab.m),
            lab[3], paste0(lab[4:8], lab.m)  
        )
        setnames(BLUP, new.names)
        cat("job done...! \n")
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
