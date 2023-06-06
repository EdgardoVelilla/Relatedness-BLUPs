prediction.S2 <- function(pedigree,
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
        Minv <- makeM( # create the inverse of M matrix adjusted by inbreeding
            pedigree = pedigree,
            ginv = TRUE,
            adj = TRUE
        ) 
        cat("Jacquard ok! ...\n")
    } else if (is.null(Fj) & type == "ide") {
        Finv <- F <- makeF.ide( # create the F matrix as an identity matrix
            pedigree = pedigree
        )
        Minv <- makeM( # create the inverse of M matrix assuming that
             pedigree = pedigree, # full-sib families are unrelated
             ginv = TRUE        
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
        Minv <- makeM( # create the inverse of the M matrix non-adjusted 
            pedigree = pedigree, 
            ginv = TRUE,
        )
        cat("ks ok! ...\n")
    } else if (is.null(Fj) & type == "") {
        stop(paste0("type=", type, " is not an option"))
    }

    if (is.null(Ainv)) { # check if Ainv matrix is NULL ?
        cat("starting to build Ainv ...\n")
        Ainv <- Ainverse(pedigree) # # if so, create it!
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

    # incidence matrices
    n <- dim(trial)[1]
    X <- ones.matrix(n, 1)
    Z <- Z.mat(pedigree)
    U <- Z.fam(trial.tmp)
    Zdom <- Z.dom(pedigree)[["Zdom"]]

    y <- as.vector(trial.tmp[, pheno])
    # number of full-sib families
    q <- dim(Zdom)[2]
    p <- dim(X)[2]
    m <- dim(pedigree)[1]

    # setting sigmas
    sigma2.add <- h2
    sigma2.dom <- H2 - h2
    sigma2.sca <- sigma2.dom / 4
    sigma2.e <- 1 - H2

    alpha1 <- sigma2.e / sigma2.add
    alpha2 <- sigma2.e / sigma2.dom
    alpha3 <- sigma2.e / sigma2.sca

    # cross-products
    XpX <- crossprod(X)
    XpZ <- crossprod(X, Z)
    XpU <- crossprod(X, U)
    Xpy <- crossprod(X, y)

    ZpX <- crossprod(Z, X)
    ZpZ <- crossprod(Z)
    ZpU <- crossprod(Z, U)
    Zpy <- crossprod(Z, y)

    UpX <- crossprod(U, X)
    UpZ <- crossprod(U, Z)
    UpU <- crossprod(U)
    Upy <- crossprod(U, y)

    if (model.base) {
        cat("starting to build & solve MME ...\n")

        # the Mixed Model
        # y ~ 1mu + Za + Usca + Zm + e

        # The mixed model equations (MME) solved are:

        # | XpX  XpZ                XpU                XpZ               | mu  | Xpy |
        # | ZpX  ZpZ + Ainv*alpha1  ZpU                ZpZ               | a   | Zpy |
        # | UpX  UpZ                UpU + Finv*alpha3  UpZ               | sca | Upy |
        # | ZpX  ZpZ                ZpU                ZpZ + Minv*alpha2 | m   | Zpy |

        # create the rows of the coefficient matrix (C)
        c1 <- cbind(XpX, XpZ, XpU, XpZ)
        c2 <- cbind(ZpX, ZpZ + Ainv * alpha1, ZpU, ZpZ)
        c3 <- cbind(UpX, UpZ, UpU + Finv * alpha3, UpZ)
        c4 <- cbind(ZpX, ZpZ, ZpU, ZpZ + Minv * alpha2)

        # the coefficient matrix
        C <- rbind(c1, c2, c3, c4) # (2m+p+q)x(2m+p+q)
        RHS <- rbind(Xpy, Zpy, Upy, Zpy)

        # solving the MME
        sol <- solve.CHMperm(C = C, b = RHS)

        # the overall mean
        mu <- sol[1]
        # retrieve additive & non-additive Mendelian effects
        BLUP_add <- as.data.table(sol[seq(
            p + 1, m + p
        )])
        BLUP_mend <- as.data.table(sol[seq(
            p + m + q + 1, 2 * m + p + q
        )])
        BLUP <- data.table(BLUP_add, BLUP_mend)
        setnames(BLUP, c("BLUP_add", "BLUP_mend"))
        # retrieve sca effects
        BLUP_sca <- data.table(BLUP_sca = sol[seq(
            p + m + 1, p + m + q
        )])
        BLUP_sca[, cross := as.vector(rownames(Finv))]
        # merge all effects
        BLUP[, TreeID := as.double(pedigree[[1]])]
        trial[, TreeID := as.double(TreeID)]
        BLUP[, cross := trial.tmp[, .(TreeID, cross)][BLUP[, .(TreeID)],
            on = .(TreeID = TreeID)
        ][, cross]]
        BLUP[, BLUP_sca := BLUP_sca[, .(cross, BLUP_sca)][BLUP[, .(cross)],
            on = .(cross = cross)
        ][, BLUP_sca]]
        BLUP[, BLUP_dom := ifelse(!is.na(BLUP_sca), BLUP_mend + BLUP_sca, NA)]
        BLUP[, BLUP_mend := NULL]

        # compute the log-likelihood (suggested only for small pedigrees...!)
        if (Loglik) {
            loglik <- loglik(
                pedigree = pedigree,
                h2 = h2,
                H2 = H2,
                X = X,
                y = y,
                Z = Z,
                F = F,
                U = U,
                C = C
            )
        } else {
            loglik <- NULL
        }
        setcolorder(
            BLUP,
            c(2, 1, 3:5)
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
        UpZf <- crossprod(U, Zf)

        fpZpX <- crossprod(Zf, X)
        fpZpZ <- crossprod(Zf, Z)
        fpZpU <- crossprod(Zf, U)
        fpZpy <- crossprod(Zf, y)

        # the Mixed Model
        # y ~ 1mu + Zfb + Za + Usca + Zm + e

        #  The mixed model equations (MME) solved are:
        # | XpX   XpZf   XpZ                XpU                XpZ                | mu | Xpy
        # | fpZpX fpZpZf fpZpZ              fpZpU              fpZpZ              | b  | fpZpy
        # | ZpX   ZpZf   ZpZ + Ainv*alpha1  ZpU                ZpZ                | a  | Zpy
        # | UpX   UpZf   UpZ                UpU + Finv*alpha3  UpZ                | sca| Upy
        # | ZpX   ZpZf   ZpZ                ZpU                ZpZ + Minv*alpha2  | m  | Zpy

        # create the rows of the coefficient matrix (C)
        c1 <- cbind(XpX, XpZf, XpZ, XpU, XpZ)
        c2 <- cbind(fpZpX, fpZpZf, fpZpZ, fpZpU, fpZpZ)
        c3 <- cbind(ZpX, ZpZf, ZpZ + Ainv * alpha1, ZpU, ZpZ)
        c4 <- cbind(UpX, UpZf, UpZ, UpU + Finv * alpha3, UpZ)
        c5 <- cbind(ZpX, ZpZf, ZpZ, ZpU, ZpZ + Minv * alpha2)

        # the coefficient matrix
        C <- rbind(c1, c2, c3, c4, c5) # (2m+p+q+1)x(2m+p+q+1)
        RHS <- rbind(Xpy, fpZpy, Zpy, Upy, Zpy)

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
        # retrieve breeding values
        BLUP_add <- as.data.table(sol[seq(
            p + 2, m + p + 1
        )])
        # retrieve BLUPs from non-additive Mendelian effects
        BLUP_mend <- as.data.table(sol[seq(
            p + m + q + 2, 2 * m + p + q + 1
        )])
        # merge additive & non-additive BLUPs
        BLUP <- data.table(BLUP_add, BLUP_mend, BLUE_id)
        setnames(BLUP, c("BLUP_add", "BLUP_mend", "id"))
        # retrieve sca effects
        BLUP_sca <- data.table(BLUP_sca = sol[seq(
            p + m + 2, p + m + q + 1
        )])
        BLUP_sca[, cross := as.vector(rownames(Finv))]

        # merge all effects
        BLUP[, TreeID := as.double(pedigree[[1]])]
        trial[, TreeID := as.double(TreeID)]
        BLUP[, cross := trial.tmp[, .(TreeID, cross)][BLUP[, .(TreeID)],
            on = .(TreeID = TreeID)
        ][, cross]]
        BLUP[, BLUP_sca := BLUP_sca[, .(cross, BLUP_sca)][BLUP[, .(cross)],
            on = .(cross = cross)
        ][, BLUP_sca]]
        BLUP[, BLUP_dom := ifelse(!is.na(BLUP_sca), BLUP_mend + BLUP_sca, NA)]
        BLUP[, BLUP_mend := NULL]

        # adjusting dominance effects by "id"
        BLUP[, BLUP_dom.adj := BLUP_dom + id]
        BLUP[, BLUP_sca.adj := BLUP_sca + id]

        if (Loglik) {
            loglik <- loglik(
                pedigree = pedigree,
                h2 = h2,
                H2 = H2,
                y = y,
                X = X,
                Z = Z,
                F = F,
                U = U,
                C = C
            )
        } else {
            loglik <- NULL
        }
        setcolorder(
            BLUP,
            c(3, 1, 4:5, 2, 8, 6:7)
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
