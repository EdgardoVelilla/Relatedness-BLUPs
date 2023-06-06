pred.ASReml <- function(pedigree,
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

    library(asreml, quietly = TRUE) # version 4.1.0.176
    asreml.options(gammaPar = TRUE)
    library(data.table, quietly = TRUE) # version 1.14.6
    library(Matrix, quietly = TRUE) # version 1.5-3
    library(stringr) # version 1.5.0

    source("forASReml.r")
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

    # check if trial has a column named "Xtype"
    if (!("Xtype" %in% names(trial.tmp))) { # if not, create it!
        trial.tmp[, Xtype := ifelse((mum != 0L | dad != 0L),
            "CP", "OP"
        )]
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
        # create the inverse of M matrix adjusted by inbreeding
        Minv <- makeM( 
            pedigree = pedigree,
            ginv = TRUE,
            adj = TRUE
        ) 
        cat("Jacquard ok! ...\n")
    } else if (is.null(Fj) & type == "ide") {
        Finv <- F <- F.ide( # create the F matrix as an identity matrix
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
        Minv <- makeM( # create the inverse of M matrix assuming that
             pedigree = pedigree, # full-sib families are unrelated
             ginv = TRUE        
         )
        cat("ks ok! ...\n")
    } else if (is.null(Fj) & type == "") {
        stop(paste0("type=", type, " is not an option"))
    }
    
    if (is.null(Ainv)) { # check if Ainv matrix is NULL ?
        cat("starting to build Ainv ...\n")
        Ainv <- Ainverse(pedigree) # if so, create it!
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

    # converting model terms into factors
    trial.tmp[, TreeID2 := TreeID]
    terms <- c("TreeID", "cross", "TreeID2", "Xtype")
    trial.tmp[, (terms) := lapply(
        .SD,
        function(x) factor(x, levels = unique(x))
    ), .SDcols = terms]

    # converting the inverse matrices (Ainv, Finv, Minv) to ASReml-R’s giv format
    Ainv.as <- forASReml(Ainv, Ginv = TRUE)
    assign("Ainv.as", Ainv.as,
        envir = globalenv()
    )
    Finv.as <- forASReml(Finv, Ginv = TRUE)
    assign("Finv.as", Finv.as,
        envir = globalenv()
    )
    Minv.as <- forASReml(Minv, Ginv = TRUE)
    assign("Minv.as", Minv.as,
        envir = globalenv()
    )

    # setting sigmas
    sigma2.add <- h2
    sigma2.dom <- H2 - h2
    sigma2.sca <- sigma2.dom / 4
    sigma2.e <- 1 - H2

    # setting gammas
    gamma.add <- sigma2.add / sigma2.e
    gamma.sca <- sigma2.sca / sigma2.e
    gamma.dom <- sigma2.dom / sigma2.e
    gammas <- c(gamma.add, gamma.sca, gamma.dom, 1)

    # setting constraints
    constraints <- c("F", "F", "F", "F")

    # ASReml-R works with data.frames
    trial.df <- as.data.frame(copy(trial.tmp))
    assign("trial.df", trial.df,
        envir = globalenv()
    )

    # setting the mixed-model
    model_pheno <- "asreml(fixed = pheno ~ "
    fixed_mu <- "1,"
    fixed_mu_f <- "1 + f,"

    # model for trial with only CP families
    random_effects_CP <-
        "random = ~vm(TreeID, Ainv.as) + ~vm(cross, Finv.as) + ~vm(TreeID2, Minv.as),
                     data = trial.df, na.action = na.method(x = 'include', y='include'),"

    # model for trial with OP & CP families
    random_effects_OP_CP <-
        "random = ~vm(TreeID, Ainv.as) + ~at(Xtype, 'CP'):vm(cross, Finv.as) + ~vm(TreeID2, Minv.as),
                     data = trial.df, na.action = na.method(x = 'include', y='include'),"

    if (!trial[, "OP" %in% Xtype]) { # trial with only CP families
        if (model.base) { # pheno ~ mu !r TreeID cross TreeID2
            model <- paste(
                model_pheno, fixed_mu,
                random_effects_CP
            )
        } else { # pheno ~ mu f !r TreeID cross TreeID2
            model <- paste(
                model_pheno, fixed_mu_f,
                random_effects_CP
            )
        }
    } else if (trial[, "OP" %in% Xtype]) { # trial with OP & CP families
        if (model.base) { # pheno ~ mu !r TreeID at(Xtype,2).cross TreeID2
            model <- paste(
                model_pheno, fixed_mu,
                random_effects_OP_CP
            )
        } else { # pheno ~ mu f !r TreeID at(Xtype,2).cross TreeID2
            model <- paste(
                model_pheno, fixed_mu_f,
                random_effects_OP_CP
            )
        }
    }

    # setting the starting values
    temp.gam <- eval(parse(text = paste(
        model,
        "start.values = TRUE)"
    )))$vparameters.table

    # setting the constraints & gammas
    temp.gam[c(1:4), "Constraint"] <- constraints
    temp.gam[c(1:4), "Value"] <- gammas

    # run the BLUP analysis
    model.asr <- eval(parse(text = paste(
        model,
        "G.param = temp.gam, R.param = temp.gam)"
    )))

    # BLUEs's extraction
    BLUEs <- coef(model.asr)$fixed
    mu <- BLUEs[rownames(BLUEs) == "(Intercept)"]

    # BLUPs's extraction
    BLUPs <- as.data.table(
        as.matrix(coef(
            model.asr
        )[["random"]]),
        keep.rownames = TRUE
    )
    BLUPs[, TreeID := str_extract(rn, "(?<=_).*")]
    BLUPs[, factor := sub("^(.*?)_.*", "\\1", rn)]
    setnames(BLUPs, "effect", "BLUP")

    # setting a template to be used for predictions
    prediction <- copy(pedigree[, c(1:3)])
    # add an cross field to the template (prediction)
    prediction[
        TreeID %in% trial[, TreeID],
        cross := trial.tmp[, cross]
    ]
    # convert TreeID to character for lookups
    prediction[, TreeID := as.character(TreeID)]

    # retrieve additive effects
    BLUP_add <- BLUPs[
        factor == "vm(TreeID, Ainv.as)",
        .(TreeID, BLUP)
    ]
    setnames(BLUP_add, "BLUP", "BLUP_add")
    # lookup BLUP_add for prediction
    prediction[, BLUP_add := BLUP_add[
        ,
        .(TreeID, BLUP_add)
    ][prediction[, .(TreeID)],
        on = .(TreeID = TreeID)
    ][, BLUP_add]]

    # retrieve Mendelian effects
    BLUP_mend <- BLUPs[
        factor == "vm(TreeID2, Minv.as)",
        .(TreeID, BLUP)
    ]
    setnames(BLUP_mend, "BLUP", "BLUP_mend")
    # lookup BLUP_mend for prediction
    prediction[, BLUP_mend := BLUP_mend[
        ,
        .(TreeID, BLUP_mend)
    ][prediction[, .(TreeID)],
        on = .(TreeID = TreeID)
    ][, BLUP_mend]]

    # retrieve sca effects
    if (!trial[, "OP" %in% Xtype]) { # trial with only CP families
        BLUP_sca <- BLUPs[
            factor == "vm(cross, Finv.as)",
            .(TreeID, BLUP)
        ]
    } else { # trial with OP & CP families
        BLUP_sca <- BLUPs[
            factor == "at(Xtype, CP):vm(cross, Finv.as)",
            .(TreeID, BLUP)
        ]
    }
    setnames(BLUP_sca, c("cross", "BLUP_sca"))
    BLUP_sca[, cross := as.factor(cross)]
    # lookup BLUP_sca for prediction
    prediction[, BLUP_sca := BLUP_sca[
        ,
        .(cross, BLUP_sca)
    ][prediction[, .(cross)],
        on = .(cross = cross)
    ][, BLUP_sca]]

    # computing dominance effects as BLUP_sca + BLUP_mend
    prediction[, BLUP_dom := ifelse(!is.na(BLUP_sca),
        BLUP_mend + BLUP_sca, NA
    )]

    # remove some redundant columns
    col.remove <- c("BLUP_mend", "mum", "dad")
    prediction[, c(col.remove) := NULL]

    # retrieve the REML log-likelihood
    if (Loglik) {
        loglik <- model.asr[["loglik"]]
    } else {
        loglik <- NULL
    }

    if (model.base) {
        setcolorder(
            prediction, c(1, 3, 2, 4:5)
        )
        lab <- names(prediction)
        lab.m <- paste0(".m", n.mod)
        new.names <- c(
            lab[1],
            paste0(lab[2], lab.m), 
            lab[3], paste0(lab[4:5], lab.m) 
        )
        setnames(prediction, new.names)
        cat("done...! \n")
        if (gen < 3L) {
            return(list(
                BLUP = prediction,
                mu = mu,
                loglik = loglik
            ))
        } else {
            return(list(
                BLUP = prediction,
                mu = mu,
                loglik = loglik,
                F = F
            ))
        }
    } else {
        # get the coefficient for inbreeding depression ("id")
        id <- BLUEs[rownames(BLUEs) == "f"]
        # estimate the effect of "id"
        Z <- Z.mat(pedigree)
        Zf <- Z %*% f
        dom0.id <- Zf * id
        n <- dim(trial)[1]
        m <- dim(pedigree)[1]
        zeros.id <- zeros.matrix(m - n, 1)
        BLUE_id <- as.data.table(as.matrix(
            rbind(zeros.id, dom0.id)
        ))
        BLUE_id[, TreeID := as.character(pedigree[, TreeID])]
        setnames(BLUE_id, 1, "id")

        # lookup BLUE_id for prediction
        prediction[, id := BLUE_id[
            ,
            .(TreeID, id)
        ][prediction[, .(TreeID)],
            on = .(TreeID = TreeID)
        ][, id]]

        # adjusting dominance effects by "id"
        prediction[, BLUP_dom.adj := BLUP_dom + id]
        prediction[, BLUP_sca.adj := BLUP_sca + id]
        setcolorder(prediction, c(1, 3, 2, 4, 6, 8, 5, 7))

        lab <- names(prediction)
        lab.m <- paste0(".m", n.mod)
        new.names <- c(
            lab[1],
            paste0(lab[2], lab.m),
            lab[3], paste0(lab[4:8], lab.m)  
        )
        setnames(prediction, new.names)
        cat("done...! \n")
        if (gen < 3L) {
            return(list(
                BLUP = prediction,
                mu = mu,
                loglik = loglik
            ))
        } else {
            return(list(
                BLUP = prediction,
                mu = mu,
                loglik = loglik,
                id = id,
                F = F
            ))
        }
    }
}
