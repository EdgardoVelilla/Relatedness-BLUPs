pred.ASRemlSA <- function(path = NULL, # provides the path to the ASReml executable file
                          workspace = NULL, # c("s7", "s8", "s9"), # 1024 Mb, 2048 Mb, 4096 Mb
                          job.file = NULL, # the name of the .as file
                          Pop.type = c("outbred", "inbred"),
                          pedigree,
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

    library(data.table, quietly = TRUE) # version 1.14.6
    library(Matrix, quietly = TRUE) # version 1.5-3
    library(stringr) # version 1.5.0
    options(scipen = 999)

    source("forASReml.r")
    source("matrixMethods.r")
    source("makeKinship.nonadj.r")
    source("makeF.r")

    if (is.null(path)) { # <-- ASReml.exe
        path <- paste0('"c:/Program Files/asreml4/bin/ASReml.exe"')
    }

    if (is.null(workspace)) {
        workspace <- "s9" # <-- 4096 Mb
    }

    # check if the job.file (.as) is provided
    # if not, create it
    if (is.null(job.file)) {
        if (Pop.type == "outbred") { # <-- "BLUP_Analysis_2023_outbred.as"
            job.file <- paste0(
                "BLUP_Analysis_", format(Sys.Date(), "%Y"),
                "_outbred", ".as"
            )
        } else if (Pop.type == "inbred") { # <-- "BLUP_Analysis_2023_inbred.as"
            job.file <- paste0(
                "BLUP_Analysis_", format(Sys.Date(), "%Y"),
                "_inbred", ".as"
            )
        }
    }

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
        Minv <- makeM( # create the inverse of M matrix adjusted by inbreeding
            pedigree = pedigree,
            ginv = TRUE,
            adj = TRUE
        ) 
        cat("Jacquard ok! ...\n")
    } else if (is.null(Fj) & type == "ide") { # full-sib families are unrelated
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

    # check if the inverse of A matrix is provided
    # if not, create it
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
    }

    # ensure that the coefficient matrix (C) is non-singular
    if (!ID || (ID && all(f == 0))) {
        model.base <- TRUE
    } else {
        model.base <- FALSE
    }

    # convert the inverse of relatedness matrices (Ainv, Finv, Minv)
    # to the ASReml format by forASReml() function
    Ainv.as <- forASReml(Ainv, # Ainv matrix
        Ginv = TRUE,
        colnames = c("row", "col", "Ainv")
    )
    Finv.as <- forASReml(Finv, # Finv matrix
        Ginv = TRUE,
        colnames = c("row", "col", "Finv")
    )
    Minv.as <- forASReml(Minv, # Minv matrix
        Ginv = TRUE, ,
        colnames = c("row", "col", "Minv")
    )

    # create an additional field for Mendelian term
    trial.tmp[, TreeID2 := TreeID]

    # select the trial's relevant columns and save it as a .csv file
    if (!ID) { # pheno ~ mu !r TreeID cross TreeID2
        trial.csv <- trial.tmp[,
            .(TreeID, TreeID2, mum, dad, gen, Xtype, cross, pheno)
        ]
    } else { # pheno ~ mu f !r TreeID cross TreeID2
        trial.csv <- trial.tmp[,
            .(TreeID, TreeID2, mum, dad, gen, f, Xtype, cross, pheno)
        ]
    }
    write.csv(trial.csv, # trial.csv file
        file = paste("trial.csv", sep = ""),
        row.names = FALSE, na = "."
    )
    # save the pedigree at the current working directory
    ped <- pedigree[, c(1:3)]
    write.csv(ped, # pedigree.csv file
        file = paste("pedigree.csv", sep = ""),
        row.names = FALSE, na = "."
    )
    # save the inverse of the relatedness matrices (Ainv, Finv, Minv) at the
    # current working directory to .giv file for 'ASReml SA' to use
    Ainv.giv <- as.data.table(as.matrix(Ainv.as))
    fwrite(Ainv.giv, # Ainv.giv file
        file = "Ainv.giv", quote = FALSE,
        sep = " ", row.names = FALSE, col.names = TRUE
    )

    Finv.giv <- as.data.table(as.matrix(Finv.as))
    fwrite(Finv.giv, # Finv.giv file
        file = "Finv.giv", quote = FALSE,
        sep = " ", row.names = FALSE, col.names = TRUE
    )

    Minv.giv <- as.data.table(as.matrix(Minv.as))
    fwrite(Minv.giv, # Minv.giv file
        file = "Minv.giv", quote = FALSE,
        sep = " ", row.names = FALSE, col.names = TRUE
    )

    if (!trial[, "OP" %in% Xtype]) { # trial with only CP families
        if (model.base) {
            part <- "1" # pheno ~ mu !r TreeID cross TreeID2
        } else {
            part <- "2" # pheno ~ mu f !r TreeID cross TreeID2
        }
    } else { # trial with OP & CP families
        if (model.base) {
            part <- "3" # pheno ~ mu !r TreeID at(Xtype,2).cross TreeID2
        } else {
            part <- "4" # pheno ~ mu f !r TreeID at(Xtype,2).cross TreeID2
        }
    }

    # the command line to run the job
    command <- paste0(
        path, " ", "-", workspace, " ", job.file, " ", part
    )
    # command <- paste0(
    #   '"c:/Program Files/asreml4/bin/ASReml.exe"', " ",
    #   "-", workspace, " ", job.file, " ", part
    # )

    # calling to 'ASReml SA' to run the BLUP analysis
    run <- system(
        command,
        wait = TRUE, intern = TRUE
    )
    # get the name of the .as file
    name <- stringr::str_extract(job.file, "(.*?)\\.")
    # read the .asr file
    asr <- readLines(paste0(name, "asr"), -1)
    loglik_line <- grep("LogL=",
        readLines(paste0(name, "asr"), -1),
        value = TRUE
    )
    # get the REML log-likelihood from asr file
    loglik <- as.double(str_extract(
        loglik_line,
        "(?<=LogL=\\s)[-+]?\\d*\\.?\\d+"
    ))
    # read the .sln file
    sln <- fread(paste0(name, "sln"))
    # get the overall mean
    mu <- sln[Model_Term == "mu", Effect]
    # get the coefficient for inbreeding depression ("id")
    if (ID) { # y ~ 1mu + Zfb + Za + Usca + Zm + e
        id <- sln[Model_Term == "f", Effect]
    } else { # y ~ 1mu + Za + Usca + Zm + e
        id <- 0
    }

    # setting a template to be used for predictions
    prediction <- copy(pedigree[, c(1:3)])
    # add an cross field to the template (prediction)
    prediction[
        TreeID %in% trial[, TreeID],
        cross := trial.tmp[, cross]
    ]
    # convert TreeID to character for lookup
    prediction[, TreeID := as.character(TreeID)]

    # retrieve additive effects
    BLUP_add <- sln[Model_Term == "TreeID", .(Level, Effect)]
    setnames(BLUP_add, c("TreeID", "BLUP_add"))
    # lookup BLUP_add for prediction
    prediction[, BLUP_add := BLUP_add[
        ,
        .(TreeID, BLUP_add)
    ][prediction[, .(TreeID)],
        on = .(TreeID = TreeID)
    ][, BLUP_add]]

    # retrieve Mendelian effects
    BLUP_mend <- sln[Model_Term == "TreeID2", .(Level, Effect)]
    setnames(BLUP_mend, c("BLUP_add", "BLUP_mend"))
    # lookup BLUP_mend for prediction
    prediction[, BLUP_mend := BLUP_mend[
        ,
        .(TreeID, BLUP_mend)
    ][prediction[, .(TreeID)],
        on = .(TreeID = TreeID)
    ][, BLUP_mend]]

    # retrieve sca effects
    if (!trial[, "OP" %in% Xtype]) { # trial with only CP families
        BLUP_sca <- sln[Model_Term == "cross", .(Effect)]
    } else { # trial with OP & CP families
        BLUP_sca <- sln[Model_Term == "at(Xtype,2).cross", .(Effect)]
    }
    BLUP_sca[, cross := F@Dimnames[[1]]]
    setnames(BLUP_sca, 1, c("BLUP_sca"))
    # lookup BLUP_sca for prediction
    prediction[, BLUP_sca := BLUP_sca[
        ,
        .(cross, BLUP_sca)
    ][prediction[, .(cross)],
        on = .(cross = cross)
    ][, BLUP_sca]]

    # computing dominance effects (BLUP_sca + BLUP_mend)
    prediction[, BLUP_dom := ifelse(!is.na(BLUP_sca),
        BLUP_mend + BLUP_sca, NA
    )]
    prediction[, BLUP_mend := NULL]
    # remove some columns
    col.remove <- c("mum", "dad")
    prediction[, c(col.remove) := NULL]

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
        cat("Job done...! \n")
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
        cat("Job done...! \n")
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
