#' makeRelatedness function                                                                                                                                                                    
#'                                                                                                                                                                                             
#' Author: Edgardo Velilla P.                                                                                                                                                                  
#' email{edgardo.velilla@cmpc.cl}                                                                                                                                                              
#' Created: 31-Ago-2022                                                                                                                                                                        
#' License: GPLv3                                                                                                                                                                              
#'                                                                                                                                                                                             
#' General description                                                                                                                                                                         
#'                                                                                                                                                                                             
#' This function generates the five relatedness matrices defined under general
#' inbreeding conditions (Harris 1964; Jacquard 1974). That is, the additive 
#' genetic relatedness matrix (A), the dominance relatedness matrix (Dr) in the
#' outbred and homozygous population (Di), the covariance relatedness matrix 
#' between additive and dominance (ADi) in the corresponding homozygous 
#' population, and the relatedness matrix due to sum of squared inbreeding 
#' depression (U) based on Jacquard’s nine condensed coefficients of identity
#' for all pairwise combinations between individuals i and j (Jacquard 1974).                                                                                                                                                        
#'                                                                                                                                                                                             
#' Argument                                                                                                                                                                                    
#'                                                                                                                                                                                             
#' @param  pedigree                                                                                                                                                                            
#'                                                                                                                                                                                             
#' A data.table/dataframe/matrix where the first three columns correspond to 
#' the identifiers for the individual, mother and father, respectively. The 
#' row giving the pedigree of an individual must appears before any row where
#' that individual appears as parent. Founders use 0 (zero) in the parental 
#' columns.                                                                         
#'                                                                                                                                                                                             
#' @return                                                                                                                                                                                     
#'                                                                                                                                                                                             
#' A list with the following matrices in sparse matrix format from a pedigree
#' frame.: A, Dr, Di, ADi, and ID                                                                                   
#'                                                                                                                                                                                             
#' References:                                                                                                                                                                                 
#'                                                                                                                                                                                             
#' Harris D.L. (1964). Genotypic covariances between inbred relatives. Genetics
#' 50:1319–1348.                                                                                                  
#'                                                                                                                                                                                             
#' Jacquard A. (1974) The Genetic Structure of Populations. Springer- Verlag, 
#' New York.                                                                                                        
                                                                                                                                                                                             
makeRelatedness <- function(pedigree,
                            init.ghost = NULL,
                            RAM = 4000L) {

    library(Matrix, quietly = TRUE)
    source("Idcoefs.r")

    if (!is.data.table(pedigree)) {
        if (is.data.frame(pedigree) || is.matrix(pedigree)) {
            pedigree <- as.data.table(pedigree)
        } else {
            stop("nothing to do...")
        }
    }

    # by default, RAM is 4 GB
    if (is.null(RAM)) RAM <- 4000L

    # get only the first three pedigree columns (TreeID, mum, dad)
    # and rename them
    ped3 <- pedigree[, c(1:3)]
    setnames(
      ped3,
      c("TreeID", "mum", "dad")
  )
    # check if there are unknown fathers (i.e., dad = 0). If so, assign them 
    # a new TreeID and add these phantom parents (i.e., TreeID, mum = 0, 
    # dad = 0) to the pedigree
    if (any(ped3[["mum"]] != 0L & ped3[["dad"]] == 0L)) {
        if (is.null(init.ghost)) {
            init.ghost <- 1e6L
        }
        ped3[
            mum != 0L & dad == 0L,
            dad := seq(init.ghost,
                length.out = .N
            )
        ]
        phantom.parents <- data.table(
            TreeID = seq(
                init.ghost,
                max(ped3[["dad"]])
            ),
            mum = 0L,
            dad = 0L
        )
        ped3 <- rbind(phantom.parents, ped3)
    }
    # idcoefs() function requires pedigree's columns must be integers
    cols <- c("TreeID", "mum", "dad")
    ped3[, (cols) := lapply(
        .SD,
        as.integer
    ), .SDcols = cols]

    # computing Jacquard’s condensed coefficients of identity
    # by using function idcoefs()

    output <- idcoefs(ped3,
        RAM,
        verbose = FALSE
    )
    # discard crosses with ghost fathers
    output <- output[x < 1e6, ]
    # build vectors index of unique values of x and y
    ux <- as.vector(unique(output[, x]))
    idx <- seq(1:length(ux))
    index.x <- data.table(
        x = ux,
        idx = idx
    )
    output <- index.x[output,
        on = .(x = x)
    ]

    uy <- as.vector(unique(output[, y]))
    idy <- seq(1:length(uy))
    index.y <- data.table(
        y = uy,
        idy = idy
    )
    output <- index.y[output,
        on = .(y = y)
    ]
    
    # convert vectors ux and uy to integers for lookups
    cols <- c("idx", "idy")
    output[, (cols) := lapply(.SD, as.double),
        .SDcols = cols
    ]

    m <- max(output[, idx])

    # bulding the numerator relationship matrix from
    # condensed coefficients of identity (Jacquard 1974)
    ksp <- 2 * (with(
        output,
        Matrix::sparseMatrix(
            i = idx,
            j = idy, x = J1 + 0.5 * (J3 + J5 + J7)
                + 0.25 * J8,
            dims = c(m, m),
            dimnames = list(
                pedigree[[1]],
                pedigree[[1]]
            ),
            triangular = FALSE,
            check = TRUE
        )
    ))
    As <- forceSymmetric(new(
        "dgCMatrix",
        i = ksp@i,
        p = ksp@p,
        Dim = ksp@Dim,
        x = ksp@x,
        Dimnames = list(
            ksp@Dimnames[[1]],
            ksp@Dimnames[[2]]
        )
    ))
    As <- drop0(As,
        tol = 1e-15,
        is.Csparse = NA
    )

    # bulding the dominance relationship matrix (Dr) in the outbred population 
    # J7 represents the expected fraternity coefficient in terms of condensed
    # identity coefficients between two non-inbred individuals
    Dr <- (with(
        output,
        Matrix::sparseMatrix(
            i = idx,
            j = idy, x = J7,
            dims = c(m, m),
            dimnames = list(
                pedigree[[1]],
                pedigree[[1]]
            ),
            triangular = FALSE,
            check = TRUE
        )
    ))
    Drs <- forceSymmetric(new(
        "dgCMatrix",
        i = Dr@i, p = Dr@p,
        Dim = Dr@Dim,
        x = Dr@x,
        Dimnames = list(
            Dr@Dimnames[[1]],
            Dr@Dimnames[[2]]
        )
    ))
    Drs <- drop0(Drs,
        tol = 1e-15,
        is.Csparse = NA
    )

    # bulding the dominance relationship matrix
    # in the completely inbred population (DI)
    Di <- (with(
        output,
        Matrix::sparseMatrix(
            i = idx,
            j = idy,
            x = J1,
            dims = c(m, m),
            dimnames = list(
                pedigree[[1]],
                pedigree[[1]]
            ),
            triangular = FALSE,
            check = TRUE
        )
    ))
    Dis <- forceSymmetric(
        new("dgCMatrix",
            i = Di@i,
            p = Di@p,
            Dim = Di@Dim,
            x = Di@x,
            Dimnames = list(
                Di@Dimnames[[1]],
                Di@Dimnames[[2]]
            )
        )
    )
    Dis <- drop0(Dis,
        tol = 1e-15,
        is.Csparse = NA
    )

    # bulding the covariance relationship matrix
    # between additive and dominance effects in
    # the homozygous population (ADI)
    ADi <- (with(
        output,
        Matrix::sparseMatrix(
            i = idx,
            j = idy,
            x = 4 * J1 + J3 + J5,
            dims = c(m, m),
            dimnames = list(
                pedigree[[1]],
                pedigree[[1]]
            ),
            triangular = FALSE,
            check = TRUE
        )
    ))
    ADis <- forceSymmetric(
        new("dgCMatrix",
            i = ADi@i,
            p = ADi@p,
            Dim = ADi@Dim,
            x = ADi@x,
            Dimnames = list(
                ADi@Dimnames[[1]],
                ADi@Dimnames[[2]]
            )
        )
    )
    ADis <- drop0(ADis,
        tol = 1e-15,
        is.Csparse = NA
    )

    # calculate the inbreeding coefficients
    # for each individual (fx, fy)
    output[, fx := J1 + J2 + J3 + J4]
    output[, fy := J1 + J2 + J5 + J6]

    # bulding the relationship matrix
    # for the square of the inbreeding depression (ID)
    ID <- (with(
        output,
        Matrix::sparseMatrix(
            i = idx,
            j = idy,
            x = J1 + J2 - fx * fy,
            dims = c(m, m),
            dimnames = list(
                pedigree[[1]],
                pedigree[[1]]
            ),
            triangular = FALSE,
            check = TRUE
        )
    ))
    IDs <- forceSymmetric(
        new("dgCMatrix",
            i = ID@i,
            p = ID@p,
            Dim = ID@Dim,
            x = ID@x,
            Dimnames = list(
                ID@Dimnames[[1]],
                ID@Dimnames[[2]]
            )
        )
    )
    IDs <- drop0(IDs,
        tol = 1e-15,
        is.Csparse = NA
    )

    return(list(
        A = As,
        Dr = Drs,
        Di = Dis,
        ADi = ADis,
        ID = IDs
    ))
}
