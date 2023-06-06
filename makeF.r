#' makeF() function
#'
#' Author: Edgardo Velilla P.
#' email{edgardo.velilla@cmpc.cl}
#' Created: 15-March-2022
#' Modified: 22-April-2023
#' License: GPLv3
#'
#' General description:
#'
#' This function builds the relatedness matrix among full-sib families (F)
#' based on the expected fraternity coefficients obtained from Jacquard’s
#' nine condensed coefficients of identity (Jacquard 1974). These coefficients
#' can be computed by using the recursive algorithm proposed by Karigl (1981)
#' and implemented by Abney's (2009) graphical algorithm. However, the critical
#' assumption of this approach is that both parents of each individual must be
#' known. Consequently, this procedure cannot be directly applied to pedigrees
#' that include open-pollinated families. To overcome this issue, an
#' alternative pedigree version is defined ("ped.jacquard"). That is, unknown
#' fathers (denoted by zeros) will be replaced with a big integer number (e.g.,
#' starting from 1e6 to avoid potential coding conflicts) to identify each
#' individual's father. Furthermore, these unknown fathers will be added as
#' phantom parents at the beginning of this updated pedigree.
#'
#' Usage:
#' makeF(pedigree, init.ghost = NULL, gen = NULL, RAM = 4000, eig.tol = NULL)
#'
#' Arguments:
#'
#' @param pedigree
#'
#' A data.table/data.frame/matrix with at least three columns: TreeID, mum, and
#' dad. The first column (TreeID) is a unique identifier for each individual.
#' The second column (mum) is the mother's TreeID. The third column (dad) is
#' the father's TreeID. If the father is unknown, the value of this column
#' must be zero.
#'
#' @param init.ghost
#'
#' Integer indicating the initial TreeID code for unknown (ghost) parents for
#' open-pollinated's selection in the first-generation. The default is
#' @init.ghost = NULL.
#'
#' @param gen: generation
#'
#' Integer value corresponding to the last breeding program's generation. The
#' default is @gen = NULL.
#'
#' @param RAM
#'
#' Specify the amount of RAM (MB) the executable program uses (“idcoefs.exe”).
#' The default is @RAM = 4000 (4 GB), but a higher allocation is recommended.
#'
#' @param eig.tol
#'
#' set the limits or tolerance to define which eigenvalues of the F matrix are
#' "numerically" zero.
#'
#'
#' @return
#'
#' A list with the following components:
#'
#' F
#'
#' The relateddness matrix between full-sib families effect in sparse matrix
#' format.
#'
#' ped.jacquard
#'
#' a suitable pedigree that can be reused to compute Jacquard coefficients via
#'  the idcoef function.
#'
#' pedfam
#'
#' a family version of the pedigree "ped.jacquard". That is, only one
#' individual represents each family.
#'
#' cond
#'
#' the condition number of the F matrix using the “spectral norm”. A finite
#' large condition number means that the matrix is close to being singular
#' (Anderson et al. 1994).
#'
#' status
#'
#' the resulting condition (positive-definite, positive semi-definite,
#' indefinite, and so on) of the F matrix based on their eigenvalues and the
#' tolerance ("eig.tol") defined.
#'
#'
#' References:
#'
#' Anderson. E. and ten others (1999) LAPACK Users' Guide. Third Edition. SIAM.
#' Available on-line at https://www.netlib.org/lapack/lug/lapack_lug.html.
#'
#' Abney M. (2009) A graphical algorithm for fast computation of identity
#' coefficients and generalized kinship coefficients. Bioinformatics.
#' 25(12):1561-3.
#'
#' Jacquard A. (1974) The Genetic Structure of Populations. Springer-
#' Verlag, New York.
#'
#' Karigl G. (1981) A recursive algorithm for the calculation of
#' identity coefficients. Ann Hum Genet 45:299–305.
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
#' output <- makeF(pedigree)
#'
#' > names(output)
#' [1] "F"  "ped.jacquard" "pedfam" "cond" "status"

makeF <- function(pedigree,
                  init.ghost = NULL,
                  RAM = 4000L, # default value is 4 GB
                  eig.tol = NULL) { # specify the tolerance to define which  
                                    # eigenvalues are "numerically" zero  
                    
    library(Matrix, quietly = TRUE)
    library(data.table, quietly = TRUE)

    source("matrixMethods.r")
    source("pedFam.r")
    source("Idcoefs.r")
    source("condNumber.r")

    # check if pedigree is a data.table
    # if not, convert it
    if (!is.data.table(pedigree)) {
        if (is.data.frame(pedigree) || is.matrix(pedigree)) {
            pedigree <- as.data.table(pedigree)
        } else {
            stop("nothing to do...")
        }
    }

    # get only the first three pedigree columns (TreeID, mum, dad)
    # and rename them
    ped3 <- pedigree[, c(1:3)]
    setnames(ped3, c("TreeID", "mum", "dad"))

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

    # convert the full pedigree into a family pedigree with function pedFam()
    pedfam <- pedFam(ped3)
    # idcoefs() function requires pedigree's columns must be integers
    cols <- c("TreeID", "mum", "dad")
    pedfam[, (cols) := lapply(
        .SD,
        as.integer
    ), .SDcols = cols]

    # computing Jacquard’s condensed coefficients by using function idcoefs()
    output <- idcoefs(pedfam,
        RAM,
        verbose = FALSE
    )

    # identity labels ("cross") for each full-sib family
    crosses <- makeFam(pedfam)[!is.na(cross), ]
    setkey(crosses, cross)

    if (!is.null(init.ghost)) {
        # identify which tree represents its own family
        crosses <- crosses[, .SD[1L],
            by = key(crosses)
        ][dad < init.ghost, ]
        # discard identity coefficients from ghost fathers
        output <- output[x < init.ghost & y < init.ghost, ]
    } else {
        # discard crosses with ghost fathers
        crosses <- crosses[dad < 1e6, ]
        crosses <- crosses[, .SD[1L],
            by = key(crosses)
        ]
    }

    # Convert the columns x and y to integers to be able to perform lookups
    cols2 <- c("x", "y")
    output[, (cols2) := lapply(
        .SD,
        as.integer
    ), .SDcols = cols2]
    # update output to match the dimension of the F matrix
    out2 <- output[x %in% crosses[, TreeID], ]
    output.update <- out2[y %in% crosses[, TreeID], ]

    # build vectors of index values for x and y
    ux <- as.vector(unique(output.update[, x]))
    idx <- seq(1L:length(ux))
    index.x <- data.table(
        x = ux,
        idx = idx
    )
    output.update <- index.x[output.update,
        on = .(x = x)
    ]

    uy <- as.vector(unique(output.update[, y]))
    idy <- seq(1L:length(uy))
    index.y <- data.table(
        y = uy,
        idy = idy
    )
    output.update <- index.y[output.update,
        on = .(y = y)
    ]
    # convert vectors ux and uy to integers for lookups
    cols3 <- c("idx", "idy")
    output.update[, (cols3) := lapply(.SD, as.integer),
        .SDcols = cols3
    ]
    # J7 represents the expected fraternity coefficients 
    # from condensed identity coefficients (Jacquard 1974)
    output.update[, J7F := ifelse((x == y), J7, J7 * 4)]

    # label rows/cols of the F matrix
    cross.label <- as.vector(crosses[[1L]])
    # set the dimension of F
    m <- length(cross.label)

    # build a sparse relatedness matrix (F) from J7F
    F <- (with(
        output.update,
        Matrix::sparseMatrix(
            i = idx,
            j = idy, x = J7F,
            dims = c(m, m),
            dimnames = list(
                cross.label,
                cross.label
            ),
            triangular = FALSE,
            check = TRUE
        )
    ))

    # the F matrix is given in sparse format, i.e., it has three columns per
    # line (row, column, value; "dgC" Matrix class). Hence, it is converted
    # to a symmetric matrix
    Fs <- forceSymmetric(new(
        "dgCMatrix",
        i = F@i,
        p = F@p,
        Dim = F@Dim,
        x = F@x,
        Dimnames = list(
            F@Dimnames[[1L]],
            F@Dimnames[[2L]]
        )
    ))
    Fs <- drop0(Fs,
        tol = 1e-15,
        is.Csparse = NA
    )

    # compute and check the condition number of the F matrix
    cnd <- condNumber(Fs, norm = "2")
    if (cnd > 4.5e15) {
        warning("Matrix 'F' is ill-conditioned")
    }
    # setting limits to define which eigenvalues are "numerically" zero
    if (!is.null(eig.tol)) {
        if ((!is.numeric(eig.tol) || length(eig.tol) != 1 || eig.tol < 0)) {
            stop("Argument 'eig.tol' must be a non-negative numeric scalar")
        }
    } else {
        # compute the eigenvalues of the F matrix
        eigen <- eigen(Fs, symmetric = TRUE)[["values"]]
        eig.tol <- length(eigen) * eps(eigen) # for details, see eps() function
    }
    # generate a mask with positive relative eigenvalues
    eigen[abs(eigen) < eig.tol] <- 0

    # examing whether the matrix F is positive definite, negative definite,
    # and so on
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
    if (any(abs(eigen) < eig.tol)) {
        warning("Matrix 'F' is singular or nearly singular")
    }

    return(list(
        F = Fs,
        pedfam = pedfam,
        cond = cnd,
        status = status
    ))
}
