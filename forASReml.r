#' forASReml function
#'
#' Author: Edgardo Velilla P.
#' email{edgardo.velilla@cmpc.cl}
#' Created: 12-Oct-2020
#' Modified (simplified!): 02-May-2023
#'
#' References:
#'
#' Butler, D. G., Cullis, B.R., A. R. Gilmour, Gogel, B.G. and Thompson, R.
#' 2017. ASReml-R Reference Manual Version 4. VSN International Ltd, Hemel
#' Hempstead, HP1 1ES, UK.
#'
#' Gilmour, A. R., Gogel, B. J., Cullis, B. R., Welham, S. J. and Thompson,
#' R. (2021). ASReml User Guide Release 4.2 Functional Specification, VSN
#' International Ltd, Hemel Hempstead, HP2 4TP, UK, www.vsni.co.uk
#'
#' General description:
#'
#' Converts a matrix object into a suitable format to be used by ASReml
#' (Gilmour et al. 2021) or ASReml-R (Butler at al. 2017) software.
#'
#' Arguments
#'
#' @param G
#'
#' An object from "matrix" or "Matrix" class representing a relatedness matrix
#' or its inverse in sparse format (three column per line).
#'
#'
#' @param Ginv
#'
#' Logical indicating if the relatedness matrix represents its inverse; the
#' default is Ginv = TRUE.
#'
#' @param rowNames
#'
#' A character vector set as the "rowNames" attribute of the G sparse matrix.
#' The default is the dimnames(G)[[1].
#'
#' @param colnames
#'
#' A character vector with the names of the columns of the G sparse matrix. The
#' default is c("row", "column", "Ginv").
#'
#'
#' @return
#'
#' The inverse of the G matrix as a lower triangular matrix in sparse format 
#' (three numbers per line, i.e., row, column, value) in row-major order 
#' (columns sorted within rows) with attribute "INVERSE" set to TRUE and 
#' attribute rowNames (see details in Gilmour et al. 2021, p. 166 - 169).
#' 
#' @example  
#'                                                                                                                                                                                   
#' Example 12.1 (pag 206) from "Linear Models for the Prediction of Animal 
#' Breeding Values" by Raphael A. Mrode, 3rd Edition (2014).      
#'
#' library(data.table)
#' mum <- c(0,0,0,0,2,4,5,5,8,8,8,8) # mothers
#' dad <- c(0,0,0,0,1,3,6,0,3,3,6,6) # fathers
#' ID <- c(seq(1:length(dad))) # individuals
#' pedigree <- data.table(ID, mum, dad)
#' F <- makeF(pedigree)[["F"]]
#' Finv.as <- forASReml(F, Ginv = FALSE)

forASReml <- function(G,
                      Ginv = TRUE,
                      rowNames = NULL,
                      colnames = c("row", "column", "Ginv")) {

    library(Matrix, quietly = TRUE) # version 1.5-3
    library(data.table, quietly = TRUE) # version 1.14.6
    source("matrixMethods.r")
    if (!inherits(G, "dgTMatrix")) {
        if (inherits(G, "dgCMatrix") || inherits(G, "dsCMatrix") ||
            inherits(G, "dgRMatrix")) {
            G <- as(G, "dgTMatrix")
        } else if (inherits(G, "matrix")) {
            G <- as(as(G, "sparseMatrix"), "dgTMatrix")
        } else {
            stop(
                substitute(G),
                " must be an object from 'matrix' or 'Matrix' class in sparse format"
            )
        }
    }
    if (Ginv) { # G is an inverse relatedness matrix (.giv in ASReml parlance)
        ginv <- data.table(
            row = G@i + 1,
            col = G@j + 1,
            Ginv = G@x
        )
        # ASReml requires a lower triangular matrix in row-major order
        # (i.e., row >= col)
        setorder(ginv, row, col)
        ginv <- as.matrix(ginv[row >= col])
        # attribute INVERSE is used in ASReml-R
        # to check if the matrix is an inverse
        attr(ginv, "INVERSE") <- TRUE
        colnames(ginv) <- colnames
        if (is.null(rowNames)) {
            if (length(G@Dimnames[[1]]) != 0) {
                rowNames <- as.character(G@Dimnames[[1]])
                attr(ginv, "rowNames") <- rowNames
            } else {
                warning(
                    substitute(G),
                    "  has no 'rowNames' attribute"
                )
            }
        } else {
            attr(ginv, "rowNames") <- as.character(rowNames)
        }
        ginv[]
    } else { # G is a relatedness matrix (.grm in ASReml parlance)
        # create the inverse of G by using function solve.CHMperm()
        ginv <- solve.CHMperm(G)
        ginv <- as(ginv, "dgTMatrix")
        ginv <- data.table(
            row = ginv@i + 1,
            col = ginv@j + 1,
            Ginv = ginv@x
        )
        setorder(ginv, row, col)
        ginv <- as.matrix(ginv[row >= col])
        attr(ginv, "INVERSE") <- TRUE
        colnames(ginv) <- colnames
        if (is.null(rowNames)) {
            if (length(G@Dimnames[[1]]) != 0) {
                rowNames <- as.character(G@Dimnames[[1]])
                attr(ginv, "rowNames") <- rowNames
            } else {
                warning(
                    substitute(G),
                    "  has no 'rowNames' attribute"
                )
            }
        } else {
            attr(ginv, "rowNames") <- as.character(rowNames)
        }
    }
    ginv[]
}
