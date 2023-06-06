#' idcoefs function
#'                                                                                                                                                                                            
#' Author: Edgardo Velilla P.                                                                                                                                                                 
#' email{edgardo.velilla@cmpc.cl}                                                                                                                                                             
#' Created: 31-Ago-2022                                                                                                                                                                       
#' License: GPLv3                                                                                                                                                                             
#'                                                                                                                                                                                            
#' General description                                                                                                                                                                        
#'                                                                                                                                                                                            
#' This function generates Jacquard’s nine condensed coefficients of identity
#' for all pairwise combinations between individuals i and j (Jacquard 1974)
#' based on the recursive algorithm      
#' proposed by Karigl (1981) and implemented by Abney’s (2009) graphical 
#' algorithm in the software IdCoefs (C++ codes), availability in:
#' https://code.google.com/archive/p/idcoefs/downloads  
#' See also https://home.uchicago.edu/abney/abney_web/Identity_Coefficients.html                                                                                                              
#'                                                                                                                                                                                            
#' Note: This function is a wrapper for the original C ++ software Idcoefs 
#' (version 2.1.1) written by Mark Abney (2009). Code was compiled for (64-bit)
#' Windows machine; for this to work, the executable (idcoefs.exe) must 
#' be placed in either the current working directory or a folder included in 
#' the PATH variable.                                                             
#'                                                                                                                                                                                           
#' References:                                                                                                                                                                                
#'                                                                                                                                                                                            
#' Abney M. (2009) A graphical algorithm for fast computation of identity 
#' coefficients and generalized kinship coefficients. Bioinformatics. 
#' 25(12):1561-3.                                   
#'                                                                                                                                                                                            
#' Jacquard A. (1974) The Genetic Structure of Populations. Springer- Verlag,
#' New York.                                                                                                       
#'                                                                                                                                                                                            
#' Karigl G. (1981) A recursive algorithm for the calculation of identity 
#' coefficients. Ann Hum Genet 45:299–305.      
#'                                                                        

idcoefs <- function(pedigree,
                    RAM = 4000,
                    verbose = FALSE) {
                      
  library(data.table, quietly = TRUE)
  if (!is.data.table(pedigree)) {
    if (is.data.frame(pedigree) || is.matrix(pedigree)) {
      pedigree <- as.data.table(pedigree)
    } else {
      stop("nothing to do...")
    }
  }
  if (is.null(RAM)) RAM <- 2000L
  pedigree <- pedigree[, c(1L:3L)]
  setnames(pedigree, c("TreeID", "mum", "dad"))
  cols <- c("TreeID", "mum", "dad")
  pedigree[, (cols) := lapply(
    .SD,
    as.integer
  ), .SDcols = cols]

  fwrite(pedigree,
    file = "pedi.csv", quote = FALSE,
    sep = " ", row.names = FALSE, col.names = FALSE
  )
  fwrite(pedigree[, .(TreeID)],
    file = "sample.csv",
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  command <- paste0(
    "idcoefs -p pedi.csv -s sample.csv -o output.csv -r ",
    paste(RAM)
  )
  res <- suppressWarnings(system(
    command,
    wait = TRUE, intern = TRUE
  ))
  if (verbose) message(res)

  output <- fread("output.csv")
  setnames(output,
    old = colnames(output),
    new = c(
      "x", "y",
      c(paste0("J", seq(1L:9L)))
    )
  )
  output[]
}
