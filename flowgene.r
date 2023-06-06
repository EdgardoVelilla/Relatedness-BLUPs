#' flowgene function
#'                                                                                                                                                                                             
#' Author: Edgardo Velilla P.                                                                                                                                                                  
#' email{edgardo.velilla@cmpc.cl}                                                                                                                                                              
#' Created:  21-Jan-2021                                                                                                                                                                       
#' Modified: 08-Jun-2021                                                                                                                                                                       
#' License: GPLv3                                                                                                                                                                              
#' References:                                                                                                                                                                                 
#'                                                                                                                                                                                             
#' Gansner, E.R., Koutsofios, E. and North, S. (2015) Dot User’s Manual.
#' https://www.graphviz.org/pdf/dotguide.pdf
#'                                                                                                                                                                                             
#' Richard Iannone (2020). DiagrammeR: Graph/Network Visualization. R package
#' version 1.0.6.1. https://CRAN.R-project.org/package=DiagrammeR                                                                                                                                               
#'                                                                                                                                                                                             
#' General description:
#'
#' Given a pedigree, its relationship structure are returned as a dot file 
#' suitable for plotting.                                                                                              
#'                                                                                                                                                                                             
#' The resulting dot file can be edited prior to rendering with Graphviz 
#' (https://www.rdocumentation.org/packages/Rgraphviz/versions/2.16.0) or 
#' DiagrammeR packages, or conversion to a graphics file format with the dot 
#' application (see https://graphviz.org/).                                                                                                                  
#'
#'
#' Arguments                                                                                                                                                                                   
#'                                                                                                                                                                                             
#' @param  pedigree                                                                                                                                                                            
#'                                                                                                                                                                                             
#' A data.table/dataframe/matrix where the first three columns correspond to 
#' the identifiers for the individual, mother and father, respectively. The 
#' row giving the pedigree of an individual must appears before any row where
#' that individual appears as parent. Founders use 0 (zero) in the parental 
#' columns.                                                              
#'                                                                                                                                                                                             
#' @param keep                                                                                                                                                                                 
#'                                                                                                                                                                                             
#' A logical vector identifying the rows of pedigree to retain. The default is
#' rep(TRUE, nrow(pedigree)).                                                                                      
#'                                                                                                                                                                                             
#' @param filename                                                                                                                                                                             
#'                                                                                                                                                                                             
#' A character indicating the primary name corresponding to the output file 
#' (dotfile); default filename is "ped" and the suffix ".dot" is appended to 
#' the file name (file´s extension).        
#'                                                                                                                                                                                             
#' @param form                                                                                                                                                                                 
#'                                                                                                                                                                                             
#' Logical indicating if the node´s shape should be returned according to 
#' generation (i.e., rectangle for founders, circle for one unknown parent, 
#' ellipse for progeny). The default is form=FALSE (i.e., all node´s shape
#' are represented as ellipses).                                                                                                                            
#'                                                                                                                                                                                             
#' @param numeric                                                                                                                                                                              
#'                                                                                                                                                                                             
#' Convert pedigree to integer values for graphing. The default is FALSE. That
#' is, mother and father labels are used in the graph.                                                             
#'                                                                                                                                                                                             
#' @param cols                                                                                                                                                                                 
#'                                                                                                                                                                                             
#' A character vector which defines fill colours of nodes for each generation, 
#' starting from founders up to the last generation. The default vector is 
#' c("chocolate4", "olivedrab", "orange2", "gold", "khaki1", "aquamarine", 
#' "bisque1") for the first seven generations.                                                   
#'                                                                                                                                                                                             
#' @param url                                                                                                                                                                                  
#'                                                                                                                                                                                             
#' If not "" (the default), then url is linked to the resulting graph (see Dot
#' User’s Manual for details).                                                                                     
#'                                                                                                                                                                                             
#' @param height                                                                                                                                                                               
#'                                                                                                                                                                                             
#' Node height (see Dot User’s Manual for details).                                                                                                                                            
#'                                                                                                                                                                                             
#' @param width                                                                                                                                                                                
#'                                                                                                                                                                                             
#' Node width (see Dot User’s Manual for details).                                                                                                                                             
#'                                                                                                                                                                                             
#' @param rotate                                                                                                                                                                               
#'                                                                                                                                                                                             
#' If rotate=90 landscape mode is selected; the default is 0.                                                                                                                                  
#'                                                                                                                                                                                             
#'                                                                                                                                                                                             
#' @example                                                                                                                                                                                    
#' Example 12.1 (pag 206) from "Linear Models for the Prediction of Animal 
#' Breeding Values" by Raphael A. Mrode, 3rd Edition (2014)                                                            
#'                                                                                                                                                                                             
#' # pedigre information                                                                                                                                                                       
#' s <- c(0,0,0,0,1,3,6,0,3,3,6,6) # sire                                                                                                                                                      
#' d <- c(0,0,0,0,2,4,5,5,8,8,8,8) # dam                                                                                                                                                       
#' ID <- c(seq(1:12)) # animals                                                                                                                                                                
#' library(data.table)                                                                                                                                                                         
#' pedm9 <- data.table(ID=ID, dam=d, sire=s)                                                                                                                                                   
#' flowgene(pedm9, filename="pedMrode9", form=FALSE, numeric=FALSE, ratio=0.7)                                                                                                                 
#' library(DiagrammeR)                                                                                                                                                                         
#' grViz("pedMrode9.dot")

flowgene <- function(pedigree,
                     keep = rep(TRUE, nrow(pedigree)),
                     filename = "",
                     form = FALSE,
                     numeric = FALSE,
                     cols = c(
                       "chocolate4",
                       "olivedrab",
                       "orange2",
                       "gold",
                       "khaki1",
                       "aquamarine",
                       "bisque1"
                     ),
                     url = "",
                     height,
                     width,
                     rotate,
                     size,
                     ratio) {
  if (filename == "") {
    warning(substitute(filename), "ped.dot is the default name of dotfile")
  }
  if (missing(width)) {
    width <- 0.75
  }
  if (missing(height)) {
    height <- 0.5
  }
  if (missing(rotate)) {
    rotate <- 0
  }
  if (missing(size)) {
    x <- 7.5
    y <- 10
    size <- paste(x, y, sep = ",")
  }
  if (missing(ratio)) ratio <- 1
  if (is.data.frame(pedigree)) {
    pedigree <- as.data.table(pedigree)
  } else if (is.matrix(pedigree)) {
    pedigree <- as.data.table(pedigree)
  } else if (!is.data.table(pedigree)) {
    stop("nothing to do...")
  }
  pedigree <- pedigree[, c(1L:3L)]
  setnames(pedigree, c("TreeID", "mum", "dad"))
  pedigree <- pedigree[keep, ]
  ssize <- dim(pedigree)[1]
  if (numeric) {
    me <- 1:ssize
    mum <- match(pedigree[, mum], pedigree[, TreeID], nomatch = NA)
    dad <- match(pedigree[, dad], pedigree[, TreeID], nomatch = NA)
  } else {
    me <- pedigree[, TreeID]
    mum <- pedigree[, mum]
    mum[mum == 0] <- NA
    dad <- pedigree[, dad]
    dad[dad == 0] <- NA
  }
  generOld <- gener <- rep(ssize + 100, ssize)
  gener[is.na(mum) & is.na(dad)] <- 0
  i <- 0
  while (!all(generOld == gener)) {
    generOld <- gener
    gener[mum %in% me[gener == i]] <- i + 1
    gener[dad %in% me[gener == i]] <- i + 1
    i <- i + 1
    if (i > ssize + 10) break
  }
  color <- rep(NA, nrow(pedigree))
  color[generOld == 0] <- cols[1]
  color[generOld == 1] <- cols[2]
  color[generOld == 2] <- cols[3]
  color[generOld == 3] <- cols[4]
  color[generOld == 4] <- cols[5]
  color[generOld == 5] <- cols[6]
  color[generOld == 6] <- cols[7]
  if (form) { # form: rectangle for founders, circle for one unknown parent
              # and elipse for progeny
    shape <- rep("ellipse", nrow(pedigree)) 
    founders <- as.numeric(is.na(mum)) + as.numeric(is.na(dad))
    shape[generOld == 0] <- "box"
    shape[founders == 1] <- "circle"
  } else {
    shape <- rep("ellipse", nrow(pedigree))
  }
  style <- rep("filled", nrow(pedigree))
  dotfile <- ifelse(filename == "", "ped", filename)
  cat(paste("[", dotfile, "]", sep = ""))
  sink(paste(dotfile, ".dot", sep = ""))
  cat(paste("digraph ped_", dotfile, sep = ""), "{\n")
  cat("fontname=Helvetica;\n")
  cat(paste("rotate=", rotate, sep = ""), " ;\n")
  cat(paste("size=", "\"", size, "\"", sep = ""), " ;\n")
  cat(paste("ratio=", ratio, sep = ""), " ;\n")
  if (url != "") {
    cat(paste("URL=\"", url, "\"", sep = ""), " ;\n")
  }
  cat(paste("\"", me, "\" [style=", style,
    ", color=", color,
    ", shape=", shape,
    ", height=", height,
    ", width=", width,
    ", fontname = Helvetica",
    "] ;",
    sep = ""
  ), sep = "\n")
  cat(paste("\"", dad, "\"", " -> ", "\"", me, "\"", ";",
    sep = ""
  )[!is.na(dad)], sep = "\n")
  cat(paste("\"", mum, "\"", " -> ", "\"", me, "\"", ";",
    sep = ""
  )[!is.na(mum)], sep = "\n")
  cat("}\n")
  sink()
  cat("\n")
  invisible()
}
