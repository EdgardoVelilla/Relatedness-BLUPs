
# Convert a full pedigree into a family pedigree keeping key links between families across generations

pedFam <- function(pedigree) {

  library(data.table, quietly = TRUE)
  source("justCheck.r")

  if (!is.data.table(pedigree)) {
    if (is.data.frame(pedigree) || is.matrix(pedigree)) {
      pedigree <- as.data.table(pedigree)
    } else {
      stop("nothing to do...")
    }
  }
  ped3 <- pedigree[, c(1:3)]
  setnames(
      ped3,
      c("TreeID", "mum", "dad")
  )
  ped <- makeFam(ped3)[!is.na(cross), ]
  setkey(ped, cross)
  ped <- ped[, .SD[1L], by = key(ped)]
  id.keep <- ped[, TreeID]

  # function to prune the pedigree but retaining the key ancestors
  cutPed <- function(ped3, id.keep) {
    ped.tmp <- copy(pedigree)
    n.id <- length(id.keep) + 1L
    while (length(id.keep) != n.id) {
      n.id <- length(id.keep)
      id.keep <- union(na.omit(c(unlist(ped.tmp[
        ,
        .(mum, dad)
      ][match(id.keep, ped.tmp[, TreeID]), ]))), id.keep)
    }
    ped.tmp <- ped.tmp[sort(match(
      id.keep,
      ped.tmp[, TreeID]
    )), ]
    # add the generation ("gen") field
    ped.tmp[, gen := gen.add(ped.tmp)[, gen]]
    ped.tmp[]
  }

  pedfam0 <- cutPed(ped3, id.keep)
  setorder(pedfam0, gen)
  pedfam0[, gen := NULL]

  # check if some mothers or fathers are missing
  justCheck(pedfam0)
  pedfam0[]
}
