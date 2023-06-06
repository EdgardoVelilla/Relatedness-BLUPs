selection <- function(base.pop,
                      nfounder,
                      p = NULL,
                      sel = c("rnd", "egv", "tbv", "tdv", "pheno"),
                      selectTop = TRUE) {

    if (is.data.frame(base.pop)) {
        base.pop <- as.data.table(base.pop)
    } else if (is.matrix(base.pop)) {
        base.pop <- as.data.table(base.pop)
    } else if (is.data.table(base.pop)) {
        base.pop <- base.pop
    } else {
        stop("nothing to do...")
    }
	
    nmothers <- ceiling(p * nfounder)
    if (missing(sel)) {
        stop(" must be specify the selection method: rnd, egv, tbv, tdv, pheno) ")
    }
    sel <- tolower(sel)
    if (sel == "rnd") {
        base.pop[, response := rep(0, dim(base.pop)[1])]
        parent.g.tmp <- sample(base.pop[, TreeID], 
        size = nmothers, replace = FALSE)
        mean.pheno <- mean(base.pop[parent.g.tmp, pheno])
        mean.geno <- mean(base.pop[parent.g.tmp, egv])
    } else if (sel == "egv") {
        base.pop[, response := egv]
        mother <- base.pop[order(response, 
        decreasing = selectTop), TreeID]
        parent.g.tmp <- mother[1:nmothers]
        mean.pheno <- mean(base.pop[parent.g.tmp, pheno])
        mean.geno <- mean(base.pop[parent.g.tmp, egv])
    } else if (sel == "tbv") {
        base.pop[, response := tbv]
        mother <- base.pop[order(response, 
        decreasing = selectTop), TreeID]
        parent.g.tmp <- mother[1:nmothers]
        mean.pheno <- mean(base.pop[parent.g.tmp, pheno])
        mean.geno <- mean(base.pop[parent.g.tmp, egv])
    } else if (sel == "tdv") {
        base.pop[, response := tdv]
        mother <- base.pop[order(response, 
        decreasing = selectTop), TreeID]
        parent.g.tmp <- mother[1:nmothers]
        mean.pheno <- mean(base.pop[parent.g.tmp, pheno])
        mean.geno <- mean(base.pop[parent.g.tmp, egv])
    } else if (sel == "pheno") {
        base.pop[, response := pheno]
        mother <- base.pop[order(response, 
        decreasing = selectTop), TreeID]
        parent.g.tmp <- mother[1:nmothers]
        mean.pheno <- mean(base.pop[parent.g.tmp, pheno])
        mean.geno <- mean(base.pop[parent.g.tmp, egv])
    } else {
        stop(paste0("sel=", sel, " is not an option"))
    }
    mean.values <- vector(length = 2, mode = "list")
    mean.values[[1]] <- mean.pheno
    mean.values[[2]] <- mean.pheno
    base.pop[, response := NULL]
    return(list(
        parents = parent.g.tmp,
        mean.values = mean.values
    ))
}
