makeBasepop <- function(nfounder,
                        init.ghost,
                        n.halfsib,
                        p = NULL,
                        h2 = NULL,
                        H2 = NULL,
                        mu = NULL) {
    library(data.table, quietly = TRUE)

    if (is.null(p)) p <- 0.05
    if (is.null(h2)) h2 <- 0.2 # narrow-sense heritability
    if (is.null(H2)) H2 <- 0.3 # broad-sense heritability
    if (is.null(mu)) mu <- 4.36018

    sigma2.add <- h2 
    sigma2.dom <- H2 - h2 
    sigma2.sca <- sigma2.dom / 4 
    sigma2.e <- 1 - H2 

    founders <- data.table(
        TreeID = seq(1L:nfounder),
        mum = rep(0L, nfounder),
        dad = rep(0L, nfounder),
        gen = rep(0L, nfounder),
        stringsAsFactors = FALSE
    )
    ###      sampling founder's breeding values     ###
    tbv <- rnorm(nfounder) * sqrt(sigma2.add)  # add ~ N(0, sigma2.add)
    # scaling breeding values
    tbv <- as.vector(scale(tbv) * sqrt(sigma2.add)) 

    ###      sampling founder's dominance deviation     ###
    tdv <- rnorm(nfounder) * sqrt(sigma2.dom) # dom ~ N(0, sigma2.dom)
    # scaling dominance deviations
    tdv <- as.vector(scale(tdv) * sqrt(sigma2.dom)) 

    ## building the base population: true breeding values (tbv) and true 
    # dominance deviations (tdv)     
    data.founder <- data.table(cbind(tbv, tdv), stringsAsFactors = FALSE)
    base.pop <- cbind(founders, data.founder)
    base.pop[, tdv.adj := tdv]

    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    ###   building the first-generation open-pollinated progeny trial       ###
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

    # individual additive deviation from its family average
    within.add <- rnorm(nfounder * n.halfsib) * sqrt(0.75 * sigma2.add) 
    within.add <- as.vector(scale(within.add) * sqrt(0.75 * sigma2.add))
    # individual dominance deviation from its family average
    within.dom <- rnorm(nfounder * n.halfsib) * sqrt(sigma2.dom) 
    within.dom <- as.vector(scale(within.dom) * sqrt(sigma2.dom))
    # individual-specific effects
    resid <- rnorm(nfounder * n.halfsib) * sqrt(sigma2.e) 
    resid <- as.vector(scale(resid) * sqrt(sigma2.e))

    base.gen1 <- data.table(cbind(within.add, within.dom, resid))
    base.gen1[, TreeID := seq(
        max(base.pop[, TreeID]) + 1L,
        max(base.pop[, TreeID]) + nfounder * n.halfsib
    )]

    base.gen1[, mum := rep(1L:nfounder, each = n.halfsib)]
    base.gen1[, dad := rep(rep(0L, times = nfounder), each = n.halfsib)]
    base.gen1 <- base.pop[, .(TreeID, tbv)][base.gen1, on = .(TreeID = mum)]
    setnames(base.gen1, c(1L:2L, 6L), c("mum", "bv_mum", "TreeID"))
    base.gen1[, amongHS := 0.5 * bv_mum]

    # compute true breeding value (tbv)
    base.gen1[, tbv := amongHS + within.add] 
    # compute true dominance value (tdv)
    base.gen1[, tdv := within.dom] 
    base.gen1[, egv := tbv + tdv]
    base.gen1[, pheno := mu + egv + resid]

    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    ###   ranking and selection the best half-sib families for first        ###
    ###   generation open-pollinated families                               ###
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

    # within-family "genetic" effect for half-sib
    base.gen1[, withinHS := within.add + within.dom] 
    # within-family effects
    base.gen1[, pheno.withinHS := withinHS + resid] 
    # mean within-family phenotypic effect by each half-sib
    mean.withinHS <- base.gen1[, mean(pheno.withinHS), by = mum]
    setnames(mean.withinHS, 2, c("mean.wHS"))
    # ranking half-sib families according to half-sib family means (HS.mean)
    HSfamily.means <- base.gen1[, .(mum, amongHS)][mean.withinHS,
        on = .(mum = mum)
    ][, .SD[1L], by = mum]
    HSfamily.means[, HS.mean := amongHS + mean.wHS]
    setorder(HSfamily.means, -HS.mean, mum)
    # number of half-sib families to select
    n.mother <- ceiling(p * nfounder) 
    # selection the best half-sib families
    select.mother <- HSfamily.means[1L:n.mother, .(mum)] 
    # individuals from the progeny trial that belong to families selected
    base.sel <- base.gen1[, .(TreeID, mum, tbv)][select.mother,
        on = .(mum = mum)
    ]

    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    ###   forward selection based on tbv (i.e., selection of the best       ###
    ###   individual from each best family)                                 ###
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

    select <- base.sel[base.sel[, .I[tbv == max(tbv)],
        by = mum
    ]$V1][order(TreeID), ]
    # ghost's fathers to match pedigree for selection from 
    # 1st open-pollination population
    select[, dad := seq.int(init.ghost,
        length.out = dim(select)[1L]
    )]
    # assign father's name to slection from 1st open-pollinated trial
    base.gen1[
        TreeID %in% select[, TreeID],
        dad := select[, dad]
    ]
    ghost <- data.table(
        TreeID = select[, dad],
        mum = rep(0L, dim(select)[1L]),
        dad = rep(0L, dim(select)[1L]),
        gen = rep(0L, dim(select)[1L]),
        stringsAsFactors = FALSE
    )
    # base population including ghost fathers
    base.founder <- rbind(ghost, founders)
    base.founder[, Xtype := rep("OP",
        times = dim(base.founder)[1L]
    )]

    # building additional fields to match base from 1st generation with 2nd 
    # and subsequent generations's trials
    base.gen1[, block := rep(1L:n.halfsib,
        times = nfounder
    )]
    base.gen1[, gen := rep(1L,
        times = nfounder * n.halfsib
    )]
    base.gen1[, f := rep(0,
        times = dim(base.gen1)[1]
    )]
    base.gen1[, Xtype := rep("OP",
        times = nfounder * n.halfsib
    )]
    base.gen1[, cross := as.factor(rep(NA,
        times = nfounder * n.halfsib
    ))]
    base.gen1[, among.dom := rep(0,
        times = nfounder * n.halfsib
    )]
    base.gen1[, among.dom.adj := rep(0,
        times = nfounder * n.halfsib
    )]
    base.gen1[, tdv.adj := tdv]
    base.gen1[, id := rep(0,
        times = nfounder * n.halfsib
    )]
    cols.delete <- c("bv_mum", "pheno.withinHS", "withinHS")
    base.gen1[, (cols.delete) := NULL]
    setnames(base.gen1, 7, c("among.add"))
    setcolorder(base.gen1, c(
        "TreeID", "mum", "dad", "block", "gen", "f", "Xtype", "cross",
        "among.add", "within.add", "tbv", "among.dom", "among.dom.adj",
        "within.dom", "tdv", "id", "tdv.adj", "egv", "resid", "pheno"
    ))
    setorder(base.gen1, block)

    return(list(
        base.gen1 = base.gen1,
        base.pop = base.pop,
        base.founder = base.founder,
        sel = select
    ))
}
  
  

  
  