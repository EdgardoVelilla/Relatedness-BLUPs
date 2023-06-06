basicSim <- function(seed,
                     h2 = 0.25,
                     H2 = 0.5,
                     mu = 4.36018, # the overall mean of the trait
                     gen = 4, # number of generations
                     nfam = 2, # number of families per generation
                     nprog = 2, # number of progenies per family
                     OP = TRUE) { # 1st generation is open-pollinated (OP)

    library(data.table, quietly = TRUE) # version 1.14.6
    library(Matrix, quietly = TRUE) # version 1.5-3
    library(nadiv)
    library(pedigree)

    source("matrixMethods.r")
    source("makeRelatedness.r")
    source("rmvChol.r")

    if (is.null(seed)) {
        seed <- 123
    }
    set.seed(seed)
    # simulating the pedigree
    pedigree <- simPedMCN(
        pedTemp = NULL,
        g = gen,
        Nfam = nfam,
        noff = nprog
    )
    pedigree <- numPed(pedigree[, c(1:3)])
    pedigree[pedigree == -998] <- 0
    class(pedigree) <- c(class(pedigree), "matrix")
    pedigree <- as.data.table(pedigree)
    setnames(
        pedigree,
        c("TreeID", "mum", "dad")
    )

    # add generation ("gen") field
    pedigree[, gen := gen.add(pedigree)[, gen]] 

    if (OP) { # 1st generation is OP
        pedigree[gen == 1, dad := 0L]
        # Identify phantom parents
        phantom_parents <- pedigree[mum == 0L & dad == 0L, ][, TreeID]
        # Identify all individuals that appear as either a mother or a father
        parents <- unique(c(
            pedigree[mum != 0, mum],
            pedigree[dad != 0, dad]
        ))
        # Check if phantom parents appear in parents vector
        if (!all(phantom_parents %in% parents)) {
            # remove the pedigree's rows where phantom parents don't appear
            pedigree <- pedigree[!(mum == 0L & dad == 0L) |
                TreeID %in% parents]
        }
    }

    # add inbreeding coefficient ("f") field
    pedigree[, f := calcInbreeding(pedigree)]
    # build the trial from the pedigree
    trial <- pedigree[gen > 0, ]
    # adding Xtype field to trial
    trial[, Xtype := ifelse((
        mum != 0L & dad != 0L),
    "CP", "OP"
    )]
    # adding cross field
    trial[, cross := makeFam(trial)[, cross]]
    # adding family field
    trial[, family := ifelse((mum != 0L & dad != 0L),
        cross, # for full-sib families family = cross
        mum # for half-sib families family = mum
    )]

    n <- dim(trial)[1] # number of individuals in the trial
    m <- dim(pedigree)[1] # number of individuals in the pedigree
    # mu <- 4.36018

    # setting sigmas
    sigma2.add <- h2
    sigma2.dom <- H2 - h2
    sigma2.e <- 1 - H2

    # make the relatedness matrices
    relatedness <- makeRelatedness(pedigree)
    D <- relatedness[["Dr"]]
    A <- relatedness[["A"]]

    # fraction of additive (An) and dominance (Dn) relatedness matrix
    # corresponding to individuals in the trial
    An <- A[seq(m - n + 1, m), seq(m - n + 1, m)]
    Dn <- D[seq(m - n + 1, m), seq(m - n + 1, m)]

    # additive (G.add) and dominance (G.dom) genetic covariance matrix
    G.add <- sigma2.add * An
    G.dom <- sigma2.dom * Dn

    # sampling additive effects
    add <- rmvChol(
        n = 1,
        C = G.add,
        adj = FALSE
    )$MVN.col
    setnames(add, c("TreeID", "bv"))

    # sampling dominance effects
    dom <- rmvChol(
        n = 1,
        C = G.dom,
        adj = FALSE
    )$MVN.col
    setnames(dom, c("TreeID", "dd"))

    # merge additive and dominance effects
    BLUPs <- cbind(add[, .(bv)], dom[, .(dd)])
    trial <- cbind(trial, BLUPs)

    # sampling environment effects
    resid <- rnorm(n) * sqrt(sigma2.e) # e ~ N(0, sigma2.e)
    resid <- as.vector(scale(resid) * sqrt(sigma2.e))
    trial[, resid := resid]

    #  computing records (pheno) as y = mu + tbv + tdv + e
    trial[, mu := rep(mu, n)]
    trial[, pheno := mu + bv + dd + resid]

    # before applying correction by inbreeding depression ("id")
    # ensure that the phenotypic variance of the trait is 1
    trial[, pheno := mu + scale(pheno,
        center = TRUE, scale = TRUE
    ),
    by = .(gen)
    ]
    #  adjusting phenotypic values by inbreding depression ("id")
    trial[, id := -2.4375 * f]
    trial[, pheno := pheno + id]

    return(list(
        pedigree = pedigree,
        trial = trial
    ))
}
