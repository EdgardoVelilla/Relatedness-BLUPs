#' Function diallelSim                                                                                                                                                                         
#'                                                                                                                                                                                             
#' Author: Edgardo Velilla P.                                                                                                                                                                  
#' email{edgardo.velilla@cmpc.cl}                                                                                                                                                              
#' Created: 11-Apr-2022                                                                                                                                                                        
#' Modified: 23-Jan-2023                                                                                                                                                                       
#' License: GPLv3                                                                                                                                                                              
#'                                                                                                                                                                                             
#' General description:                                                                                                                                                                        
#'                                                                                                                                                                                             
#' This function enables a breeding program's simulation for a quantitative 
#' trait using stochastic simulation based on an infinitesimal model (Fisher
#' 1918; Bulmer 1980) with both additive and dominance genetic effects, but 
#' assuming that the epistatic sources of variance are of negligible 
#' significance. 
#' 
#' In the initial phase, founder individuals ("nfounder") are sampled from 
#' an ideal random mating base population —hence, one with no inbreeding and
#' Hardy-Weinberg linkage equilibrium. Thus, the founder's true additive and
#' dominance effects are taken at random from a normal distribution with mean
#' zero and additive (sigma2.add) and dominance (sigma2.dom) genetic variances
#' corresponding to the additive and non-additive genetic effects. In the base 
#' population, the phenotypic variance is set to be 1, so narrow-sense 
#' heritability (h2) equals the additive genetic variance (sigma2.add) and 
#' broad-sense heritability (H2) equals the dominance genetic variance 
#' (sigma2.dom) plus narrow-sense heritability (h2). Hence, non-genetic or 
#' environment variance (sigma2.e) equals 1 - H2. 
#' 
#' The first-generation simulation starts generating open-pollinated (OP) 
#' progenies representing "nfounder" half-sib families, each having "n.halfsib"
#' offsprings. In this first-generation trial, an experimental design was 
#' ignored mainly for simplicity. Then, selection proceeds by selecting the 
#' best m1 half-sib families in terms of their half-sib family means, where m1
#' depends on the intensity of selection ("p", ranging from 0 to 1) and the 
#' number of founder individuals ("nfounder"). 
#' 
#' Next, forward selection is carried out (i.e., selecting the best individual
#' from each of top m1 families). In the second and subsequent breeding 
#' generations, parents in the previous population are mated under a half-
#' diallel design ("HDD")(i.e., neither selfed nor reciprocal crosses). An 
#' individual's additive genetic effect (a) of these offspring and offspring
#' born in subsequent generations is computed as the average of its parental 
#' values plus an additive Mendelian sampling deviation taken at random from a
#' normal distribution with mean zero and variance equals 
#' 0.5(1 - ave.f)*sigma2.add, where "ave.f" is the average of the parents’ 
#' coefficient of inbreeding ("f", ranging from 0 to 1), which is computing
#' using the fast algorithm given by Meuwissen and Luo (1992). 
#' 
#' On the other hand, an individual's dominance genetic effect corrected for
#' the average effect of inbreeding on the mean is computed as d* = d + id, 
#' where d is the dominance effect ignoring inbreeding which is decomposed as
#' described by Hoeschele & VanRaden (1991), so d = df + dw, where df is the 
#' among-family effects between full-sib families drawn from a "multivariate 
#' normal distribution" (MVN) with mean zero and variance-covariance equals
#' 0.25sigma2.dom*F, with F representing the relatedness matrix between these 
#' full-sib families. In similar way, dw is the non-additive Mendelian sampling
#' deviation or the so-called within-family effect, taken from a normal 
#' distribution with mean zero and variance equals (1 - f)0.75*sigma2.dom, 
#' where "f" is the individual's inbreeding coefficient. The term, id, 
#' represents the inbreeding depression defined as id = b*f (de Boer and 
#' Hoeschele 1993), where b is the regression coefficient between the 
#' standardized phenotype, i.e., in units of phenotypic standard deviations
#' (Borralho 1994), which is based on results from Costa e Silva et al. (2011)
#' for DBH data and CMPC's measurements for 3rd generation E. globulus full-sib
#' progeny trials, so b equal 0.024375 units of phenotypic standard deviations 
#' per 1% of inbreeding, which equals the complete inbreeding depression taken
#' as at 0.57166% (Costa e Silva et al. 2001) of the mean DBH per 1% inbreeding.                                                                                                                                       
#'                                                                                                                                                                                             
#' Hence, a tree's record or phenotypic value ("pheno") is simulated as    
#'                                                                                                                     
#'                    pheno = mu + a + d + id + e                                                                                        
#'                                                                                                                                                                                             
#' where "mu" is the expected mean phenotypic value in the first generation, "a" 
#' and "d" are the total additive and dominance effects, "id" is the inbreeding
#' depression, and "e" is the environmental effect.                                                                                                                                                                       
#'                                                                                                                                                                                             
#' Selection proceeds by selecting individuals in a given generation to become 
#' parents of the next generation. Truncation selection is used by default, 
#' that is, the best-performing individuals are selected. The selection 
#' criteria can be based on individual's true breeding value ("tbv"), true 
#' dominance value ("tdv"), expected genetic value ("egv", i.e. "tbv" plus 
#' "tdv"), phenotypic value ("pheno") or randomly ("rnd").                                                                                                                                             
#'                                                                                                                                                                                             
#' References:                                                                                                                                                                                 
#'                                                                                                                                                                                             
#' Borralho N (1994). Heterogeneous selfing rates and dominance efects in 
#' estimating heritabilities from openpollinated progeny. Can J For Res 
#' 24:1079–1082.                                   
#'                                                                                                                                                                                             
#' Bulmer, M. G. (1980). The mathematical theory of quantitative genetics.
#' Oxford Univ. Press, Oxford, UK.                                                                                     
#'                                                                                                                                                                                             
#' Costa e Silva, J., Hardner, C., Tilyard, P. et al. (2011). The effects of
#' age and environment on the expression of inbreeding depression in Eucalyptus
#' globulus. Heredity 107, 50–60.       
#'                                                                                                                                                                                             
#' de Boer IJM, Hoeschele I (1993). Genetic evaluation methods for populations
#' with dominance and inbreeding. Theor Appl Genet 86: 245–258.                                                    
#'                                                                                                                                                                                             
#' Fisher, R. A. (1918). The correlation between relatives on the supposition 
#' of Mendelian inheritance. Trans. R. Soc. Edin. 52, 399-433.                                                      
#'                                                                                                                                                                                             
#' Hoeschele, I. and VanRaden, P.M. (1991). Rapid inverse of dominance 
#' relationship matrices for noninbred populations by including sire by dam 
#' subclass effects. Journal of Dairy Science 74, 557–569.                                                                                                                                                                                
#'                                                                                                                                                                                             
#' Meuwissen, T., Luo, Z. (1992). Computing inbreeding coefficients in large
#' populations. Genet Sel Evol 24, 305.                                                                              
#'  
#'                                                                                                                                                                                             
#' Arguments                                                                                                                                                                                   
#'                                                                                                                                                                                             
#' @param  nfounder                                                                                                                                                                            
#'                                                                                                                                                                                             
#' Integer number of founder individuals which conform the base population.                                                                                                                    
#'                                                                                                                                                                                             
#' @param n.halfsib                                                                                                                                                                            
#'                                                                                                                                                                                             
#' Integer number of offspring per half-sib family at first-generation.                                                                                                                        
#'                                                                                                                                                                                             
#' @param h2                                                                                                                                                                                   
#'                                                                                                                                                                                             
#' Numerical value for the narrow-sense heritability defined as the ratio of
#' additive genetic variance to the phenotypic variance.                                                             
#'                                                                                                                                                                                             
#' @param H2                                                                                                                                                                                   
#'                                                                                                                                                                                             
#' Numerical value for the broad-sense heritability defined as the ratio of the
#' sum of additive and dominance genetic variances to the phenotypic variance.                                    
#'                                                                                                                                                                                             
#' @mu                                                                                                                                                                                         
#'                                                                                                                                                                                             
#' Numerical value for the expected mean phenotypic value in the first 
#' generation.                                                                                                             
#'                                                                                                                                                                                             
#' @seed                                                                                                                                                                                       
#'                                                                                                                                                                                             
#' Integer number that specifies the initial value of the random-number seed 
#' that allows reproduction results. The seed is set by set.seed(seed, kind = 
#' "L'Ecuyer-CMRG").                      
#'                                                                                                                                                                                             
#' @selectTop                                                                                                                                                                                  
#'                                                                                                                                                                                             
#' Logical indicating if the highest values should be selected. Selects the 
#' lowest values if FALSE. The default is @selectTop=TRUE.                                                            
#'                                                                                                                                                                                             
#' @p                                                                                                                                                                                          
#'                                                                                                                                                                                             
#' Numerical value representing selection intensity expressed as the proportion
#' of individuals selected according to parameter @selectTop. For instance, 
#' p=0.05 means that the top 5% of the population is selected if @selectTop=TRUE.                                                                                                                                                  
#'                                                                                                                                                                                             
#' @sel                                                                                                                                                                                        
#'                                                                                                                                                                                             
#' A character indicating which criteria of selection will be used. Select on 
#' additive genetic values "tbv", dominance genetic deviation "tdv", expected 
#' genetic values "egv" (egv=tbv plus tdv), phenotypes "pheno", or randomly 
#' "rnd".                                                                                                                                                
#'                                                                                                                                                                                             
#' @HDD                                                                                                                                                                                        
#'                                                                                                                                                                                             
#' Logical indicating whether an incomplete diallel design should be created. 
#' If FALSE, only a random "proportion of crosses" from all combinations of the
#' incomplete diallel design will be done. The proportions are 0.25, 0.3 and 
#' 0.4 for the second, third and subsequent generations, respectively. The 
#' default is @HDD=TRUE.                                                       
#'                                                                                                                                                                                             
#' @ngen                                                                                                                                                                                       
#'                                                                                                                                                                                             
#' Integer number of generations (non-overlapping) to simulate.                                                                                                                                
#'                                                                                                                                                                                             
#' @adj                                                                                                                                                                                        
#'                                                                                                                                                                                             
#' Logical indicating if the properties of the sampled data-set from a 
#' multivariate normal distribution (MNV) match means and covariances "exactly".
#' adj= FALSE implies that the means and covariance values represent the 
#' population values. Hence, the sampled data-set most likely does not match 
#' these values exactly. The default is @adj= FALSE.                                 
#'                                                                                                                                                                                             
#' @n.blocks                                                                                                                                                                                   
#'                                                                                                                                                                                             
#' Integer number that specifies the number of blocks (replicates) of the field
#' design. If adj=TRUE, @n.blocks must be equal to the number of full-sib 
#' families in each generation plus one, 
#' i.e., @n.blocks >= n.parents*(n.parents - 1L)/2L + 1L. The default value is
#' @n.blocks= NULL, so the number of blocks is a function of parameters
#' @nfounder and @p.                          
#'                                                                                                                                                                                             
#'                                                                                                                                                                                             
#' @return                                                                                                                                                                                     
#'                                                                                                                                                                                             
#' A list with the following components:                                                                                                                                                       
#'                                                                                                                                                                                             
#' pedigree                                                                                                                                                                                    
#'                                                                                                                                                                                             
#' A data.table with dimension (mx5) and columns corresponding to:                                                                                                                             
#'                                                                                                                                                                                             
#' TreeID Integer for each individual's unique identifying code.                                                                                                                               
#' mum    Integer indicating each individual's mother.                                                                                                                                         
#' dad    Integer indicating each individual's father.                                                                                                                                         
#' f      Numeric value of each individual's inbreeding coefficient.                                                                                                                           
#' gen    Integer value of the generation in which each individual was born.                                                                                                                   
#'                                                                                                                                                                                             
#'                                                                                                                                                                                             
#' trial                                                                                                                                                                                       
#'                                                                                                                                                                                             
#' A data.table with with dimension (nx19) and columns corresponding to:                                                                                                                       
#'                                                                                                                                                                                             
#' TreeID      Same as in the pedigree.                                                                                                                                                        
#' mum         Same as in the pedigree.                                                                                                                                                        
#' dad         Same as in the pedigree.                                                                                                                                                        
#' block       Integer indicating the replication code (ranging from 1 to 
#'             n.blocks) according to the experimental design in each 
#'             generation.                                                   
#' gen         Same as in the pedigree.                                                                                                                                                        
#' f           Same as in the pedigree.                                                                                                                                                        
#' Xtype       A character indicating the type of cross used, "OP" and "CP"
#'             to denote open-pollinated and controlled-pollinated.                                                               
#' cross       A character identifying each crossing, such as cross = mumxdad
#'             = dadxmum.                                                                                                       
#' among.add   Numerical value for the additive among-family effect.                                                                                                                           
#' within.add  Numeric value of each individual’s additive within-family effect.                                                                                                               
#' tbv         Numeric value of each individual’s total additive genetic effect 
#'             (among.add + within.add).                                                                                      
#' among.dom   Numeric value of each individual’s non-additive among-family 
#'             effect.                                                                                                            
#' within.dom  Numeric value of each individual’s non-additive within-family 
#'             effect.                                                                                                           
#' tdv         Numeric value of each individual’s total dominance genetic 
#'             effect (among.dom + within.dom).                                                                                     
#' id          Numeric value of each individual’s inbreeding depression.                                                                                                                       
#' tdv.adj     Numeric value of each individual’s true dominance deviation 
#'             adjusted by inbreeding depression (tdv + id).                                                                       
#' egv         Numeric value of each individual’s expected genetic effect
#'             (tbv + tdv).                                                                                                         
#' resid       Numeric value of each individual’s residual (environmental)
#'             deviation.                                                                                                          
#' pheno       Numerical value of each individual’s phenotypic value (tbv + tdv
#'             + id + resid).                                                                                                 
#'                                                                                                                                                                                                                                                                                                                                                                                         
#' base.tgv                                                                                                                                                                                    
#'                                                                                                                                                                                             
#' A data.table with dimension (mx7) containing the "true genetic values" 
#' (including founder individuals) that will be used to calculate the BIAS and
#' RMSE according to different prediction approaches. It includes the following 
#' columns: TreeID, mum, dad, gen, tbv, tdv and tdv.adj that were previously 
#' described.                                                                  
#'                                                                                                                                                                                             
#'                                                                                                                                                                                             
#' F                                                                                                                                                                                           
#'                                                                                                                                                                                             
#' The relateddness matrix between full-sib families effect in sparse matrix 
#' format.                                                                                                           
#'                                                                                                                                                                                             
#'                                                                                                                                                                                             
#' initghost                                                                                                                                                                                   
#'                                                                                                                                                                                             
#' Integer indicating the initial ID code for ghost parents from open-
#' pollinated trial that will be used to build the relatedness matrix between
#' full-sib (F) on subsequent predictions.                                                                                       
#'                                                                                                                                                                                             
#'                                                                                                                                                                                             
#' inState & outState                                                                                                                                                                          
#'                                                                                                                                                                                             
#' Integer vectors containing the random number generator (RNG) state for 
#' random number generation in each simulation to reproduce a particular 
#' simulation.                                    
#'                                                                                                                                                                                             
#'                                                                                                                                                                                             
#' @example single simulation acroos six breeding cycles for a quantitative 
#' trait (DBH)                                                                                                        
#'                                                                                                                                                                                             
#' system.time(d1 <- diallelSim(nfounder= 300, n.halfsib= 100, h2= 0.2, 
#'                              H2= 0.4, mu= 4.36018, seed= 123, p= 0.05,                                                                              
#'                              sel= "egv", HDD= TRUE, selectTop= TRUE, 
#'                              ngen= 6, n.blocks= NULL))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            

diallelSim <- function(nfounder,
                       n.halfsib = NULL,
                       h2 = NULL,
                       H2 = NULL,
                       mu = NULL,
                       seed = NULL,
                       p = NULL,
                       sel = c("tbv", "tdv", "egv", "pheno", "rnd"),
                       HDD = TRUE,
                       selectTop = TRUE,
                       adj = TRUE,
                       ngen,
                       n.blocks = NULL) {
                        
    options(scipen = 999)
    Sys.setenv(LANG = "en")

    library(Matrix, quietly = TRUE)
    requireNamespace("MatrixModels", quietly = TRUE)
    library(data.table, quietly = TRUE)
    library(pedigree, quietly = TRUE)

    source("makeBasepop.r")
    source("selection.r")
    source("makeCross.r")
    source("rmvChol.r")
    source("makeF.r")
    source("matrixMethods.r")

    ## selection intensity should be between 0 and 1    
    if (!is.numeric(p) || p < 0 || p > 1) {
        stop("The selection intensity should be between 0 and 1.")
    }
    ## number of founders should be between 10 and 300    
    if (!is.numeric(nfounder) || nfounder < 10 || nfounder > 300) {
        stop("The number of founders should be between 10 and 300.")
    }
    if ((p == 0.01 & nfounder < 400) | (p == 0.05 & nfounder < 100)) {
        warning("the number of founders is too small to select parents for the next generations")
    }
    if (nfounder < 20) {
        stop("the number of founders is not enough")
    }

    ## Check and set initial seed    
    if (is.numeric(seed)) {
        seed <- as.integer(seed)
        set.seed(seed, kind = "L'Ecuyer-CMRG")
    } else if (!is.null(seed)) {
        set.seed(seed, kind = "L'Ecuyer-CMRG")
    } else {
        seed <- 123
        set.seed(seed, kind = "L'Ecuyer-CMRG")
    }

    ## saving the seed for reproducible results    
    inState <- .Random.seed
    #  cat('... seed inState ...', inState ,' ...\n')

    ## set ghost fathers for 1st open-pollinated generation's selection   
    n.parents <- ceiling(p * nfounder)
    n.cross <- n.parents * (n.parents - 1L) / 2L
    if (is.null(n.blocks)) n.blocks <- n.cross + 1L
    if (is.null(n.halfsib)) n.halfsib <- n.blocks
    n.trial <- nfounder * n.halfsib + n.blocks * n.cross * (ngen - 1L)
    if (n.trial < 1e5) {
        init.ghost <- 1e6L
    } else {
        init.ghost <- 1e8L
    }

    # Assumptions & Parameters:
    # founder individuals were derived from an unrelated and unselected base 
    # population
    # 100% survival is assumed so that all individuals have records
    # environmental correlation = 0 (no environment gradients within the trial)
    # phenotypic variance = 1. Hence, additive variance = h2 and dominance 
    # variance = H2 - h2

    if (is.null(h2)) h2 <- 0.2
    if (is.null(H2)) H2 <- 0.3
    if (is.null(mu)) mu <- 4.36018
    sigma2.add <- h2
    sigma2.dom <- H2 - h2
    sigma2.sca <- sigma2.dom / 4
    sigma2.e <- 1 - H2

    ## create the Base Population and the first-generation of 
    ## open-pollinated progeny trial   
    out.basegen1 <- makeBasepop(
        nfounder,
        init.ghost,
        n.halfsib,
        p,
        h2,
        H2,
        mu
    )

    base.gen1 <- out.basegen1$base.gen1
    basegen.g <- copy(base.gen1)
    ped.gen1 <- base.gen1[, .(TreeID, mum, dad, Xtype)]

    ##  base population's true genetic values    
    base.tgv <- out.basegen1$base.pop
    founders.tmp <- out.basegen1$base.founder
    founder.ghost <- founders.tmp[TreeID >= init.ghost, ]
    founder.ghost[, f := rep(0, times = dim(founder.ghost)[1])]
    founders.tmp <- founders.tmp[TreeID < init.ghost, ]
    cat("... generation ...", 1, " ... completed ...\n")

    ## loop through breeding generations    

    for (g in 2L:ngen) {
        ## parents selected for 2nd generation of controlled-pollinated trial
        if (g == 2L) {
            parent.g.tmp <- out.basegen1$sel[, TreeID]

            ## select parents for the subsequent generations	 
        } else {
            basegen.g <- whole.trial
            ped.old <- ped.jacquard
            parent.g.tmp <- selection(
                basegen.g,
                nfounder,
                p,
                sel,
                selectTop
            )[['parents']]
        }

        ## create an Incomplete diallel design – no selfed, 
        ## no reciprocal crosses    
        crossPlan.g <- makeCross(
            parent.g.tmp,
            g,
            HDD
        )

        ## create the pedigree for the current generation    
        crossPlan.oper <- crossPlan.g[,
            list(freq = rep(1L, n.blocks)),
            by = c("mum", "dad")
        ]
        crossPlan.oper[, freq := NULL]
        n.trees.g <- dim(crossPlan.oper)[1L]
        ped.g <- crossPlan.oper
        ped.g[, TreeID := as.integer(
            seq(
                max(basegen.g[, TreeID]) + 1L,
                max(basegen.g[, TreeID]) + n.trees.g
            )
        )]
        ped.g[, Xtype := rep("CP", times = n.trees.g)]
        setcolorder(ped.g, c("TreeID", "mum", "dad", "Xtype"))

        if (g == 2L) {
            ped.tmp <- rbind(
                founders.tmp[, c(1L:3L, 5L)],
                ped.gen1, ped.g
            )[order(TreeID), ]

            ## Inbreeding coefficient calculation (Meuwissen and Luo 1992)  
            ped.tmp[, f := calcInbreeding(ped.tmp[, c(1L:3L)])]
            ped.jacquard <- rbind(
                founder.ghost[, c(1L:3L, 5L:6L)],
                ped.tmp
            )
        } else {
            ped.tmp <- rbind(
                ped.old[, c(1L:4L)],
                ped.g
            )[TreeID < init.ghost, ][order(TreeID), ]
            ped.tmp[, f := calcInbreeding(ped.tmp[, c(1L:3L)])]
            ped.jacquard <- rbind(
                founder.ghost[, c(1L:3L, 5L:6L)],
                ped.tmp
            )
            setattr(ped.jacquard, "Jacquard", TRUE)
        }

        ## starting to build the 2nd-generation and subsequent CP generations
        design.g <- ped.g
        
         ## create the field "block" in the design
        design.g[, cross := makeFam(design.g)[, cross]]

        ## create the field "block"
        design.g[, block := rep(1:n.blocks,
            times = length(unique(cross))
        )]

        ## assign parent's breeding values (bv)    
        design.g[, bv_mum := basegen.g[, .(TreeID, tbv)][design.g[, .(mum)],
            on = .(TreeID = mum)
        ][, tbv]]
        design.g[, bv_dad := basegen.g[, .(TreeID, tbv)][design.g[, .(dad)],
            on = .(TreeID = dad)
        ][, tbv]]

        ## computing average parent's inbreeding coefficient (ave.f )   
        design.g[, f_mum := ped.jacquard[, .(TreeID, f)][design.g[, .(mum)], 
        on = .(TreeID = mum)][, f]]
        design.g[, f_dad := ped.jacquard[, .(TreeID, f)][design.g[, .(dad)], 
        on = .(TreeID = dad)][, f]]
        f_cols <- c("f_mum", "f_dad")
        design.g[, ave.f := rowMeans(.SD), by = TreeID, .SDcols = f_cols]

        ## sampling additive Mendelian Sampling term   
        ave.f <- design.g[, ave.f]
        ## w1 ~ N(0, 0.5*sigma2.add*(1 - ave.f))
        within.add <- rnorm(n.trees.g) * sqrt(0.5 * sigma2.add * (1 - ave.f)) 
        within.add <- as.vector(scale(within.add) * sqrt(0.5 * sigma2.add * (1 - ave.f)))
        design.g[, within.add := within.add]

        ## sampling non-additive Mendelian Sampling term    
        f <- ped.jacquard[, .(TreeID, f)][design.g, on = .(TreeID = TreeID)][, f]
        ## w2 ~ N(0, 0.75*sigma2.dom*(1 - f))
        within.dom <- rnorm(n.trees.g) * sqrt(0.75 * sigma2.dom * (1 - f)) 
        within.dom <- as.vector(scale(within.dom) * sqrt(0.75 * sigma2.dom * (1 - f)))
        design.g[, within.dom := within.dom]

        ## sampling environment effects 
        ## e ~ N(0, sigma2.e)     
        res <- rnorm(n.trees.g) * sqrt(sigma2.e) 
        resid <- as.vector(scale(res) * sqrt(sigma2.e))
        design.g[, resid := resid]

        ## calculation the "true breeding values" as 
        ## (tbv):= 0.5*add.mum + 0.5*add.dad + add Mendelian Sampling term   
        design.g[, among.add := 0.5 * bv_mum + 0.5 * bv_dad]
        design.g[, tbv := among.add + within.add]

        ## building the relatedness matrix between full-sib families, F, from
        ## Jacquard's identity coefficients    

        if (g < 3L) {
            F <- makeF(pedigree = ped.jacquard,
                init.ghost = init.ghost,
                gen = g
            )
        } else {
            F.out <- makeF(pedigree = ped.jacquard,
                init.ghost = init.ghost,
                gen = g
            )
            F <- F.out[["F"]]

            ###   get conditional number for matrix inversion   ###
            cond <- F.out[["cond"]]

            ###   examing whether the matrix F is positive definite, negative definite, and so on    
            status <- F.out[["status"]]
        }
        cat("job done...! \n")

        ###    number of full-sib families at the current("n.fs") and cumulative ("ntot.fs") generations    ###
        n.fs <- length(unique(design.g[
            !is.na(design.g[, cross]), cross
        ]))
        ntot.fs <- dim(F)[1L]

        ###    F.gen equals to F matrix at the current generation (at gen=2 all covariances of F are zero)    ###
        F.gen <- F[
            seq(
                from = ntot.fs - n.fs + 1L,
                to = ntot.fs
            ),
            seq(
                from = ntot.fs - n.fs + 1L,
                to = ntot.fs
            )
        ]

        ###    genetic variance-covariance matrix among full-sib families    ###
        G.gen <- sigma2.sca * F.gen

        ###    sampling the vector of sca effects    ###
        sca.g <- rnorm(n.fs) * sqrt(sigma2.sca)
        sca.g <- as.vector(scale(sca.g) * sqrt(sigma2.sca))

        if (g == 2L) {
            sca.old <- sca.g
        } else {
            sca.old <- c(sca.old, sca.g)
        }

        ###    sampling the among-family effects (sca) from a multivariate normal distribution    ###
        ##     sca ~ MVN(sca.g, sigma2.sca*F.gen)   ##

        among.g <- rmvChol(
            n = n.blocks,
            mu = sca.g,
            C = G.gen,
            adj = adj
        )$MVN.col
        setnames(among.g, c("cross", "among.dom"))

        ###    lookup sca effects by TreeID (= key)    ###
        setorder(among.g, cross)
        setorder(design.g, cross)
        among.g[, TreeID := design.g[, TreeID]]
        design.g[, among.dom := among.g[
            ,
            .(TreeID, among.dom)
        ][design.g[, .(TreeID)],
            on = .(TreeID = TreeID)
        ][, among.dom]]

        ###    estimate the effect of inbreeding depression ("id")    ###
        design.g[, f := ped.jacquard[, .(TreeID, f)][design.g[, .(TreeID)],
            on = .(TreeID = TreeID)
        ][, f]]
        design.g[, id := -2.4375 * f]

        ###    calculation the true dominance deviation ("tdv")    ###
        design.g[, tdv := among.dom + within.dom]

        ###    computing records (pheno) as y = mu + tbv + tdv + e    ###
        whole.trial <- design.g[
            ,
            pheno := mu + tbv + tdv + resid
        ]
        whole.trial[
            ,
            c("bv_mum", "bv_dad", "ave.f", "f_mum", "f_dad") := NULL
        ]

        # before applying correction by inbreeding depression ("id")
        # ensure that the phenotypic variance of the trait is 1
        mu.gen <- mean(whole.trial[, pheno])
        whole.trial[, pheno := mu.gen + scale(pheno,
            center = TRUE, scale = TRUE
        )]

        # adjusting phenotypic values by "id"
        whole.trial[, pheno := pheno + id]

        ###    adjusting dominance effects by "id"    ####
        design.g[, tdv.adj := tdv + id]
        design.g[, among.dom.adj := among.dom + id]

        ###    computing the expected genetic values ("egv")    ###
        setkey(design.g, TreeID)
        design.g[, egv := tbv + tdv + id]

        whole.trial[, gen := g]
        cols.gen1 <- names(base.gen1)
        setcolorder(whole.trial, cols.gen1)

        if (g == 2L) {
            whole.trial.old <- whole.trial
        } else {
            whole.trial.old <- rbind(
                whole.trial.old,
                whole.trial
            )
        }
        cat("be patient...\n")
        cat("... generation ...", g, " ... completed ...\n")
    }

    ###     end for g - 1 loop through generations      ###

    gen.add2 <- function(pedigree) {
        pedigree[, gen := 0L]
        for (i in 0:nrow(pedigree)) {
            pedigree[
                mum %in% TreeID[gen == i] | dad %in% TreeID[gen == i],
                gen := gen + 1L
            ]
            if (.Last.updated == 0L) break
        }
        pedigree[]
    }

    ped.jacquard <- gen.add2(ped.jacquard)

    # reconstruction the traditional pedigree
    pedigree <- ped.jacquard[TreeID < init.ghost, ]
    pedigree[, dad := ifelse(dad < init.ghost,
        dad,
        0L
    )]

    ##  ???                      ????
    setattr(pedigree, "Jacquard", NULL)

    ###     building the base of true genetic values (tgv)      ####
    trial.tgv <- whole.trial.old[, .(
        TreeID, mum, dad, gen, cross, among.add, within.add,
        tbv, among.dom, among.dom.adj, within.dom, tdv, id, tdv.adj
    )]
    fun <- function(x) {
        as.list(rep(0L, times = 6L))
    }
    cols <- c(
        "cross", "among.add", "within.add", "among.dom",
        "among.dom.adj", "within.dom", "id"
    )
    base.tgv <- base.tgv[, (cols) := fun(), by = TreeID][]
    cols.tgv <- names(trial.tgv)
    setcolorder(base.tgv, cols.tgv)

    basegen1.tgv <- base.gen1[, .(
        TreeID, mum, dad, gen, cross, among.add, within.add, tbv,
        among.dom, among.dom.adj, within.dom, tdv, id, tdv.adj
    )]
    
    basefull.tgv <- rbind(base.tgv, basegen1.tgv, trial.tgv)
    trial <- rbind(base.gen1, whole.trial.old)
    param <- data.table(
        nfounder = nfounder,
        n.halfsib = n.halfsib,
        h2 = h2,
        H2 = H2,
        mu = mu,
        p = p,
        sel = sel,
        HDD = HDD,
        selectTop = selectTop,
        adj = adj,
        ngen = ngen,
        n.blocks = n.blocks,
        n.cross = n.cross,
        n.parents = n.parents
    )
    outState <- .Random.seed
    # cat('... seed outState ...', outState ,' ...\n')

    return(list(
        trial = trial,
        base.tgv = basefull.tgv,
        ped.jacquard = ped.jacquard,
        pedigree = pedigree,
        F = F,
        cond = cond,
        status = status,
        sca = sca.old,
        init.ghost = init.ghost,
        parameters = param,
        inState = inState,
        outState = outState
    ))
}
