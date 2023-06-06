# !WORKSPACE 120000 !ARGS 4
BLUP_Analysis_2023_inbred

 TreeID !P
 TreeID2 !P
 mum !I
 dad !I
 gen
 f
 Xtype !A
 cross !A
 pheno

pedigree.csv !SKIP 1
Ainv.giv !SKIP 1
Finv.giv !SKIP 1
Minv.giv !SKIP 1
trial.csv  !SKIP 1 !CSV !MAXIT 1 !DOPART $A !NOREORDER !DISPLAY 2

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
#          control-pollinated families from an inbred population
#
# sigma2.e <- 0.5
# sigma2.add <- 0.25
# sigma2.dom  <- 0.25
# sigma2.sca <- sigma2.dom / 4
#
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
!PART 1 # BLUP analysis with gamma parametrization
pheno ~ mu !r TreeID cross TreeID2
1 1 3
0 0 IDEN !s2== 1
TreeID 1
TreeID 0 GIV1 0.5 !GF
cross 1
cross 0 GIV2 0.125 !GF
TreeID2 1
TreeID2 0 GIV3 0.5 !GF


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
#          control-pollinated families from an inbred population
#
# sigma2.e <- 0.5
# sigma2.add <- 0.25
# sigma2.dom  <- 0.25
# sigma2.sca <- sigma2.dom / 4
#
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
!PART 2 # BLUP analysis with gamma parametrization
pheno ~ mu f !r TreeID cross TreeID2
1 1 3
0 0 IDEN !s2== 1
TreeID 1
TreeID 0 GIV1 0.5 !GF
cross 1
cross 0 GIV2 0.125 !GF
TreeID2 1
TreeID2 0 GIV3 0.5 !GF


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
#    open- and control-pollinated families from an inbred population
#
# sigma2.e <- 0.5
# sigma2.add <- 0.25
# sigma2.dom  <- 0.25
# sigma2.sca <- sigma2.dom / 4
#
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
!PART 3 # BLUP analysis with gamma parametrization
pheno ~ mu !r TreeID at(Xtype,2).cross TreeID2
1 1 3
0 0 IDEN !s2== 1
TreeID 1
TreeID 0 GIV1 0.5 !GF
at(Xtype,2).cross 1
at(Xtype,2).cross 0 GIV2 0.125 !GF
TreeID2 1
TreeID2 0 GIV3 0.5 !GF

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
#    open- and control-pollinated families from an inbred population
#
# sigma2.e <- 0.5
# sigma2.add <- 0.25
# sigma2.dom  <- 0.25
# sigma2.sca <- sigma2.dom / 4
#
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
!PART 4 # BLUP analysis gamma with parametrization
pheno ~ mu f !r TreeID at(Xtype,2).cross TreeID2
1 1 3
0 0 IDEN !s2== 1
TreeID 1
TreeID 0 GIV1 0.5 !GF
at(Xtype,2).cross 1
at(Xtype,2).cross 0 GIV2 0.125 !GF
TreeID2 1
TreeID2 0 GIV3 0.5 !GF



