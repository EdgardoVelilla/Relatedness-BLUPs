!WORKSPACE 120000 !ARGS  1
BLUP_Analysis_2023_oubred

 TreeID !P
 TreeID2 !P
 mum !I
 dad !I
 gen *
 Xtype !A
 cross !A
 pheno

pedigree.csv !SKIP 1
Ainv.giv !SKIP 1
Finv.giv !SKIP 1
Minv.giv !SKIP 1
trial.csv  !SKIP 1 !CSV !MAXIT 1 !DOPART $A !NOREORDER !DISPLAY 2

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
#         control-pollinated families from an outbred population
#
# sigma2.e <- 120
# sigma2.add <- 90
# sigma2.dom  <- 80
# sigma2.sca <- sigma2.dom / 4
#
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

!PART 1 # BLUP analysis with gamma parametrization
pheno ~ mu !r TreeID cross TreeID2
1 1 3
0 0 IDEN !s2== 1
TreeID 1
TreeID 0 GIV1 0.75 !GF
cross 1
cross 0 GIV2 0.1666667 !GF
TreeID2 1
TreeID2 0 GIV3 0.6666667 !GF


!PART 3 # BLUP analysis with gamma parametrization
pheno ~ mu !r TreeID at(Xtype,2).cross TreeID2
1 1 3
0 0 IDEN !s2== 1
TreeID 1
TreeID 0 GIV1 0.75 !GF
at(Xtype,2).cross 1
at(Xtype,2).cross 0 GIV2 0.1666667 !GF
TreeID2 1
TreeID2 0 GIV3 0.6666667 !GF

