makeCross <- function(parents,
					  gen, 
					  HDD=TRUE){
							  
###     create an Incomplete Diallel Design – no selfed, no reciprocal crosses	  ###	
males <- females <- parents  
diallel.tmp <- data.table(
    mum = rep(females, times= 1),
    dad = rep(males, each = length(females)))	
d0 <- diallel.tmp[, .(mum, dad)]
d2 <- d0[which(d0[, mum] >= d0[, dad]), ]  
d3 <- d2[order(d2[, mum], d2[, dad]), ]
parent.g <- d3[which(d3[, mum] != d3[, dad]), ] # ok...!

 if(HDD) {
     crossPlan.g <- parent.g   
	 return(crossPlan.g)
 } else {
     if(gen == 1) {
        weight.g <- 0.25
     } else if(gen == 2) { 
	    weight.g <- 0.3
     } else if(gen == 3) {
	    weight.g <- 0.4 
     } else {
	    weight.g <- 0.4
     } 
   npar.tmp <- length(unique(unlist(parent.g[, c(1, 2)])))	 
   nCrosses <- floor(weight.g*npar.tmp*(npar.tmp - 1L)/2L)
   crossPlan.g <- parent.g[sample(
	                1:dim(parent.g)[1], 
                    size = nCrosses), ] 									
   return(crossPlan.g)  				
   }
}



