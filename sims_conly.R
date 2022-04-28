

rm(list = ls(all=TRUE))

set.seed(4)


library(MASS)
source('MRest.R')

make_geno <- function(nid, nsnp, af)
{
  return(matrix(rbinom(nid * nsnp, 2, af), nid, nsnp))
}


reps = 2
n = 150000 #number of individuals
l = 150    #number of SNPs (for exposures - total)
lo = 50 #number of SNPs for outcome

results = NULL

#define effects outside the repetitions so they are consistent across the simulations
effs_g <- rnorm(l,0,sqrt(0.15/l))  
effs_g2 <- 0.3*effs_g + rnorm((l),0,sqrt(0.15/l))
effs_out <- rnorm(lo,0,sqrt(0.3/l))
effs_c1 <- 0.5
effs_c2 <- 0.5

for(i in 1:reps){
  
  
  #variables
  g <- make_geno(n,(l+lo),0.5)
  Sigma = matrix(c(1, 0.8,0.8,0.8,1,0.8,0.8,0.8,1), nrow = 3)
  u <- mvrnorm(n, c(1,1,1), Sigma)
  ua = u[,1]
  u2a = u[,2]
  u3a = u[,3]
  
  gb <- make_geno(n,(l+lo),0.5)
  u <- mvrnorm(n, c(1,1,1), Sigma)
  ub = u[,1]
  u2b= u[,2]
  u3b = u[,3]
  
  #c. model where y -> x2
  
  #exposures
  L1 <- g[,1:l]%*%effs_g 
  L2 <- g[,1:l]%*%effs_g2
  
  x1 <- L1 + effs_c1*ua + rnorm(n,0,1)
  ya <- 0.2*x1 +g[,(l+1):(l+lo)]%*%effs_out + effs_c1*ua + effs_c2*u2a + rnorm(n,0,1)
  x2 <- L2 + effs_c2*u2a + 0.1*x1 + 0.5*ya + 0.5*rnorm(n,0,1) 
  
  L1b <- gb[,1:l]%*%effs_g 
  x1b <- L1b + effs_c1*ub + rnorm(n,0,1)
  L2b <- gb[,1:l]%*%effs_g2 
  
  y <- 0.2*x1b +gb[,(l+1):(l+lo)]%*%effs_out + effs_c1*ub + effs_c2*u2b + rnorm(n,0,1)
  x2b <- L2b + effs_c2*u2b + 0.1*x1b + 0.5*y + 0.5*rnorm(n,0,1) 
  
  res <- MRest()
  resc <- data.frame("c", res)
  colnames(resc)[1] <- ("sim")
  
  resc$beta1_u <- 0.2
  resc$beta2_u <- 0 + cor(effs_g, effs_g2)*0.2 
  resc$beta1_m <- 0.2 
  resc$beta2_m <- 0
  
  results <- rbind(results,resc)
  
}


#save(results, file="sim_output_c.Rda")

