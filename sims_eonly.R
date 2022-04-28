

rm(list = ls(all=TRUE))

set.seed(4)


library(MASS)
library(ppcor)
source('MRest_three.R')

make_geno <- function(nid, nsnp, af)
{
  return(matrix(rbinom(nid * nsnp, 2, af), nid, nsnp))
}


reps = 1000
n = 150000 #number of individuals
l = 150    #number of SNPs (for exposures - total)
lo = 10 #number of SNPs for outcome

results = NULL

beta1 = 0.3
beta2 = 0.2
beta3 = 0.2

#define effects outside the repetitions so they are consistent across the simulations
var1 <- 0.10/l
cor12 <- 0.3
cor13 <- 0.05
cor23 <- 0.05
mu <- as.vector(c(0,0,0))
sig <- matrix(c(var1, cor12*var1, cor13*var1, cor13*var1, var1, cor23*var1, cor13*var1, cor23*var1, var1), 3, 3)
effects <- mvrnorm(l, mu, sig)


effs_g <- (effects[,1])
effs_g2 <- (effects[,2])
effs_g3 <- (effects[,3])


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
  
  L1 <- g[,1:l]%*%effs_g 
  L2 <- g[,1:l]%*%effs_g2
  L3 <- g[,1:l]%*%effs_g3
  
  x1 <- L1 + effs_c1*ua + rnorm(n,0,1)
  x2 <- L2 + effs_c2*u2a + 0.1*x1 + 0.9*rnorm(n,0,1) 
  x3 <- L3 + effs_c2*u3a + 0.1*x2 + 0.9*rnorm(n,0,1) 
  
  L1b <- gb[,1:l]%*%effs_g 
  L2b <- gb[,1:l]%*%effs_g2 
  L3b <- gb[,1:l]%*%effs_g3
  
  x1b <- L1b + effs_c1*ub + rnorm(n,0,1)
  x2b <- L2b + effs_c2*u2b + 0.1*x1b + 0.9*rnorm(n,0,1) 
  x3b <- L3b + effs_c2*u3b + 0.1*x2b + 0.9*rnorm(n,0,1) 
  
  #outcome
  y <- beta1*x1b + beta2*x2b + beta3*x3b + gb[,(l+1):(l+lo)]%*%effs_out + 0.3*effs_c1*ub + 0.3*effs_c2*u2b + 0.3*effs_c2*u3b
  
  res <- MRest_three()
  rese <- data.frame("e", res)
  colnames(rese)[1] <- ("sim")
  
  
  rese$beta1_u <- beta1 + 0.1*beta2 + 0.01*beta3 + rese$d1cor12*(sqrt(rese$d1var2)/sqrt(rese$d1var1))*beta2 + rese$d1cor13*(sqrt(rese$d1var3)/sqrt(rese$d1var1))*beta3 
  rese$beta2_u <- beta2 + rese$d2cor12*(sqrt(rese$d2var1)/sqrt(rese$d2var2))*(beta1) + rese$cor23*(sqrt(rese$d2var3)/sqrt(rese$d2var2))*beta3 
  
  rese$beta1_m <- beta1 + rese$pcor13_2*(sqrt(rese$var3)/sqrt(rese$var1))*beta3
  rese$beta2_m <- beta2 + rese$pcor23_1*(sqrt(rese$var3)/sqrt(rese$var2))*beta3
  
  
  results <- rbind(results,rese)
  
}


save(results, file="sim_output_e.Rda")

