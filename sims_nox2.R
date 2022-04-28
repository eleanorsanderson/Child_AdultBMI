

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
n = 150000 #number of individuals (should be 150000)
l = 150    #number of SNPs (for exposures - total)
lo = 10 #number of SNPs for outcome

beta1 = 0.3
beta2 = 0
beta3 = 0.2

results = NULL

#define effects outside the repetitions so they are consistent across the simulations
var1 <- 0.10/l
cor12 <- 0.3
cor13 <- 0.1
cor23 <- 0.3
mu <- as.vector(c(0,0,0))
sig <- matrix(c(var1, cor12*var1, cor13*var1, cor13*var1, var1, cor23*var1, cor13*var1, cor23*var1, var1), 3, 3)
effects <- mvrnorm(l, mu, sig)

#effs_g <- rnorm(l,0,sqrt(0.15/l))  
#effs_g2 <- 0.3*effs_g + rnorm((l),0,sqrt(0.15/l))
#effs_g3 <- 0.1*effs_g + 0.25*effs_g2 + 0.65*rnorm(l,0,sqrt(0.15/l))

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
  
  #d. model where there are three periods but only two included.
  
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
  res_nox2 <- data.frame("nox2", res)
  colnames(res_nox2)[1] <- ("sim")
  
  
  res_nox2$beta1_u <- beta1 + 0.1*beta2 + 0.01*beta3 + res_nox2$cor12*(sqrt(res_nox2$var2)/sqrt(res_nox2$var1))*beta2 + res_nox2$cor13*(sqrt(res_nox2$var3)/sqrt(res_nox2$var1))*beta3 
  res_nox2$beta2_u <- beta2 + res_nox2$cor12*(sqrt(res_nox2$var1)/sqrt(res_nox2$var2))*(beta1) + res_nox2$cor23*(sqrt(res_nox2$var3)/sqrt(res_nox2$var2))*beta3 
  
  res_nox2$beta1_m <- beta1 + res_nox2$pcor13_2*(sqrt(res_nox2$var3)/sqrt(res_nox2$var1))*beta3
  res_nox2$beta2_m <- beta2 + res_nox2$pcor23_1*(sqrt(res_nox2$var3)/sqrt(res_nox2$var2))*beta3
  
  
  results <- rbind(results,res_nox2)
  
}


save(results, file="sim_output_nox2.Rda")

