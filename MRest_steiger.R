
MRest_steiger <- function()
{
  
  mvdat <- data.frame()
  for(s in 1:(l+lo)){
    x1_res <- summary(lm(x1~g[,s]))
    x2_res <- summary(lm(x2~g[,s]))
    y_res <- summary(lm(y~gb[,s]))
    mvdat[s,1] <- (x1_res$coefficient[2,1])
    mvdat[s,2] <- (x2_res$coefficient[2,1])
    mvdat[s,3] <- (y_res$coefficient[2,1])
    
    mvdat[s,4] <- (x1_res$coefficient[2,2])
    mvdat[s,5] <- (x2_res$coefficient[2,2])
    mvdat[s,6] <- (y_res$coefficient[2,2])
    
    mvdat[s,7] <- (x1_res$coefficient[2,4])
    mvdat[s,8] <- (x2_res$coefficient[2,4])
    mvdat[s,9] <- (y_res$coefficient[2,4])
    
  }
  
  colnames(mvdat) <- c("exposure_beta.x1", "exposure_beta.x2", "outcome_beta","exposure_se.x1", "exposure_se.x2", 
                       "outcome_se","exposure_pval.x1", "exposure_pval.x2", "outcome_pval" )
  
 #steiger filtering
  mvdat$diff <- as.numeric(mvdat$outcome_pval < mvdat$exposure_pval.x2)
  mvdat$x2small <- as.numeric(mvdat$exposure_pval.x2 <  5e-08)
  mvdat$both <- as.numeric(mvdat$diff==1 & mvdat$x2small ==1)
  mvdat <- mvdat[(mvdat$both==0),]

 dat1 <- subset(mvdat, mvdat$exposure_pval.x1 < 5e-08)

 unix1 <- summary(lm(dat1$outcome_beta ~ -1 + dat1$exposure_beta.x1, weights = 1/(dat1$outcome_se^2)))
 
 F_1 <- (dat1$exposure_beta.x1/dat1$exposure_se.x1)^2
 
 dat2 <- subset(mvdat, mvdat$exposure_pval.x2 < 5e-08)

 unix2 <- summary(lm(dat2$outcome_beta ~ -1 + dat2$exposure_beta.x2, weights = 1/(dat2$outcome_se^2)))
 F_2 <- (dat2$exposure_beta.x2/dat2$exposure_se.x2)^2
 
 mvdat$minp <- apply(mvdat[,c("exposure_pval.x1","exposure_pval.x2")], 1, min)
  mv <- subset(mvdat, mvdat$minp < 5e-08)
  mvmr <- summary(lm(mv$outcome_beta ~ -1 + mv$exposure_beta.x1 + mv$exposure_beta.x2, weights = 1/(mv$outcome_se^2)))
  
  rho = cor(x1,x2)
  sig12 = as.vector(rho)*mv$exposure_se.x1*mv$exposure_se.x2
  
  delta1 <- lm(mv$exposure_beta.x1~ -1 + mv$exposure_beta.x2)$coefficients[1]
  delta2 <- lm(mv$exposure_beta.x2~ -1 + mv$exposure_beta.x1)$coefficients[1]
  
  Qind_1 <- ((mv$exposure_beta.x1 - delta1*mv$exposure_beta.x2)^2)/(mv$exposure_se.x1^2 + (delta1^2)*mv$exposure_se.x2^2 - 2*delta1*sig12)
  Qind_2 <- ((mv$exposure_beta.x2 - delta2*mv$exposure_beta.x1)^2)/(mv$exposure_se.x2^2 + (delta2^2)*mv$exposure_se.x1^2 - 2*delta2*sig12)
  
  snps_outx1 <- as.numeric(sum((mvdat$exposure_pval.x1[(l+1):(l+lo)]<5e-08)))
  snps_outx2 <- as.numeric(sum((mvdat$exposure_pval.x2[(l+1):(l+lo)]<5e-08)))
  
  out <- data.frame(rho, unix1$coefficients[1,1],unix1$coefficients[1,2], nrow(dat1), mean(F_1), 
                    unix2$coefficients[1,1],unix2$coefficients[1,2],nrow(dat2), mean(F_2), 
                    mvmr$coefficients[1,1],mvmr$coefficients[1,2], sum(Qind_1)/(nrow(mv)-1), 
                    mvmr$coefficients[2,1],  mvmr$coefficients[2,2], sum(Qind_2)/(nrow(mv)-1), nrow(mv),
                    snps_outx1, snps_outx2)
  
 
colnames(out) <- c("rho", "uni_x1_b","uni_x1_se", "uni_x1_nsnp","F_x1" , 
                   "uni_x2_b","uni_x2_se", "uni_x2_nsnp", "F_x2", 
                       "mv_x1_b",     "mv_x1_se", "CF_x1", 
                   "mv_x2_b","mv_x2_se",  "CF_x2","mv_nsnp",
                   "snps_outx1", "snps_outx2")
  
return(out)

}

