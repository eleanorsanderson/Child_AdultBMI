
MRest_x1x3 <- function()
{
  
  mvdat <- data.frame()
  for(s in 1:(l+lo)){
    x1_res <- summary(lm(x1~g[,s]))
    x2_res <- summary(lm(x2~g[,s]))
    x3_res <- summary(lm(x3~g[,s]))
    y_res <- summary(lm(y~gb[,s]))
    mvdat[s,1] <- (x1_res$coefficient[2,1])
    mvdat[s,2] <- (x3_res$coefficient[2,1])
    mvdat[s,3] <- (y_res$coefficient[2,1])
    
    mvdat[s,4] <- (x1_res$coefficient[2,2])
    mvdat[s,5] <- (x3_res$coefficient[2,2])
    mvdat[s,6] <- (y_res$coefficient[2,2])
    
    mvdat[s,7] <- (x1_res$coefficient[2,4])
    mvdat[s,8] <- (x3_res$coefficient[2,4])
    mvdat[s,9] <- (y_res$coefficient[2,4])
    
    mvdat[s,10] <- (x2_res$coefficient[2,1])
    
  }
  
  colnames(mvdat) <- c("exposure_beta.x1", "exposure_beta.x3", "outcome_beta","exposure_se.x1", "exposure_se.x3", 
                       "outcome_se","exposure_pval.x1", "exposure_pval.x3", "outcome_pval",  "exposure_beta.x2" )
 
 dat1 <- subset(mvdat, mvdat$exposure_pval.x1 < 5e-08)
 unix1 <- summary(lm(dat1$outcome_beta ~ -1 + dat1$exposure_beta.x1, weights = 1/(dat1$outcome_se^2)))
 
 F_1 <- (dat1$exposure_beta.x1/dat1$exposure_se.x1)^2
 
 dat2 <- subset(mvdat, mvdat$exposure_pval.x3 < 5e-08)
 unix2 <- summary(lm(dat2$outcome_beta ~ -1 + dat2$exposure_beta.x3, weights = 1/(dat2$outcome_se^2)))
 F_2 <- (dat2$exposure_beta.x3/dat2$exposure_se.x3)^2
 
 mvdat$minp <- apply(mvdat[,c("exposure_pval.x1","exposure_pval.x3")], 1, min)
  mv <- subset(mvdat, mvdat$minp < 5e-08)
  mvmr <- summary(lm(mv$outcome_beta ~ -1 + mv$exposure_beta.x1 + mv$exposure_beta.x3, weights = 1/(mv$outcome_se^2)))
  
  rho = cor(x1,x3)
  sig12 = as.vector(rho)*mv$exposure_se.x1*mv$exposure_se.x3
  
  delta1 <- lm(mv$exposure_beta.x1~ -1 + mv$exposure_beta.x3)$coefficients[1]
  delta2 <- lm(mv$exposure_beta.x3~ -1 + mv$exposure_beta.x1)$coefficients[1]
  
  Qind_1 <- ((mv$exposure_beta.x1 - delta1*mv$exposure_beta.x3)^2)/(mv$exposure_se.x1^2 + (delta1^2)*mv$exposure_se.x3^2 - 2*delta1*sig12)
  Qind_2 <- ((mv$exposure_beta.x3 - delta2*mv$exposure_beta.x1)^2)/(mv$exposure_se.x3^2 + (delta2^2)*mv$exposure_se.x1^2 - 2*delta2*sig12)
  
  snps_outx1 <- as.numeric(sum((mvdat$exposure_pval.x1[(l+1):(l+lo)]<5e-08)))
  snps_outx3 <- as.numeric(sum((mvdat$exposure_pval.x3[(l+1):(l+lo)]<5e-08)))
  
  #pcor12 <- pcor.test(mv$exposure_beta.x2, mv$exposure_beta.x1, mv$exposure_beta.x3)$estimate
  #pcor23 <- pcor.test(mv$exposure_beta.x2, mv$exposure_beta.x3, mv$exposure_beta.x1)$estimate
  
out <- data.frame(rho, cor(mv$exposure_beta.x1,mv$exposure_beta.x2),   cor(mv$exposure_beta.x1,mv$exposure_beta.x3),   cor(mv$exposure_beta.x2,mv$exposure_beta.x3),
                  pcor.test(mv$exposure_beta.x1,mv$exposure_beta.x2,mv$exposure_beta.x3)$estimate, pcor.test(mv$exposure_beta.x3,mv$exposure_beta.x2,mv$exposure_beta.x1)$estimate, 
                  var(mv$exposure_beta.x1), var(mv$exposure_beta.x2),var(mv$exposure_beta.x3),
                   unix1$coefficients[1,1],unix1$coefficients[1,2], nrow(dat1), mean(F_1), 
                    unix2$coefficients[1,1],unix2$coefficients[1,2],nrow(dat2), mean(F_2), 
                    mvmr$coefficients[1,1],mvmr$coefficients[1,2], sum(Qind_1)/(nrow(mv)-1), 
                    mvmr$coefficients[2,1],  mvmr$coefficients[2,2], sum(Qind_2)/(nrow(mv)-1), nrow(mv),
                    snps_outx1, snps_outx3)
  
  

colnames(out) <- c("rho", "cor12", "cor13", "cor23",
                   "pcor12_3", "pcor23_1",
                   "var1", "var2", "var3",
                   "uni_x1_b","uni_x1_se", "uni_x1_nsnp","F_x1" , 
                   "uni_x3_b","uni_x3_se", "uni_x3_nsnp", "F_x3", 
                   "mv_x1_b", "mv_x1_se", "CF_x1", 
                   "mv_x3_b","mv_x3_se", "CF_x3","mv_nsnp", 
                   "snps_outx1", "snps_outx3")
  
return(out)

}

