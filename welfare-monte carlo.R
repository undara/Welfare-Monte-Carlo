# This code runs the Monte Carlo simulations in 
# Identification and Inference for Welfare Gains without Unconfoundedness
# Written by: Undral Byambadalai, Boston University  

rm(list=ls())

# Set seed
set.seed(123)

# Setup 
n <- 5000 # sample size
N <- 5000 # number of replications
y_under <- 0 # a priori lower bound on y
y_bar <- 160000 # a priori upper bound on y
c <- 5 # strength of the instrument
gamma1 <-2000 #  coefficient for u in m1
gamma0 <-1000 # coefficient for u in m0

# Functions to be used
source("fns/gen_data.R")
source("fns/population_wc.R") 

#Population objects
eta0 <-population_wc(c)$eta0
eta1 <-population_wc(c)$eta1
p <-population_wc(c)$p
beta_0 <- population_wc(c)$beta_0

# Preallocate for simulation
CI_or <- matrix(0, nrow =N, ncol = 6)
CI_de <- matrix(0, nrow =N, ncol = 6)
CI_or_cf <- matrix(0, nrow =N, ncol = 6)
CI_de_cf <- matrix(0, nrow =N, ncol = 6)

# Run simulation
for (boot_iter in 1:N){
  #generate data
  jtpa <- gen_data(n)
  
  #estimate first-step nuisance functions
  eyd1x <- mean(jtpa$earnings[jtpa$edu==12 & jtpa$D==1])
  eyd0x <- mean(jtpa$earnings[jtpa$edu==12 & jtpa$D==0])
  pdx <- mean(jtpa$D[jtpa$edu==12])
  
  # prep
  term1 <- (eyd1x -y_bar)*pdx+ (y_under-eyd0x)*(1-pdx)
  term2 <- (eyd1x+eyd0x-(y_under+y_bar))*(1-pdx)
  term3 <- (eyd1x+eyd0x-(y_under+y_bar))*(0-pdx)
  term4 <- (jtpa$earnings[jtpa$edu==12 & jtpa$D==1]-eyd1x)
  term5 <- (jtpa$earnings[jtpa$edu==12 & jtpa$D==0]-eyd0x)
  
  # ------------
  # ORIGINAL 
  # ------------
  # beta_hat with original moment condition
  beta_hat_or <- 1/n*sum(jtpa$edu==12)* term1
  # variance with original moment condition
  V_hat_or <- 1/n* sum(jtpa$edu==12)*(term1-beta_hat_or)^2
  # Values
  CI_or[boot_iter,1] <-  round(beta_hat_or, digits = 0)
  CI_or[boot_iter,2] <- round(V_hat_or, digits = 0)
  CI_or[boot_iter,3] <- round(beta_hat_or-1.96*sqrt(V_hat_or)/sqrt(n), digits = 0)
  CI_or[boot_iter,4] <- round(beta_hat_or+1.96*sqrt(V_hat_or)/sqrt(n), digits = 0)               
  CI_or[boot_iter,5] <- round(CI_or[boot_iter,4]-CI_or[boot_iter,3], digits = 0)
  if (CI_or[boot_iter,3] < beta_0 & beta_0 < CI_or[boot_iter,4]) {
    CI_or[boot_iter,6] <- 1
  } else {
    CI_or[boot_iter,6] <- 0
  }
  
  # ------------
  # DEBIASED
  # ------------
  # beta_hat with debiased moment condition
  beta_hat_de <- beta_hat_or + 1/2*1/n*sum(term2+term4) + 1/2*1/n*sum(term3+term5)
  # variance with debiased moment condition
  V_hat_de <- 1/4*1/n*(sum((2*term1-2*beta_hat_de+term2+term4)^2)+sum((2*term1-2*beta_hat_de+term3+term5)^2))
  # Values
  CI_de[boot_iter,1] <- round(beta_hat_de, digits = 0)
  CI_de[boot_iter,2] <- round(V_hat_de, digits = 0)
  CI_de[boot_iter,3] <- round(beta_hat_de-1.96*sqrt(V_hat_de)/sqrt(n), digits = 0)
  CI_de[boot_iter,4] <- round(beta_hat_de+1.96*sqrt(V_hat_de)/sqrt(n), digits = 0)               
  CI_de[boot_iter,5] <- round(CI_de[boot_iter,4]-CI_de[boot_iter,3], digits = 0)
  if (CI_de[boot_iter,3] < beta_0 & beta_0 < CI_de[boot_iter,4]) {
    CI_de[boot_iter,6] <- 1
  } else {
    CI_de[boot_iter,6] <- 0
  }
  
  # ------------
  # ORIGINAL WITH CROSS-FITTING
  # ------------
  # Cross-fitting with L=2 
  L1=sort(sample(1:n,n/2))
  L2=setdiff(1:n,L1)
  data_L1 <-jtpa[L1,]
  data_L2 <-jtpa[L2,]
  
  #estimate first-step nuisance functions
  ## L1
  L1_eyd1x <- mean(data_L1$earnings[data_L1$edu==12 & data_L1$D==1])
  L1_eyd0x <- mean(data_L1$earnings[data_L1$edu==12 & data_L1$D==0])
  L1_pdx <- mean(data_L1$D[data_L1$edu==12])
  
  ## L2
  L2_eyd1x <- mean(data_L2$earnings[data_L2$edu==12 & data_L2$D==1])
  L2_eyd0x <- mean(data_L2$earnings[data_L2$edu==12 & data_L2$D==0])
  L2_pdx <- mean(data_L2$D[data_L2$edu==12])
  
  # prep - L1
  L1_term1 <- (L1_eyd1x -y_bar)*L1_pdx+ (y_under-L1_eyd0x)*(1-L1_pdx)
  L1_term2 <- (L1_eyd1x+L1_eyd0x-(y_under+y_bar))*(1-L1_pdx)
  L1_term3 <- (L1_eyd1x+L1_eyd0x-(y_under+y_bar))*(0-L1_pdx)
  L1_term4 <- (data_L2$earnings[data_L2$edu==12 & data_L2$D==1]-L1_eyd1x)
  L1_term5 <- (data_L2$earnings[data_L2$edu==12 & data_L2$D==0]-L1_eyd0x)
  
  #prep - L2
  L2_term1 <- (L2_eyd1x -y_bar)*L2_pdx+ (y_under-L2_eyd0x)*(1-L2_pdx)
  L2_term2 <- (L2_eyd1x+L2_eyd0x-(y_under+y_bar))*(1-L2_pdx)
  L2_term3 <- (L2_eyd1x+L2_eyd0x-(y_under+y_bar))*(0-L2_pdx)
  L2_term4 <- (data_L1$earnings[data_L1$edu==12 & data_L1$D==1]-L2_eyd1x)
  L2_term5 <- (data_L1$earnings[data_L1$edu==12 & data_L1$D==0]-L2_eyd0x)
  
  # beta_hat with original moment condition
  beta_hat_or_cf <- 1/n*(sum(data_L2$edu==12)*L1_term1 + sum(data_L1$edu==12)*L2_term1)
  
  # variance with original moment condition
  V_hat_or_cf <- 1/n*(sum(data_L2$edu==12)*(L1_term1-beta_hat_or_cf)^2 + sum(data_L1$edu==12)*(L2_term1-beta_hat_or_cf)^2)

  # Values
  CI_or_cf[boot_iter,1] <-  round(beta_hat_or_cf, digits = 0)
  CI_or_cf[boot_iter,2] <- round(V_hat_or_cf, digits = 0)
  CI_or_cf[boot_iter,3] <- round(beta_hat_or_cf-1.96*sqrt(V_hat_or_cf)/sqrt(n), digits = 0)
  CI_or_cf[boot_iter,4] <- round(beta_hat_or_cf+1.96*sqrt(V_hat_or_cf)/sqrt(n), digits = 0)               
  CI_or_cf[boot_iter,5] <- round(CI_or_cf[boot_iter,4]-CI_or_cf[boot_iter,3], digits = 0)
  if (CI_or_cf[boot_iter,3] < beta_0 & beta_0 < CI_or_cf[boot_iter,4]) {
    CI_or_cf[boot_iter,6] <- 1
  } else {
    CI_or_cf[boot_iter,6] <- 0
  }
  
  # ------------
  # DEBIASED WITH CROSS-FITTING
  # ------------
  # beta_hat with debiased moment condition
  beta_hat_de_cf <- 1/n* (sum(L1_term1+1/2*L1_term2+1/2*L1_term4) + sum(L1_term1+1/2*L1_term3+ 1/2*L1_term5) 
                          +sum(L2_term1+1/2*L2_term2+1/2*L2_term4) + sum(L2_term1+1/2*L2_term3+1/2*L2_term5))
  
  # variance with debiased moment condition
  V_hat_de_cf <- 1/4*1/n*(sum((2*L1_term1-2*beta_hat_de_cf+L1_term2+L1_term4)^2)+sum((2*L1_term1-2*beta_hat_de_cf+L1_term3+L1_term5)^2)
                      + sum((2*L2_term1-2*beta_hat_de_cf+L2_term2+L2_term4)^2)+sum((2*L2_term1-2*beta_hat_de_cf+L2_term3+L2_term5)^2))
  
  # Values
  CI_de_cf[boot_iter,1] <- round(beta_hat_de_cf, digits = 0)
  CI_de_cf[boot_iter,2] <- round(V_hat_de_cf, digits = 0)
  CI_de_cf[boot_iter,3] <- round(beta_hat_de_cf-1.96*sqrt(V_hat_de_cf)/sqrt(n), digits = 0)
  CI_de_cf[boot_iter,4] <- round(beta_hat_de_cf+1.96*sqrt(V_hat_de_cf)/sqrt(n), digits = 0)
  CI_de_cf[boot_iter,5] <- round(CI_de_cf[boot_iter,4]-CI_de_cf[boot_iter,3], digits = 0)
  if (CI_de_cf[boot_iter,3] < beta_0 & beta_0 < CI_de_cf[boot_iter,4]) {
    CI_de_cf[boot_iter,6] <- 1
  } else {
    CI_de_cf[boot_iter,6] <- 0
  }
}

# Calculate coverage probabilities
coverage_prob_or <- round(mean(CI_or[,6]), digits = 3)
coverage_prob_de <- round(mean(CI_de[,6]), digits = 3)
coverage_prob_or_cf <-round( mean(CI_or_cf[,6]), digits = 3)
coverage_prob_de_cf <-round( mean(CI_de_cf[,6]), digits = 3)

# Calculate average length of the confidence intervals
average_length_or <- round( mean(CI_or[,5]), digits = 0)
average_length_de <- round( mean(CI_de[,5]), digits = 0)
average_length_or_cf <- round(mean(CI_or_cf[,5]), digits = 0)
average_length_de_cf <- round(mean(CI_de_cf[,5]), digits = 0)
