## Title: Identification of welfare differences of policies
## Written by: Undral Byambadalai
## Date: 02/26/2019

## this code calculates 
## benchmark treatment rule
## upper and lower bounds on welfare regret
## under mean-independence IV
## when X=edu is used as a covariate

################################ read data #######################################################
## data from Abadie, Angrist, Imbens (2002)
jtpa_aai = read.table('jtpa.tab', header = TRUE)
colnames(jtpa_aai)[1] <- "recid" #record ids
colnames(jtpa_aai)[3] <- "Z" # random assignment
colnames(jtpa_aai)[4] <- "D"  # program enrollment
jtpa_aai <- jtpa_aai[-c(2,5:19)]  # remove other covariates

## data from Kitagawa, Tetenov (2018)
jtpa_kt= read.table('jtpa_kt.tab', header = TRUE)
jtpa_kt <- jtpa_kt[-c(2)] # remove "D"
jtpa <- merge(jtpa_aai, jtpa_kt, by="recid") #merge them

#################################################################################################
# set ybar (upper bound on outcome)
ybar <- 160000

## Use X=edu as a covariate and find benchmark treatment rule under Manski assumptions
## E[Y|D=d, X=x, Z=z]
yhat_iv <- function(d,x,z){ if (sum(jtpa$earnings[jtpa$D==d & jtpa$Z==z & jtpa$edu==x])==0) return (0)
  else return(mean(jtpa$earnings[jtpa$D==d & jtpa$Z==z & jtpa$edu==x]))}

## P[D=1|X=x, Z=z]=E[D|X=x, Z=z]
phat_iv <- function(x,z){mean(jtpa$D[jtpa$Z==z &jtpa$edu==x], na.rm = TRUE)}

#benchmark rule:
delta_iv <- function(x){max(yhat_iv(1,x,1)*phat_iv(x,1),yhat_iv(1,x,0)*phat_iv(x,0))-max(yhat_iv(0,x,1)*(1-phat_iv(x,1)), yhat_iv(0,x,0)*(1-phat_iv(x,0)))}

# upper bound on welfare:
up_pos_iv <-function(x){ mean(jtpa$edu==x)*min((yhat_iv(1,x,1)*phat_iv(x,1) + (ybar - yhat_iv(0,x,1))*(1-phat_iv(x,1))), (yhat_iv(1,x,0)*phat_iv(x,0) + (ybar - yhat_iv(0,x,0))*(1-phat_iv(x,0))))} #used when delta*=1, delta=0
up_neg_iv <- function(x){-mean(jtpa$edu==x)*max(((yhat_iv(1,x,1)-ybar)*phat_iv(x,1)+(-yhat_iv(0,x,1))*(1-phat_iv(x,1))),((yhat_iv(1,x,0)-ybar)*phat_iv(x,0)+(-yhat_iv(0,x,0))*(1-phat_iv(x,0))))} #used when delta*=0, delta=1

# lower bound on welfare:
low_pos_iv <- function(x){mean(jtpa$edu==x)*max(((yhat_iv(1,x,1)-ybar)*phat_iv(x,1) + (- yhat_iv(0,x,1))*(1-phat_iv(x,1))),((yhat_iv(1,x,0)-ybar)*phat_iv(x,0) + (- yhat_iv(0,x,0))*(1-phat_iv(x,0))))} #used when delta*=1, delta=0
low_neg_iv <-function(x){-mean(jtpa$edu==x)*min((yhat_iv(1,x,1)*phat_iv(x,1)+(ybar-yhat_iv(0,x,1))*(1-phat_iv(x,1))),(yhat_iv(1,x,0)*phat_iv(x,0)+(ybar-yhat_iv(0,x,0))*(1-phat_iv(x,0))))} #used when delta*=0, delta=1

## benchmark: treat if x leq 12 (until high school)
## treatment: treat if x leq 16 (until college)
up_neg_iv(13)+up_neg_iv(14)+up_neg_iv(15)+up_neg_iv(16)
low_neg_iv(13)+low_neg_iv(14)+low_neg_iv(15)+low_neg_iv(16)

## 
edu <- c(7:18)
lapply(up_pos_iv, edu)
