## Title: Identification of welfare differences of policies
## Written by: Undral Byambadalai
## Date: 03/04/2019

## this code calculates 
## benchmark treatment rule
## upper and lower bounds on welfare regret
## under 
  # (1) worst-case
  # (2) MTR
  # (3) IV-worst case
  # (4) IV-MTR
  # (5) MIV-worst case
  # (6) MIV-MTR

## when Y=post-program earnings;
## D=program enrollment;
## Z=random assignement;
## X=education;

library(stargazer)

##################################################################################################
################################ read data #######################################################
##################################################################################################
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

#################################################################################################
################################ moments ########################################################
#################################################################################################
## E[Y|D=d, X=x]
eydx <- function(d,x){ mean(jtpa$earnings[jtpa$D==d & jtpa$edu==x])} 

## P[D=1|X=x]=E[D|X=x]
pdx <- function(x){mean(jtpa$D[jtpa$edu==x])} 

## E[Y|X=x]
eyx <- function(x){mean(jtpa$earnings[jtpa$edu==x])} 

## E[Y|D=d, X=x, Z=z]
eydxz <- function(d,x,z){ if (sum(jtpa$earnings[jtpa$D==d & jtpa$Z==z & jtpa$edu==x])==0) return (0)
  else return(mean(jtpa$earnings[jtpa$D==d & jtpa$Z==z & jtpa$edu==x]))}

## P[D=1|X=x, Z=z]=E[D|X=x, Z=z]
pdxz <- function(x,z){mean(jtpa$D[jtpa$Z==z &jtpa$edu==x], na.rm = TRUE)}

## E[Y|X=x, Z=z]
eyxz <- function(x,z){ if (sum(jtpa$earnings[jtpa$Z==z & jtpa$edu==x])==0) return (0)
  else return(mean(jtpa$earnings[jtpa$Z==z & jtpa$edu==x]))}

## P(Z=1)
pz <- mean(jtpa$Z==1)

#################################################################################################
################################ (1) worst-case #################################################
#################################################################################################
#benchmark rule (treat if delta(x)>0; not treat otherwise)
delta <- function(x){eydx(1,x)*pdx(x)-eydx(0,x)*(1-pdx(x))} 


## upper bound on welfare:
## P(X)*(E[Y|D=1, X]P(D=1|X)+ (ybar-E[Y|D=0,X])*P(D=0|X))
up_pos<-function(x){ mean(jtpa$edu==x)*(eydx(1,x)*pdx(x) + (ybar - eydx(0,x))*(1-pdx(x)))} #used when delta=1, delta*=0

## -P(X)*((E[Y|D=1, X]-ybar)P(D=1|X)+ (-E[Y|D=0,X])*P(D=0|X))
up_neg<- function(x){-mean(jtpa$edu==x)*((eydx(1,x)-ybar)*pdx(x)+(-eydx(0,x))*(1-pdx(x)))} #used when delta=0, delta*=1

## lower bound on welfare:
## P(X)*((E[Y|D=1, X]-ybar)P(D=1|X)+ (-E[Y|D=0,X])*P(D=0|X))
low_pos <- function(x){mean(jtpa$edu==x)*((eydx(1,x)-ybar)*pdx(x) + (- eydx(0,x))*(1-pdx(x)))} #used when delta=1, delta*=0

## -P(X)*(E[Y|D=1, X]P(D=1|X)+ (ybar-E[Y|D=0,X])*P(D=0|X))
low_neg <-function(x){-mean(jtpa$edu==x)*(eydx(1,x)*pdx(x)+(ybar-eydx(0,x))*(1-pdx(x)))} #used when delta=0, delta*=1

#################################################################################################
################################ (2) MTR  #######################################################
#################################################################################################
#benchmark rule (treat if delta(x)>0; not treat otherwise)
delta_mtr <- function(x){eyx(x)-eydx(0,x)*(1-pdx(x))} 

## upper bound on welfare:
## P(X)*(E[Y|D=1, X]P(D=1|X)+ (ybar-E[Y|D=0,X])*P(D=0|X))
up_pos_mtr<-function(x){ mean(jtpa$edu==x)*(eydx(1,x)*pdx(x) + (ybar - eydx(0,x))*(1-pdx(x)))} #used when delta=1, delta*=0

## lower bound on welfare:
## -P(X)*(E[Y|D=1, X]P(D=1|X)+ (ybar-E[Y|D=0,X])*P(D=0|X))
low_neg_mtr <-function(x){-mean(jtpa$edu==x)*(eydx(1,x)*pdx(x)+(ybar-eydx(0,x))*(1-pdx(x)))} #used when delta=0, delta*=1

#################################################################################################
################################ (3) IV-worst-case  #############################################
#################################################################################################
#benchmark rule:
delta_iv <- function(x){max(eydxz(1,x,1)*pdxz(x,1),eydxz(1,x,0)*pdxz(x,0))-max(eydxz(0,x,1)*(1-pdxz(x,1)), eydxz(0,x,0)*(1-pdxz(x,0)))}

## upper bound on welfare:
up_pos_iv <-function(x){ mean(jtpa$edu==x)*min((eydxz(1,x,1)*pdxz(x,1) + (ybar - eydxz(0,x,1))*(1-pdxz(x,1))), (eydxz(1,x,0)*pdxz(x,0) + (ybar - eydxz(0,x,0))*(1-pdxz(x,0))))} #used when delta=1, delta*=0
up_neg_iv <- function(x){-mean(jtpa$edu==x)*max(((eydxz(1,x,1)-ybar)*pdxz(x,1)+(-eydxz(0,x,1))*(1-pdxz(x,1))),((eydxz(1,x,0)-ybar)*pdxz(x,0)+(-eydxz(0,x,0))*(1-pdxz(x,0))))} #used when delta=0, delta*=1

## lower bound on welfare:
low_pos_iv <- function(x){mean(jtpa$edu==x)*max(((eydxz(1,x,1)-ybar)*pdxz(x,1) + (- eydxz(0,x,1))*(1-pdxz(x,1))),((eydxz(1,x,0)-ybar)*pdxz(x,0) + (- eydxz(0,x,0))*(1-pdxz(x,0))))} #used when delta=1, delta*=0
low_neg_iv <-function(x){-mean(jtpa$edu==x)*min((eydxz(1,x,1)*pdxz(x,1)+(ybar-eydxz(0,x,1))*(1-pdxz(x,1))),(eydxz(1,x,0)*pdxz(x,0)+(ybar-eydxz(0,x,0))*(1-pdxz(x,0))))} #used when delta=0, delta*=1

#################################################################################################
################################ (4) IV-MTR  ####################################################
#################################################################################################
#benchmark rule:
delta_iv_mtr <- function(x){max(eyxz(x,1),eyxz(x,0))-max(eydxz(0,x,1)*(1-pdxz(x,1)), eydxz(0,x,0)*(1-pdxz(x,0)))}

## upper bound on welfare:
up_pos_iv_mtr <-function(x){ mean(jtpa$edu==x)*min((eydxz(1,x,1)*pdxz(x,1) + (ybar - eydxz(0,x,1))*(1-pdxz(x,1))), (eydxz(1,x,0)*pdxz(x,0) + (ybar - eydxz(0,x,0))*(1-pdxz(x,0))))} #used when delta=1, delta*=0

## lower bound on welfare:
low_neg_iv_mtr <-function(x){-mean(jtpa$edu==x)*min((eydxz(1,x,1)*pdxz(x,1)+(ybar-eydxz(0,x,1))*(1-pdxz(x,1))),(eydxz(1,x,0)*pdxz(x,0)+(ybar-eydxz(0,x,0))*(1-pdxz(x,0))))} #used when delta=0, delta*=1

#################################################################################################
################################ (5) MIV-worst-case  ############################################
#################################################################################################
## upper bound on welfare:
up_pos_miv <- function(x){mean(jtpa$edu==x)*
    (pz*(eydxz(1,x,1)*pdxz(x,1)+ybar*(1-pdxz(x,1))-max(eydxz(0,x,1)*(1-pdxz(x,1)),eydxz(0,x,0)*(1-pdxz(x,0))))
  +(1-pz)*(min(eydxz(1,x,1)*pdxz(x,1)+ybar*(1-pdxz(x,1)),eydxz(1,x,0)*pdxz(x,0)+ybar*(1-pdxz(x,0)))-eydxz(0,x,0)*(1-pdxz(x,0))))}

up_neg_miv <- function(x){-mean(jtpa$edu==x)*
    (pz*(max(eydxz(1,x,1)*pdxz(x,1), eydxz(1,x,0)*pdxz(x,0))-(ybar*pdxz(x,1)+eydxz(0,x,1)*pdxz(x,1)))
     +(1-pz)*(eydxz(1,x,0)*pdxz(x,0))-min(ybar*pdxz(x,1)+eydxz(0,x,1)*pdxz(x,1),ybar*pdxz(x,0)+eydxz(0,x,0)*pdxz(x,0)))}

## lower bound on welfare:
low_pos_miv <- function(x){mean(jtpa$edu==x)*
    (pz*(max(eydxz(1,x,1)*pdxz(x,1),eydxz(1,x,0)*pdxz(x,0))-(ybar*pdxz(x,1)+eydxz(0,x,1)*(1-pdxz(x,1))))
     +(1-pz)*(eydxz(1,x,0)*pdxz(x,0)-min(ybar*pdxz(x,1)+eydxz(0,x,1)*(1-pdxz(x,1)),ybar*pdxz(x,0)+eydxz(0,x,0)*(1-pdxz(x,0)))))}


low_neg_miv <- function(x){-mean(jtpa$edu==x)*
      (pz*(eydxz(1,x,1)*pdxz(x,1)+ybar*(1-pdxz(x,1))-max(eydxz(0,x,1)*(1-pdxz(x,1)),eydxz(0,x,0)*(1-pdxz(x,0))))
       +(1-pz)*(min(eydxz(1,x,1)*pdxz(x,1)+ybar*(1-pdxz(x,1)),eydxz(1,x,0)*pdxz(x,0)+ybar*(1-pdxz(x,0)))-eydxz(0,x,0)*(1-pdxz(x,0))))}

#################################################################################################
################################ (5) MIV-MTR ############################################
#################################################################################################
## upper bound on welfare:
up_pos_miv_mtr <- function(x){mean(jtpa$edu==x)*
    (pz*(eydxz(1,x,1)*pdxz(x,1)+ybar*(1-pdxz(x,1))-max(eydxz(0,x,1)*(1-pdxz(x,1)),eydxz(0,x,0)*(1-pdxz(x,0))))
     +(1-pz)*(min(eydxz(1,x,1)*pdxz(x,1)+ybar*(1-pdxz(x,1)),eydxz(1,x,0)*pdxz(x,0)+ybar*(1-pdxz(x,0)))-eydxz(0,x,0)*(1-pdxz(x,0))))}

up_neg_miv_mtr <- function(x){-mean(jtpa$edu==x)*
    (pz*(max(eyxz(x,1),eyxz(x,0))-eyxz(x,1))
     +(1-pz)*(eyxz(x,0)-min(eyxz(x,1),eyxz(x,0))))}

## lower bound on welfare:
low_pos_miv_mtr <- function(x){mean(jtpa$edu==x)*(pz*(max(eyxz(x,1),eyxz(x,0))-eyxz(x,1))+(1-pz)*(eyxz(x,0)-min(eyxz(x,1),eyxz(x,0))))}

low_neg_miv_mtr <- function(x){-mean(jtpa$edu==x)*
    (pz*(eydxz(1,x,1)*pdxz(x,1)+ybar*(1-pdxz(x,1))-max(eydxz(0,x,1)*(1-pdxz(x,1)),eydxz(0,x,0)*(1-pdxz(x,0))))
     +(1-pz)*(min(eydxz(1,x,1)*pdxz(x,1)+ybar*(1-pdxz(x,1)),eydxz(1,x,0)*pdxz(x,0)+ybar*(1-pdxz(x,0)))-eydxz(0,x,0)*(1-pdxz(x,0))))}

#################################################################################################
################################ estimation results  ############################################
#################################################################################################
## (1) benchmark: treat if x leq 12; treatment: treat if x leq 16
#welfare_low1 <-  low_pos(13)+low_pos(14)+low_pos(15)+low_pos(16)
#welfare_up1 <-  up_pos(13)+up_pos(14)+up_pos(15)+up_pos(16)

#welfare_low2 <-  0
#welfare_up2 <- up_pos_mtr(13)+up_pos_mtr(14)+up_pos_mtr(15)+up_pos_mtr(16)

#welfare_low3 <-  low_pos_iv(13)+low_pos_iv(14)+low_pos_iv(15)+low_pos_iv(16)
#welfare_up3 <-  up_pos_iv(13)+up_pos_iv(14)+up_pos_iv(15)+up_pos_iv(16)

#welfare_low4 <- 0 
#welfare_up4 <-  up_pos_iv_mtr(13)+up_pos_iv_mtr(14)+up_pos_iv_mtr(15)+up_pos_iv_mtr(16)

#welfare_low5 <-  low_pos_miv(13)+low_pos_miv(14)+low_pos_miv(15)+low_pos_miv(16)
#welfare_up5 <-  up_pos_miv(13)+up_pos_miv(14)+up_pos_miv(15)+up_pos_miv(16)

#welfare_low6 <-  low_pos_miv_mtr(13)+low_pos_miv_mtr(14)+low_pos_miv_mtr(15)+low_pos_miv_mtr(16)
#welfare_up6 <-  up_pos_miv_mtr(13)+up_pos_miv_mtr(14)+up_pos_miv_mtr(15)+up_pos_miv_mtr(16)

#welfare <- rbind(welfare_low1, welfare_up1, welfare_low2, welfare_up2, welfare_low3, welfare_up3, welfare_low4, welfare_up4, welfare_low5, welfare_up5, welfare_low6, welfare_up6)


## (2) benchmark: treat if x leq 12; treatment: treat if x geq 12
#welfare_low1 <-  low_neg(7)+low_neg(8)+low_neg(9)+low_neg(10) + low_neg(11) + low_neg(12) +low_pos(13)+low_pos(14)+low_pos(15)+low_pos(16)+low_pos(17)+low_pos(18)
#welfare_up1 <-  up_neg(7)+up_neg(8)+up_neg(9)+up_neg(10) + up_neg(11) + up_neg(12) +up_pos(13)+up_pos(14)+up_pos(15)+up_pos(16)+up_pos(17)+up_pos(18)

#welfare_low2 <-  low_neg_mtr(7)+low_neg_mtr(8)+low_neg_mtr(9)+low_neg_mtr(10) + low_neg_mtr(11) + low_neg_mtr(12)
#welfare_up2 <- up_pos_mtr(13)+up_pos_mtr(14)+up_pos_mtr(15)+up_pos_mtr(16)+up_pos_mtr(17)+up_pos_mtr(18)

#welfare_low3 <-  low_neg_iv(7)+low_neg_iv(8)+low_neg_iv(9)+low_neg_iv(10) + low_neg_iv(11) + low_neg_iv(12) +low_pos_iv(13)+low_pos_iv(14)+low_pos_iv(15)+low_pos_iv(16)+low_pos_iv(17)+low_pos_iv(18)

#welfare_up3 <-  up_neg_iv(7)+up_neg_iv(8)+up_neg_iv(9)+up_neg_iv(10) + up_neg_iv(11) + up_neg_iv(12) +up_pos_iv(13)+up_pos_iv(14)+up_pos_iv(15)+up_pos_iv(16)+up_pos_iv(17)+up_pos_iv(18)

#welfare_low4 <- low_neg_iv_mtr(7)+low_neg_iv_mtr(8)+low_neg_iv_mtr(9)+low_neg_iv_mtr(10) + low_neg_iv_mtr(11) + low_neg_iv_mtr(12)
#welfare_up4 <- up_pos_iv_mtr(13)+up_pos_iv_mtr(14)+up_pos_iv_mtr(15)+up_pos_iv_mtr(16)+up_pos_iv_mtr(17)+up_pos_iv_mtr(18)

#welfare_low5 <-  low_neg_miv(7)+low_neg_miv(8)+low_neg_miv(9)+low_neg_miv(10) + low_neg_miv(11) + low_neg_miv(12) +low_pos_miv(13)+low_pos_miv(14)+low_pos_miv(15)+low_pos_miv(16)+low_pos_miv(17)+low_pos_miv(18)
#welfare_up5 <- up_neg_miv(7)+up_neg_miv(8)+up_neg_miv(9)+up_neg_miv(10) + up_neg_miv(11) + up_neg_miv(12)+up_pos_miv(13)+up_pos_miv(14)+up_pos_miv(15)+up_pos_miv(16)+up_pos_miv(17)+up_pos_miv(18)

#welfare_low6 <-  low_neg_miv_mtr(7)+low_neg_miv_mtr(8)+low_neg_miv_mtr(9)+low_neg_miv_mtr(10) + low_neg_miv_mtr(11) + low_neg_miv_mtr(12)+low_pos_miv_mtr(13)+low_pos_miv_mtr(14)+low_pos_miv_mtr(15)+low_pos_miv_mtr(16)+low_pos_miv_mtr(17)+low_pos_miv_mtr(18)
#welfare_up6 <-  up_neg_miv_mtr(7)+up_neg_miv_mtr(8)+up_neg_miv_mtr(9)+up_neg_miv_mtr(10) + up_neg_miv_mtr(11) + up_neg_miv_mtr(12)+up_pos_miv_mtr(13)+up_pos_miv_mtr(14)+up_pos_miv_mtr(15)+up_pos_miv_mtr(16)+up_pos_miv_mtr(17)+up_pos_miv_mtr(18)

#welfare <- rbind(welfare_low1, welfare_up1, welfare_low2, welfare_up2, welfare_low3, welfare_up3, welfare_low4, welfare_up4, welfare_low5, welfare_up5, welfare_low6, welfare_up6)


## (3) benchmark: treat if x leq 11; treatment: treat if x leq 12
welfare_low1 <-  low_pos(12)
welfare_up1 <-  up_pos(12)

welfare_low2 <-  0
welfare_up2 <- up_pos_mtr(12)

welfare_low3 <-  low_pos_iv(12)
welfare_up3 <-  up_pos_iv(12)

welfare_low4 <- 0 
welfare_up4 <-  up_pos_iv_mtr(12)

welfare_low5 <-  low_pos_miv(12)
welfare_up5 <-  up_pos_miv(12)

welfare_low6 <-  low_pos_miv_mtr(12)
welfare_up6 <-  up_pos_miv_mtr(12)

welfare <- rbind(welfare_low1, welfare_up1, welfare_low2, welfare_up2, welfare_low3, welfare_up3, welfare_low4, welfare_up4, welfare_low5, welfare_up5, welfare_low6, welfare_up6)


