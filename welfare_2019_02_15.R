## Title: Identification of welfare differences of policies
## Written by: Undral Byambadalai
## Date: 02/19/2019

## this code calculates 
## benchmark treatment rule
## upper & lower bounds on welfare regret
## under Manski bound
## when X=edu is used as a covariate

################################ read data ##############################################################################
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

#########################################################################################################################

## E[Y|D=d, X=x]
yhat <- function(d,x){ mean(jtpa$earnings[jtpa$D==d & jtpa$edu==x])} 

## P[D=1|X=x]=E[D|X=x]
phat <- function(x){mean(jtpa$D[jtpa$edu==x])} 

#benchmark rule (treat if delta(x)>0; not treat otherwise)
delta <- function(x){yhat(1,x)*phat(x)-yhat(0,x)*(1-phat(x))} 

## set ybar (upper bound on outcome)
ybar <- 160000

## upper bound on welfare:
## P(X)*(E[Y|D=1, X]P(D=1|X)+ (ybar-E[Y|D=0,X])*P(D=0|X))
up_pos<-function(x){ mean(jtpa$edu==x)*(yhat(1,x)*phat(x) + (ybar - yhat(0,x))*(1-phat(x)))} #used when delta*=1, delta=0

## -P(X)*((E[Y|D=1, X]-ybar)P(D=1|X)+ (-E[Y|D=0,X])*P(D=0|X))
up_neg<- function(x){-mean(jtpa$edu==x)*((yhat(1,x)-ybar)*phat(x)+(-yhat(0,x))*(1-phat(x)))} #used when delta*=0, delta=1

## lower bound on welfare:
## P(X)*((E[Y|D=1, X]-ybar)P(D=1|X)+ (-E[Y|D=0,X])*P(D=0|X))
low_pos <- function(x){mean(jtpa$edu==x)*((yhat(1,x)-ybar)*phat(x) + (- yhat(0,x))*(1-phat(x)))} #used when delta*=1, delta=0

## -P(X)*(E[Y|D=1, X]P(D=1|X)+ (ybar-E[Y|D=0,X])*P(D=0|X))
low_neg <-function(x){-mean(jtpa$edu==x)*(yhat(1,x)*phat(x)+(ybar-yhat(0,x))*(1-phat(x)))} #used when delta*=0, delta=1


################################# various treatment rules ################################################################

### benchmark: treatment if x=15, not treat otherwise
# treatment: treat if more than 12
up_neg(12)+up_neg(13)+up_neg(14)+up_neg(16)+up_neg(17)+up_neg(18)
low_neg(12)+low_neg(13)+low_neg(14)+low_neg(16)+low_neg(17)+low_neg(18)

# treatment: treat if more than 13
up_neg(13)+up_neg(14)+up_neg(16)+up_neg(17)+up_neg(18)
low_neg(13)+low_neg(14)+low_neg(16)+low_neg(17)+low_neg(18)

# treatment: treat if more than 14
up_neg(14)+up_neg(16)+up_neg(17)+up_neg(18)
low_neg(14)+low_neg(16)+low_neg(17)+low_neg(18)

# treatment: treat if more than 15
up_neg(16)+up_neg(17)+up_neg(18)
low_neg(16)+low_neg(17)+low_neg(18)

# treatment: treat if more than 16
up_pos(15)+up_neg(16)+up_neg(17)+up_neg(18)
low_pos(15)+low_neg(16)+low_neg(17)+low_neg(18)

# treatment: treat if equal to 12
up_neg(12)+up_pos(15)
low_neg(12)+low_pos(15)

# treatment: treat everyone
up_neg(7)+up_neg(8)+up_neg(9)+up_neg(10)+up_neg(11)+up_neg(12)+up_neg(13)+up_neg(14)+up_neg(16)+up_neg(17)+up_neg(18)
low_neg(7)+low_neg(8)+low_neg(9)+low_neg(10)+low_neg(11)+low_neg(12)+low_neg(13)+low_neg(14)+low_neg(16)+low_neg(17)+low_neg(18)

# treatment: treat noone
up_pos(15)
low_pos(15)

## benchmark: treat if x leq 12 (until high school)
## treatment: treat if x leq 16 (until college)
up_neg(13)+up_neg(14)+up_neg(15)+up_neg(16)
low_neg(13)+low_neg(14)+low_neg(15)+low_neg(16)
