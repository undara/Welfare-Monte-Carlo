# Function to calculate population quantities given the data generating process in
# Identification and Inference for Welfare Gains without Unconfoundedness
# Written by: Undral Byambadalai, Boston University 

population_wc <- function(c=10){
# Setup  
y_under <- 0
y_bar <- 160000
gamma1 <-2000
gamma0 <-1000
# propensity score, mean treatment response
p.z <- function(x){2/3}
p.xz <-  function (x,z){1/(1+exp(-(-4.89+0.05*x +c*z)))} 
p.x <- function(x){p.z(x)*p.xz(x,1)+(1-p.z(x))*p.xz(x,0)}
m1.xu <- function(x,u){5591+ 1027*x + gamma1*u}
m0.xu <- function(x,u){-1127 + 1389*x + gamma0*u}
m1 <- function(u){m1.xu(12,u)}
m0 <- function(u){m0.xu(12,u)}
f <- function(u){m1.xu(12,u)-m0.xu(12,u)}
#welfare gain
beta <- round(0.43*as.numeric(integrate(f,0,1)[1]), digits=0)
# E[Y|D=0, X,Z=1]
y.d0xz1 <- round(1/(1-p.xz(12,1))*as.numeric(integrate(m0, p.xz(12,1),1)[1]), digits =0)
#E [Y|D=0,X, Z=0]
y.d0xz0 <- round(1/(1-p.xz(12,0))*as.numeric(integrate(m0, p.xz(12,0),1)[1]), digits =0)
# E[Y|D=1, X,Z=1]
y.d1xz1 <- round(1/p.xz(12,1)*as.numeric(integrate(m1, 0,p.xz(12,1))[1]), digits =0)
#E [Y|D=1,X, Z=0]
y.d1xz0 <- round(1/p.xz(12,0)*as.numeric(integrate(m1, 0, p.xz(12,0))[1]), digits =0)
#p(Z=1|D=0,X)
pz1.d0x <- round((2/3*(1-p.xz(12,1)))/(1-2/3*p.xz(12,1)-1/3*p.xz(12,0)), digits = 2)
#p(Z=0|D=0,X)
pz0.d0x <- round((1/3*(1-p.xz(12,0)))/(1-2/3*p.xz(12,1)-1/3*p.xz(12,0)), digits = 2)
#p(Z=1|D=1,X)
pz1.d1x <- round((2/3*p.xz(12,1))/(2/3*p.xz(12,1)+1/3*p.xz(12,0)), digits = 2)
#p(Z=0|D=1,X)
pz0.d1x <- round((1/3*p.xz(12,0))/(2/3*p.xz(12,1)+1/3*p.xz(12,0)), digits = 2)
# E[Y|D=0,X]
y.d0x <- round(y.d0xz1*pz1.d0x + y.d0xz0*pz0.d0x, digits=0)
# E[Y|D=1,X]
y.d1x <- round(y.d1xz1*pz1.d1x + y.d1xz0*pz0.d1x, digits=0)
#P(D=1|X)
pd1.x <- round(p.x(12), digits=2)
# worst-case:
wc_low  <- 0.43*((y.d1x-y_bar)*pd1.x + (y_under-y.d0x)*(1-pd1.x))
wc_up <- 0.43*((y.d1x-y_under)*pd1.x + (y_bar-y.d0x)*(1-pd1.x))
# IV-worst case:
iv1 <- max(y.d1xz1*p.xz(12,1) + y_under*(1-p.xz(12,1)),y.d1xz0*p.xz(12,0) + y_under*(1-p.xz(12,0)))
iv2 <- min(y_bar*p.xz(12,1) + y.d0xz1*(1-p.xz(12,1)),y_bar*p.xz(12,0) + y.d0xz0*(1-p.xz(12,0)))
iv_wc_low <- round(0.43*(iv1 - iv2), digits=0)
list <- list(eta0=y.d0x, eta1=y.d1x, p=pd1.x, beta_0=wc_low, 
             eyd0xz1=y.d0xz1, eyd0xz0=y.d0xz0, eyd1xz1=y.d1xz1, eyd1xz0=y.d1xz0, 
             pdx1 = p.xz(12,1), pdx0 =p.xz(12,0),
             iv_max =iv1, iv_min=iv2, beta_iv = iv_wc_low)
return (list)
}
