# Function to generate data from the data generating process specified in 
# Identification and Inference for Welfare Gains without Unconfoundedness
# Written by: Undral Byambadalai, Boston University 

gen_data <- function(n=10){
s <- 11000 #variance of Y
X <- sample(x=7:18, prob = c(0.01, 0.06, 0.07, 0.11, 0.13, 0.43, 0.07, 0.06, 0.02, 0.02, 0.01, 0.01), size = n, replace = TRUE)
Z <- sample(x=0:1, prob=c(1/3,2/3), size = n, replace = TRUE)
U <- runif(n,0,1)
D <- vector()
for (i in 1:n){if(U[i]<1/(1+exp(-(-4.89+0.05*X[i]+c*Z[i])))){D[i] <- 1} else {D[i] <-0}}
m1 <- vector()
s1 <- vector()
for (i in 1:n){m1[i] <- log((5591+1027*X[i]+gamma1*U[i])^2/sqrt(s^2+(5591+1027*X[i]+gamma1*U[i])^2))}
for (i in 1:n){s1[i] <- sqrt(log(1+(s^2/(5591+1027*X[i]+gamma1*U[i])^2)))}
m0 <- vector()
s0 <- vector()
for (i in 1:n){m0[i] <- log((-1127+1389*X[i]+gamma0*U[i])^2/sqrt(s^2+(-1127+1389*X[i]+gamma0*U[i])^2))}
for (i in 1:n){s0[i] <- sqrt(log(1+(s^2/(-1127+1389*X[i]+gamma0*U[i])^2)))}
Y1 <- vector()
Y0 <- vector()
for (i in 1:n) {Y1[i] <- rlnorm(n=1, mean=m1[i], sd=s1[i])}
for (i in 1:n) {Y0[i] <- rlnorm(n=1, mean=m0[i], sd=s0[i])}
jtpa <- data.frame(Y1,Y0,D,Z,X)
jtpa$earnings <- jtpa$Y1*jtpa$D + jtpa$Y0*(1-jtpa$D)
jtpa$edu <- jtpa$X
return (jtpa)
}
