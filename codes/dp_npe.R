library(distr) 
# https://www.di.fc.ul.pt/~jpn/r/dp/dp.html

# emp_cdf is an object of class "DiscreteDistribution"
cdf_sample <- function(emp_cdf, n=1e3) {
  emp_cdf@r(n)
}

set.seed(101)
# make empirical cdf from some samples
X       <- rnorm(1e3)
emp_cdf <- DiscreteDistribution(X)

# get some points from this empirical cdf
someXs <- cdf_sample(emp_cdf, 1e3)

# check if it's near the original sample generator (norm(0,1))
xs <- seq(-3,3,len=100)
plot(xs, pnorm(xs), type="l", lwd=2)
points(xs, ecdf(someXs)(xs), type="l", col="red")

# Dirichlet process to draw a random distribution F from a prior F0
#  alpha determines the similarity with the prior guess
#  F0 is the cdf of the prior guess (an object of class "DiscreteDistribution")
dp <- function(alpha, F0, n=1e3) { # n should be large since it's an approx for +oo
  
  s <- cdf_sample(F0,n)            # step 1: draw from F0
  V <- rbeta(n,1,alpha)            # step 2: draw from beta(1,alpha)
  w <- c(1, rep(NA,n-1))           # step 3: compute 'stick breaking process' into w[i]
  w[2:n] <- sapply(2:n, function(i) V[i] * prod(1 - V[1:(i-1)]))
  
  # return the sampled function F which can be itself sampled 
  # this F is a probability mass function where each s[i] has mass w[i]
  function (size=1e4) {
    sample(s, size, prob=w, replace=TRUE)
  }
}

f0 <- function(n) rnorm(n, 0, 1)    # eg pdf of prior guess
F0 <- DiscreteDistribution(f0(1e4)) # make its cdf

# generate a prior from the Dirichlet process
dpF <- dp(10, F0, n=1e4)

# plot the pdf
hist(dpF(), breaks=50, prob=T, ylab="", xlab="",
     main=expression(paste("pmf of F ~ ",pi)))

# plot the cdf
plot(xs, pnorm(xs, 0, 1), type="l", col="blue", lwd=2, ylab="", xlab="",
     main=expression(paste("cdf of F ~ ",pi)))
points(xs, ecdf(dpF())(xs), type="l", col="red")
legend("topleft",c("prior guess","sampled F"), 
       col=c("blue","red"), lwd=c(2,1), bty = "n")

# X is the evidence (x1, x2, ..., xn)
dp_posterior <- function(alpha, F0, X) {
  n <- length(X)
  F_n <- DiscreteDistribution(X) # compute empirical cdf
  
  F_bar <- n/(n+alpha) * F_n + alpha/(n+alpha) * F0
  
  dp(alpha+n, F_bar)
}

f0 <- function(n) rnorm(n, 1, 1)    # the prior guess (in pdf format)
F0 <- DiscreteDistribution(f0(1e3)) # the prior guess (in cdf format)

data <- rnorm(30,3,1)               # the data

# apply Dirichlet process
runs  <- 50
xs    <- seq(-2,6,len=50)
y_hat <- matrix(nrow=length(xs), ncol=runs)
for(i in 1:runs) {
  Fpost <- dp_posterior(10, F0, data)  
  y_hat[,i] <- ecdf(Fpost())(xs)    # just keeping values of cdf(xs), not the object
}

df <- 2
f0 <- function(n) rt(n, df=df, ncp=1) # the prior guess (in pdf format)
F0 <- DiscreteDistribution(f0(1e3))   # the prior guess (in cdf format)

y_hat <- matrix(nrow=length(xs), ncol=runs)
for(i in 1:runs) {
  Fpost <- dp_posterior(10, F0, data)  
  y_hat[,i] <- ecdf(Fpost())(xs)    # just keeping values of cdf(xs), not the object
}

# plot the black area
plot(xs, pt(xs, df=df, ncp=1), ylim=c(-.1,1.1), type="n", col="blue", lwd=2, ylab="", xlab="")

# compute & plot 95% credible interval of the posterior
crebible_int <- apply(y_hat, 1, function(row) HPDinterval(as.mcmc(row), prob=0.95))
polygon(c(rev(xs), xs), c(rev(crebible_int[1,]), 
                          crebible_int[2,]), col = 'grey90')    

# plot the prior cdf
points(xs, pt(xs, df=df, ncp=1), type="l", col="blue", lwd=2)

# plot mean estimate of the posterior
means <- apply(y_hat, 1, mean)
points(xs, means, type="l", col="red", lwd=2)                  

# plot true data generator
points(xs, pnorm(xs, 3, 1), type="l", col="darkgreen", lwd=2)
legend("topleft",c("prior","posterior mean", "truth"), 
       col=c("blue","red","darkgreen"), lwd=2, bty = "n") 

