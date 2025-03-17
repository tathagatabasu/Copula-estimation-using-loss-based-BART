# code for copula estimation

install.packages(c("copula", "MASS", "coda"))
library(copula)
library(MASS)   # For multivariate normal functions
library(coda)   # For MCMC diagnostics

# data generation

set.seed(123)

# Define true copula parameter
rho_true <- 0.7

# Define the Gaussian copula
cop <- normalCopula(param = rho_true, dim = 2)

# Generate uniform data from the copula
n <- 500
u <- rCopula(n, cop)

# Convert to normal marginals
x1 <- qnorm(u[,1])
x2 <- qnorm(u[,2])

# Store in a dataframe
data <- data.frame(x1, x2)

# define likelihood function

log_likelihood <- function(rho, x1, x2) {
  if (abs(rho) >= 1) return(-Inf)  # Ensure rho is in valid range (-1,1)
  
  # Compute copula density
  n <- length(x1)
  log_lik <- -0.5 * sum(
    (x1^2 - 2 * rho * x1 * x2 + x2^2) / (1 - rho^2)
  ) - n * log(sqrt(1 - rho^2))
  
  return(log_lik)
}

# MH algorithm

metropolis_hastings <- function(log_lik_func, x1, x2, 
                                n_iter = 10000, init_rho = 0.5, proposal_sd = 0.05) {
  rho <- numeric(n_iter)
  rho[1] <- init_rho
  accept <- 0
  
  for (i in 2:n_iter) {
    # Proposal step
    rho_prop <- rnorm(1, mean = rho[i - 1], sd = proposal_sd)
    
    # Compute acceptance ratio
    log_alpha <- log_lik_func(rho_prop, x1, x2) - log_lik_func(rho[i - 1], x1, x2)
    alpha <- exp(log_alpha)
    
    # Accept or reject
    if (runif(1) < alpha) {
      rho[i] <- rho_prop
      accept <- accept + 1
    } else {
      rho[i] <- rho[i - 1]
    }
  }
  
  cat("Acceptance Rate:", accept / n_iter, "\n")
  return(rho)
}

# Run MCMC
set.seed(42)
posterior_samples <- metropolis_hastings(log_likelihood, x1, x2, n_iter = 5000)

# illustration

# Convert to MCMC object
mcmc_samples <- as.mcmc(posterior_samples)

# Plot trace and density
plot(mcmc_samples, col = "blue")

summary(mcmc_samples)
