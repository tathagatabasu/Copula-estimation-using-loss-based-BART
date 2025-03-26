# codes and packages

#install.packages(c("copula", "MASS", "coda"))
source('code/import_functions.R')
source('mclapply.R')
library(data.tree)
library(ggplot2)
library(ggpubr)
library(copula) # VineCopula
library(MASS)   # For multivariate normal functions
library(coda)   # For MCMC diagnostics

# data generation

set.seed(123)

# generate predictor
n <- 500
x_pred <- matrix(runif(2*n), ncol = 2)

# Define true copula parameter
calib_true <- 0.6 * sin(x_pred^2 %*% c(3,5))

rho_true <- 2/(abs(calib_true)+1) - 1

# Define the Gaussian copula
cop <- sapply(rho_true, function(x)BiCopSim(1, 1, x))

# Convert to normal marginals
y1 <- qnorm(cop[1,])
y2 <- qnorm(cop[2,])

# empirical cdf function
ecdf_y1 = ecdf(y1)
ecdf_y2 = ecdf(y2)

##################################################
pseudo_u1 = sapply(y1, ecdf_y1)
pseudo_u2 = sapply(y2, ecdf_y2)

triweight <- function(u){
  35/32*(1-sum(u^2))^3 *(sum(u^2)<1)
}

epanechnikov <- function(u){
  3/4*(1-sum(u^2)) *(sum(u^2)<1)
}

gaussweight <- function(x){
  exp(-sum(x^2)/2)/(sqrt(2*pi))^length(x)
}

NW_weights <- function(x, x_obs, band = 2.3 * (1/nrow(x_obs))^.2){
  wt <- apply(x_obs, 1, function(t)triweight((x-t)/band))
  
  return(wt/sum(wt))
}

rho.cond=function(uu, x, x_obs)
{
  uu1=uu[,1]
  uu2=uu[,2]
  estim=c()
  
  n <- nrow(uu)
  ww <- NW_weights(x, x_obs)
  
  # Compute uu.hat
  uu1.hat=c()
  uu2.hat=c()
  
  
  for(i in 1:n)
  {
    uu1.hat[i]=sum(ww*(uu1<=uu1[i]))
    uu2.hat[i]=sum(ww*(uu2<=uu2[i]))
  }
  estim=12*sum(ww*(1-uu1.hat)*(1-uu2.hat))-3
  # Return
  return(estim)
}


sample_rho = apply(x_pred, 1, function(x)rho.cond(cbind(pseudo_u1,pseudo_u2), x, x_pred))

plot(rho_true, sin(sample_rho * pi/6))

######################

X_pred <- x_pred

summary(X_pred)

Y_mal <- 1/2 *log((1+sample_rho)/(1-sample_rho))

plot(calib_true, Y_mal)

############################################################
#
# normalise predictors 
X_pred.norm <- as.data.frame(apply(X_pred, 2, \(x) (x - min(x))/(max(x) - min(x))))
X_pred.norm <- as.matrix(X_pred.norm)
rownames(X_pred.norm) <- 1:nrow(X_pred.norm)
# calculate correlations

Y_mal <- (Y_mal - max(Y_mal))/(max(Y_mal)-min(Y_mal))

cor(X_pred.norm, Y_mal)




n.chain_par <- 5
n.iter_par <- 500
incl.split_par <- FALSE
cont.unif_par <- TRUE
moves.prob_par <- c(0.3, 0.3, 0.3, 0.1)

#############
## DEFUALT ##
#############

lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.56, 0.62))
mcmc_lb.def <- multichain_MCMC_known_var(n.chain = n.chain_par,
                                         n.iter = n.iter_par,
                                         X = X_pred.norm,
                                         Y = Y_mal,
                                         Y.var = 1, 
                                         n.cores = 5,
                                                     mu.prior.mean = 0, 
                                         mu.prior.var = var(Y_mal), 
                                                     prior_list = lb.prior.def, 
                                                     moves.prob = moves.prob_par, 
                                         starting.tree = NULL,
                                                     cont.unif = cont.unif_par,
                                                     include.split = incl.split_par)
####################
## DEFAULT MODELS ##
####################

model.list.def <- list(
  mcmc_lb.def)

names(model.list.def) <- c(
  'LB - default')


# extract depth, number of terminal nodes, missing rate and loglik of all the trees
depth.df <- apply_fun_models(fun_ = get_depth, 
                             mcmc.list = model.list.def,
                             born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)
nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes, 
                             mcmc.list = model.list.def, 
                             born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)

hist.nl <- ggplot(nterm.df) + 
  geom_histogram(aes(x = y, y = after_stat(density)), 
                 binwidth = 1, color = 'black', fill = 'white') + 
  facet_wrap(facets = ~panel.name) + 
  xlab(~n[L]) + 
  ylab('PMF') + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,24,by = 3))  

hist.depth <- ggplot(depth.df) + 
  geom_histogram(aes(x = y, y = after_stat(density)), 
                 binwidth = 1, color = 'black', fill = 'white') + 
  facet_wrap(facets = ~panel.name) + 
  xlab('Depth') + 
  ylab('PMF') + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,10,by = 1))

hist.nl
hist.depth

x_new <- matrix(runif(2*n, min = -1), ncol = 2)

calib_true_new <- x_new %*% c(.2,1)

rho_true_new <- (exp(2*calib_true_new)-1)/ (exp(2*calib_true_new)+1)

x_new_norm <- as.data.frame(sapply(1:ncol(x_new), function(i)(x_new[,i] - min(x_pred[,i]))/(max(x_pred[,i])- min(x_pred[,i]))))
X_new_norm <- as.matrix(x_new_norm)
rownames(x_new_norm) <- 1:nrow(X_pred.norm)


pred_cond = sapply(1:length(model.list.def$`LB - default`$trees), function(i)get_value_tree(model.list.def[[1]]$trees[[i]],x_new_norm))

est_rho_link = rowMeans(pred_cond)

est_rho = (exp(2*est_rho_link)-1)/ (exp(2*est_rho_link)+1)

plot(est_rho, rho_true_new)
