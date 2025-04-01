# codes and packages

#install.packages(c("copula", "MASS", "coda"))
source('code/import_functions.R')
source('mclapply.R')
library(data.tree)
library(ggplot2)
library(ggpubr)
library(CondCopulas)
library(VineCopula)
library(ggplot2)
library(MASS)   # For multivariate normal functions
library(coda)   # For MCMC diagnostics

# data generation

set.seed(123)

# generate predictor
n <- 500
X_pred <- matrix(runif(n), ncol = 1)

# Define true kendall's tau
tau_true <- 0.6 * sin(3*X_pred) + 0.3

# Simulation of Y1 and Y2
Y1 <- rep(NA, n)
Y2 <- rep(NA, n)

for (i in 1:n) {
  simCopula <- BiCopSim(N=1 , family = 1, par = BiCopTau2Par(1 , tau_true[i])) # gaussian
  Y1[i] <- qnorm(simCopula[1])
  Y2[i] <- qnorm(simCopula[2])
}

plot(Y1, Y2)

##################################################
# empirical cdf function
ecdf_y1 = ecdf(Y1)
ecdf_y2 = ecdf(Y2)

pseudo_u1 = sapply(Y1, ecdf_y1)
pseudo_u2 = sapply(Y2, ecdf_y2)

triweight <- function(u){
  35/32*(1-sum(u^2))^3 *(sum(u^2)<1)
}

NW_weights <- function(x, x_obs, band = 1.06 * min(sd(x_obs), IQR(x_obs)/1.34) * (1/nrow(x_obs))^.2){
  wt <- apply(x_obs, 1, function(t)triweight((x-t)/band))
  
  return(wt/sum(wt))
}

tau_cond <- function(y1_obs, y2_obs, x, x_obs)
{
  n <- length(y1_obs)
  ww <- NW_weights(x, x_obs)
  
  sum_cop <- 0
  
  for(i in 1:n){
    for(j in 1:n){
      sum_cop <- sum_cop + ww[i]*ww[j]* ((y1_obs[i] < y1_obs[j]) & (y2_obs[i] < y2_obs[j]))
    }
  }
  
  estim <- (4/(1 - sum(ww^2))) * sum_cop - 1
  
  # Return
  return(estim)
}

sample_tau <- apply(X_pred, 1, function(x)tau_cond(pseudo_u1, pseudo_u2, x, X_pred))

plot(X_pred, tau_true, col = "blue")
points(X_pred, sample_tau, col = "red")


######################

summary(X_pred)

Y_mal <- sample_tau



############################################################
#
# normalise predictors 
# X_pred.norm <- as.data.frame(apply(X_pred, 2, \(x) (x - min(x))/(max(x) - min(x))))
# X_pred.norm <- as.matrix(X_pred.norm)
rownames(X_pred) <- 1:nrow(X_pred)
# # calculate correlations
# 
# Y_mal <- (Y_mal - min(Y_mal))/(max(Y_mal)-min(Y_mal))

cor(X_pred, Y_mal)

n.chain_par <- 5
n.iter_par <- 500
incl.split_par <- TRUE
cont.unif_par <- TRUE
moves.prob_par <- c(0.1, 0.4, 0.25, 0.25)

#############
## DEFUALT ##
#############

lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944))
mcmc_lb.def <- multichain_MCMC_known_var(n.chain = n.chain_par,
                                         n.iter = n.iter_par,
                                         X = X_pred,
                                         Y = Y_mal,
                                         Y.var = 0.1, 
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

get_tree_plot(mcmc_lb.def[[1]][[1]])

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

trace.nl <- ggplot(nterm.df) + 
  geom_vline(xintercept = seq(0,1250,by = 50), 
             color = 'grey', size = 0.2, alpha=0.75)+
  geom_line(aes(x, y)) + 
  facet_wrap(facets = ~panel.name, ncol = 2) + 
  xlab('Iteration') + 
  ylab(~n[L]) + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,1250,by = 50)) + 
  theme(axis.text.x = element_text(angle = 30))

trace.depth <- ggplot(depth.df) + 
  geom_vline(xintercept = seq(0,1250,by = 50),
             color = 'grey', size = 0.2, alpha=0.75)+
  geom_line(aes(x, y)) +
  facet_wrap(facets = ~panel.name, ncol = 2) + 
  xlab('Iteration') + 
  ylab('Depth') + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,1250,by = 50)) + 
  theme(axis.text.x = element_text(angle = 30))

hist.nl
hist.depth

trace.nl
trace.depth

x_new <- matrix(runif(n), ncol = 1)

rownames(x_new) <- 1:nrow(x_new)

tau_true_new <- 0.6 * sin(3*x_new) + 0.3

pred_cond = sapply(1:length(model.list.def$`LB - default`$trees), function(i)get_value_tree(model.list.def[[1]]$trees[[i]],x_new))

est_tau = rowMeans(pred_cond)

plot(x_new, as.vector(tau_true_new), col = "red")
points(x_new, est_tau, col = "blue")

##################################################################################
# install.packages("GauPro")

library(GauPro)

gp_tau <- GauPro(X = X_pred, Z = sample_tau)

est_tau_gp <- predict(gp_tau, x_new)

plot(x_new, as.vector(tau_true_new), col = "red")
points(x_new, est_tau_gp, col = "blue")


library(rpart)

tree_dat <- data.frame("Y" = Y_mal, "X" = X_pred)

rtree.fit <- rpart(Y ~ X,
                   data = tree_dat,
                   method="anova", #for regression tree
                   control=rpart.control(minsplit=30,cp=0.001))
plot(rtree.fit, uniform=TRUE, 
     main="Regression Tree for Kendall's Tau")
text(rtree.fit, use.n=TRUE, all=TRUE, cex=.8)

test_tree_dat <- data.frame("X" = x_new)
tau_cart <- predict(rtree.fit, test_tree_dat)

plot(x_new, as.vector(tau_true_new), col = "red")
points(x_new, tau_cart, col = "blue")
