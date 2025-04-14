# codes and packages
#install.packages("mc2d")
#install.packages(c("copula", "MASS", "coda"))
source('code/import_functions.R')
source('mclapply.R')
source('MCMC_BART_copula.R')
library(data.tree)
library(dplyr)
library(ggplot2)
# library(ggpubr)
library(CondCopulas)
library(VineCopula)
library(latex2exp)
library(ggplot2)
library(MASS)   # For multivariate normal functions
library(coda)   # For MCMC diagnostics
library(plot3D)

# data generation

set.seed(123)

# generate predictor
n <- 500
X_obs <- matrix(runif(n), ncol = 1)

# Define gaussian copula parameter (rho)

# rho with tree structure
tree_top <- generate_random_binary_tree_n_delta(8,4)
tree_top <- assign_node_idx(tree_top)
tree_top <- assign_split_rules(tree_top, X_obs)
tree_top <- assign_term_node_values_binary(tree_top, 7, 2)
rho_true_1 <- sample_CART(tree_top, X_obs, sigma_ = 0.001) 

# monotone
rho_true_2 <- 0.5 + 0.2 * sin(3*X_obs) + 0.3*X_obs^2
# convex
rho_true_3 <- 0.6 + 0.3 * sin(3*X_obs)
# concave
rho_true_4 <- .9 - 0.3 * sin(3*X_obs)
# non-convex
rho_true_5 <- 0.7 - 0.3 * sin(2*X_obs) + 0.2 * sin(4*X_obs) + 0.3 * X_obs^2


plot(X_obs, rho_true_1)
plot(X_obs, rho_true_2)
plot(X_obs, rho_true_3)
plot(X_obs, rho_true_4)
plot(X_obs, rho_true_5)

for (i in 1:5) {
  assign(paste0("copula_uu_",i), sapply(1:n, function(k)BiCopSim(N=1 , family = 1, par = get(paste0("rho_true_",i))[k])))
}

plot(copula_uu_1[1,], copula_uu_1[2,], xlab = "U1", ylab = "U2")
plot(copula_uu_2[1,], copula_uu_2[2,], xlab = "U1", ylab = "U2")
plot(copula_uu_3[1,], copula_uu_3[2,], xlab = "U1", ylab = "U2")
plot(copula_uu_4[1,], copula_uu_4[2,], xlab = "U1", ylab = "U2")
plot(copula_uu_5[1,], copula_uu_5[2,], xlab = "U1", ylab = "U2")

##################################################
# normalise predictors 
X_obs.norm <- as.data.frame(apply(X_obs, 2, \(x) (x - min(x))/(max(x) - min(x))))
X_obs.norm <- as.matrix(X_obs.norm)
rownames(X_obs.norm) <- 1:nrow(X_obs)

n.chain_par <- 5
n.iter_par <- 500
incl.split_par <- TRUE
cont.unif_par <- TRUE
moves.prob_par <- c(0.1, 0.4, 0.25, 0.25)

#############
## DEFUALT ##
#############

lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) # c(1.5618883, 0.6293944)

for (i in 1:5) {
  assign(paste0("mcmc_lb.def_unif_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                          n.iter = n.iter_par,
                                                          X = X_obs.norm,
                                                          U1 = get(paste0("copula_uu_",i))[1,],
                                                          U2 = get(paste0("copula_uu_",i))[2,],
                                                          mu = 0, 
                                                          sigma = .1, alpha_val = 1, beta_val = 1, 
                                                          prior_list = lb.prior.def, 
                                                          moves.prob = moves.prob_par, 
                                                          starting.tree = NULL,
                                                          cont.unif = cont.unif_par,
                                                          include.split = incl.split_par))
}

for (i in 1:5) {
  assign(paste0("mcmc_lb.def_jeff_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                               n.iter = n.iter_par,
                                                               X = X_obs.norm,
                                                               U1 = get(paste0("copula_uu_",i))[1,],
                                                               U2 = get(paste0("copula_uu_",i))[2,],
                                                               mu = 0, 
                                                               sigma = .1, alpha_val = .01, beta_val = .01, 
                                                               prior_list = lb.prior.def, 
                                                               moves.prob = moves.prob_par, 
                                                               starting.tree = NULL,
                                                               cont.unif = cont.unif_par,
                                                               include.split = incl.split_par))
}


for (i in 1:5) {
  assign(paste0("mcmc_lb.def_half_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                               n.iter = n.iter_par,
                                                               X = X_obs.norm,
                                                               U1 = get(paste0("copula_uu_",i))[1,],
                                                               U2 = get(paste0("copula_uu_",i))[2,],
                                                               mu = 0, 
                                                               sigma = .1, alpha_val = .5, beta_val = .5, 
                                                               prior_list = lb.prior.def, 
                                                               moves.prob = moves.prob_par, 
                                                               starting.tree = NULL,
                                                               cont.unif = cont.unif_par,
                                                               include.split = incl.split_par))
}

for (i in 1:5) {
  assign(paste0("mcmc_lb.def_two_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                               n.iter = n.iter_par,
                                                               X = X_obs.norm,
                                                               U1 = get(paste0("copula_uu_",i))[1,],
                                                               U2 = get(paste0("copula_uu_",i))[2,],
                                                               mu = 0, 
                                                               sigma = .1, alpha_val = 2, beta_val = 2, 
                                                               prior_list = lb.prior.def, 
                                                               moves.prob = moves.prob_par, 
                                                               starting.tree = NULL,
                                                               cont.unif = cont.unif_par,
                                                               include.split = incl.split_par))
}

####################
## DEFAULT MODELS ##
####################

model.list.def <- list(
  mcmc_lb.def_two_5)

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

# ggarrange(hist.nl, hist.depth, trace.nl, trace.depth, ncol = 2, nrow = 2)

pred_cond = sapply(1:length(model.list.def$`LB - default`$trees), function(i)get_value_tree(model.list.def[[1]]$trees[[i]],X_obs.norm))

est_par = apply(pred_cond[,-c(as.vector(sapply(0:(n.chain_par-1), function(i) 1:250 + i*n.iter_par)))], 1, mean)

est_par_95 = apply(pred_cond[,-c(as.vector(sapply(0:4, function(i) 1:250 + i*n.iter_par)))], 1, function(x)quantile(x, probs = 0.95))
est_par_05 = apply(pred_cond[,-c(as.vector(sapply(0:4, function(i) 1:250 + i*n.iter_par)))], 1, function(x)quantile(x, probs = 0.05))


copula_uu_pred <- sapply(1:n, function(i)BiCopSim(N=1 , family = 1, par = est_par[i]))

df = data.frame("X" = X_obs, "rho" = rho_true_5, "est_rho" = est_par,
                "est_rho_min" = est_par_05, "est_rho_max" = est_par_95)

ggplot(data = df, aes(x = X)) +
  geom_line(aes(y = rho), color = 2, size = 1) +
  geom_line(aes(y = est_rho), color = 3, size = 1) +
  geom_ribbon(aes(y = est_rho, ymin = est_rho_min, ymax = est_rho_max), alpha = .2) +
  xlab("Observations") +
  theme_bw() +
  theme(legend.key = element_blank()) +
  theme(legend.position = c(1.1,.6), legend.direction = "vertical") +
  theme(legend.title = element_blank())

plot(copula_uu_5[1,], copula_uu_5[2,], xlab = "U1", ylab = "U2")
plot(copula_uu_pred[1,], copula_uu_pred[2,], xlab = "Predicted U1", ylab = "Predicted U2")
