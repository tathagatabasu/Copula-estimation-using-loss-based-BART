# codes and packages
#install.packages("mc2d")
#install.packages(c("copula", "MASS", "coda"))
source('code/import_functions.R')
source('mclapply.R')
source('MCMC_BART_copula.R')
library(data.tree)
library(dplyr)
library(ggplot2)
library(ggpubr)
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
rho_true_1 <- matrix(rho_true_1, ncol = 1)
rm(tree_top)
# monotone
rho_true_2 <- 0.5 + 0.2 * sin(3*X_obs) + 0.3*X_obs^2
# convex
rho_true_3 <- 0.6 + 0.3 * sin(3*X_obs)
# concave
rho_true_4 <- .9 - 0.3 * sin(3*X_obs)
# non-convex
rho_true_5 <- 0.7 - 0.3 * sin(2*X_obs) + 0.2 * sin(4*X_obs) + 0.3 * X_obs^2


plot(X_obs, rho_true_1, xlab = "Observations", ylab = "rho")
plot(X_obs, rho_true_2, xlab = "Observations", ylab = "rho")
plot(X_obs, rho_true_3, xlab = "Observations", ylab = "rho")
plot(X_obs, rho_true_4, xlab = "Observations", ylab = "rho")
plot(X_obs, rho_true_5, xlab = "Observations", ylab = "rho")

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
                                                          log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
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
                                                               sigma = .1, alpha_val = 0, beta_val = 0, 
                                                               log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
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
                                                               log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
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
                                                              log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
                                                               prior_list = lb.prior.def, 
                                                               moves.prob = moves.prob_par, 
                                                               starting.tree = NULL,
                                                               cont.unif = cont.unif_par,
                                                               include.split = incl.split_par))
}

for (i in 1:5) {
  assign(paste0("mcmc_lb.def_LN0.8_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                              n.iter = n.iter_par,
                                                              X = X_obs.norm,
                                                              U1 = get(paste0("copula_uu_",i))[1,],
                                                              U2 = get(paste0("copula_uu_",i))[2,],
                                                              mu = 0, 
                                                              sigma = .1, alpha_val = 2, beta_val = 2, 
                                                              log_nor_mu = 0, log_nor_sigma = 0.8, prior_type = "LN",
                                                              prior_list = lb.prior.def, 
                                                              moves.prob = moves.prob_par, 
                                                              starting.tree = NULL,
                                                              cont.unif = cont.unif_par,
                                                              include.split = incl.split_par))
}

for (i in 1:5) {
  assign(paste0("mcmc_lb.def_LN1_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                n.iter = n.iter_par,
                                                                X = X_obs.norm,
                                                                U1 = get(paste0("copula_uu_",i))[1,],
                                                                U2 = get(paste0("copula_uu_",i))[2,],
                                                                mu = 0, 
                                                                sigma = .1, alpha_val = 2, beta_val = 2, 
                                                                log_nor_mu = 0, log_nor_sigma = 1, prior_type = "LN",
                                                                prior_list = lb.prior.def, 
                                                                moves.prob = moves.prob_par, 
                                                                starting.tree = NULL,
                                                                cont.unif = cont.unif_par,
                                                                include.split = incl.split_par))
}

for (i in 1:5) {
  assign(paste0("mcmc_lb.def_IG11_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                n.iter = n.iter_par,
                                                                X = X_obs.norm,
                                                                U1 = get(paste0("copula_uu_",i))[1,],
                                                                U2 = get(paste0("copula_uu_",i))[2,],
                                                                mu = 0, 
                                                                sigma = .1, alpha_val = 1, beta_val = 1, 
                                                                log_nor_mu = 0, log_nor_sigma = 0.8, prior_type = "IG",
                                                                prior_list = lb.prior.def, 
                                                                moves.prob = moves.prob_par, 
                                                                starting.tree = NULL,
                                                                cont.unif = cont.unif_par,
                                                                include.split = incl.split_par))
}

for (i in 1:5) {
  assign(paste0("mcmc_lb.def_IG22_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                              n.iter = n.iter_par,
                                                              X = X_obs.norm,
                                                              U1 = get(paste0("copula_uu_",i))[1,],
                                                              U2 = get(paste0("copula_uu_",i))[2,],
                                                              mu = 0, 
                                                              sigma = .1, alpha_val = 2, beta_val = 2, 
                                                              log_nor_mu = 0, log_nor_sigma = 1, prior_type = "IG",
                                                              prior_list = lb.prior.def, 
                                                              moves.prob = moves.prob_par, 
                                                              starting.tree = NULL,
                                                              cont.unif = cont.unif_par,
                                                              include.split = incl.split_par))
}

####################
## DEFAULT MODELS ##
####################

test_case = 1

model.list.def <- list(
  get(paste0("mcmc_lb.def_unif_",test_case)),
  get(paste0("mcmc_lb.def_half_",test_case)),
  get(paste0("mcmc_lb.def_jeff_",test_case)),
  get(paste0("mcmc_lb.def_two_",test_case)),
  get(paste0("mcmc_lb.def_LN0.8_",test_case)),
  get(paste0("mcmc_lb.def_LN1_",test_case)),
  get(paste0("mcmc_lb.def_IG22_",test_case)),
  get(paste0("mcmc_lb.def_IG11_",test_case)))

names(model.list.def) <- c(
  'LB - default - unif',
  'LB - default - half',
  'LB - default - jeff',
  'LB - default - two',
  'LB - default - LN0.8',
  'LB - default - LN1',
  'LB - default - IG11',
  'LB - default - IG22')


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
  facet_wrap(facets = ~panel.name, ncol = 2) +
  xlab(~n[L]) +
  ylab('PMF') +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,24,by = 3))

hist.depth <- ggplot(depth.df) +
  geom_histogram(aes(x = y, y = after_stat(density)),
                 binwidth = 1, color = 'black', fill = 'white') +
  facet_wrap(facets = ~panel.name, ncol = 2) +
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

pred_cond = lapply(1:nrow(X_obs.norm), function(i)apply_fun_models(fun_ = function(x)get_value_tree(x, X_obs.norm[i,,drop = FALSE]),
                                                                   mcmc.list = model.list.def,
                                                                   born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par))


pred_cond = do.call(rbind,pred_cond)

pred_cond$obs = as.vector(apply(X_obs.norm, 1, function(x)rep(x, (length(model.list.def) * 1250))))
pred_cond$rho_true = as.vector(apply(get(paste0("rho_true_",test_case)), 1, function(x)rep(x, (length(model.list.def) * 1250))))

ggplot(pred_cond) +
  geom_line(aes(obs, y)) +
  geom_line(aes(obs, rho_true), col = 2) +
  facet_wrap(facets = ~panel.name, ncol = 2) +
  xlab('X') +
  ylab('estimated rho') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30))


pred_cond_mod = pred_cond %>%
  group_by(panel.name, obs, rho_true) %>%
  summarise(rho_mean = mean(y), rho_q95 = quantile(y, .95), rho_q05 = quantile(y, .05)) 

ggplot(pred_cond_mod) +
  geom_line(aes(obs, rho_mean)) +
  geom_line(aes(obs, rho_true), col = 2) +
  geom_line(aes(obs, rho_q95), col = 3) +
  geom_line(aes(obs, rho_q05), col = 3) +
  facet_wrap(facets = ~panel.name, ncol = 2) +
  xlab('X') +
  ylab('estimated rho') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30))

copula_uu_pred <- sapply(1:nrow(pred_cond_mod), function(i)BiCopSim(N=1 , family = 1, par = pred_cond_mod$rho_mean[i]))

pred_cond_mod$U1 = copula_uu_pred[1,]
pred_cond_mod$U2 = copula_uu_pred[2,]

ggplot(pred_cond_mod) +
  geom_point(aes(U1, U2), size = 0.7) +
  facet_wrap(facets = ~panel.name, ncol = 2) +
  xlab('U1') +
  ylab('U2') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30))


data_plot_prior = data.frame(p = seq(-0.999, 0.999, length.out = 1000))

data_plot_prior = data_plot_prior %>% 
  mutate(unif = sapply(p, function(x)exp(logprior_unif(x, 1, 1)))) %>%
  mutate(jeff = sapply(p, function(x)exp(logprior_unif(x, 0, 0)))) %>%
  mutate(two = sapply(p, function(x)exp(logprior_unif(x, 2, 2)))) %>%
  mutate(half = sapply(p, function(x)exp(logprior_unif(x, 0.5, 0.5)))) %>%
  mutate(IG11 = sapply(p, function(x)exp(logprior_inv_gamma(x, 1, 1)))) %>%
  mutate(IG22 = sapply(p, function(x)exp(logprior_inv_gamma(x, 2, 2)))) %>%
  mutate(LN0.8 = sapply(p, function(x)exp(logprior_log_normal(x, 0, 0.8)))) %>%
  mutate(LN1 = sapply(p, function(x)exp(logprior_log_normal(x, 0, 1))))

p_prior = ggplot(data_plot_prior) +
  geom_line(aes(p, unif)) +
  geom_line(aes(p, jeff), col = "red") +
  geom_line(aes(p, two), col = "grey") +
  geom_line(aes(p, half), col = "orange") +
  geom_line(aes(p, IG11), col = "green") +
  geom_line(aes(p, IG22), col = "darkgreen") +
  geom_line(aes(p, LN1), col = "blue") +
  geom_line(aes(p, LN0.8), col = "darkblue") +
  xlab('X') +
  ylab('estimated rho') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30))

p_prior_zoom = p_prior +
  # xlim(c(0.95, 0.99)) +
  ylim(c(0, 20))

p_prior_zoom
