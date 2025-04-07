# codes and packages

#install.packages(c("copula", "MASS", "coda"))
source('code/import_functions.R')
source('mclapply.R')
source('MCMC_BART_copula.R')
library(data.tree)
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

# Define true kendall's tau
tau_true <- 0.6 * sin(3*X_obs) + 0.3

simCopula <- sapply(1:n, function(i)BiCopSim(N=1 , family = 1, par = BiCopTau2Par(1 , tau_true[i])))

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

lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944))
mcmc_lb.def <- multichain_MCMC_copula(n.chain = n.chain_par,
                                         n.iter = n.iter_par,
                                      X = X_obs.norm,
                                      U1 = simCopula[1,],
                                      U2 = simCopula[2,],
                                      Y.var = 0.1, 
                                      mu = 0, 
                                      sigma = 1, 
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

ggarrange(hist.nl, hist.depth, trace.nl, trace.depth, ncol = 2, nrow = 2)

x_new <- matrix(runif(n), ncol = 1)

rownames(x_new) <- 1:nrow(x_new)

tau_new <- 0.6 * sin(3*x_new) + 0.3

pred_cond = sapply(1:length(model.list.def$`LB - default`$trees), function(i)get_value_tree(model.list.def[[1]]$trees[[i]],x_new))

est_par = rowMeans(pred_cond)

plot(est_par, BiCopTau2Par(1 , tau_new), ylim = c(0,1), xlim = c(0,1), xlab = "Predicted rho", ylab = "True rho")

simCopula2 <- sapply(1:n, function(i)BiCopSim(N=1 , family = 1, par = est_par[i]))

plot(simCopula[1,], simCopula[2,], xlab = "U1", ylab = "U2")
plot(simCopula2[1,], simCopula2[2,], xlab = "Predicted U1", ylab = "Predicted U2")

###############################
# alternate
###############################

lb.prior.alt <- list(fun = joint.prior.new.tree, param = c(0.01, 0.01))
mcmc_lb.alt <- multichain_MCMC_copula(n.chain = n.chain_par,
                                      n.iter = n.iter_par,
                                      X = X_obs.norm,
                                      U1 = simCopula[1,],
                                      U2 = simCopula[2,],
                                      Y.var = 0.1, 
                                      mu = 0, 
                                      sigma = 1, 
                                      prior_list = lb.prior.alt, 
                                      moves.prob = moves.prob_par, 
                                      starting.tree = NULL,
                                      cont.unif = cont.unif_par,
                                      include.split = incl.split_par)

####################
## ALT MODELS ##
####################

model.list.alt <- list(
  mcmc_lb.alt)

names(model.list.alt) <- c(
  'LB - alt')


# extract depth, number of terminal nodes, missing rate and loglik of all the trees
depth.df.alt <- apply_fun_models(fun_ = get_depth, 
                             mcmc.list = model.list.alt,
                             born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)
nterm.df.alt <- apply_fun_models(fun_ = get_num_terminal_nodes, 
                             mcmc.list = model.list.alt, 
                             born.out.pc = 250, n.chain = n.chain_par, sample.pc = n.iter_par)

hist.nl.alt <- ggplot(nterm.df.alt) + 
  geom_histogram(aes(x = y, y = after_stat(density)), 
                 binwidth = 1, color = 'black', fill = 'white') + 
  facet_wrap(facets = ~panel.name) + 
  xlab(~n[L]) + 
  ylab('PMF') + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,24,by = 3))  

hist.depth.alt <- ggplot(depth.df.alt) + 
  geom_histogram(aes(x = y, y = after_stat(density)), 
                 binwidth = 1, color = 'black', fill = 'white') + 
  facet_wrap(facets = ~panel.name) + 
  xlab('Depth') + 
  ylab('PMF') + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,10,by = 1))

trace.nl.alt <- ggplot(nterm.df.alt) + 
  geom_vline(xintercept = seq(0,1250,by = 50), 
             color = 'grey', size = 0.2, alpha=0.75)+
  geom_line(aes(x, y)) + 
  facet_wrap(facets = ~panel.name, ncol = 2) + 
  xlab('Iteration') + 
  ylab(~n[L]) + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,1250,by = 50)) + 
  theme(axis.text.x = element_text(angle = 30))

trace.depth.alt <- ggplot(depth.df.alt) + 
  geom_vline(xintercept = seq(0,1250,by = 50),
             color = 'grey', size = 0.2, alpha=0.75)+
  geom_line(aes(x, y)) +
  facet_wrap(facets = ~panel.name, ncol = 2) + 
  xlab('Iteration') + 
  ylab('Depth') + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,1250,by = 50)) + 
  theme(axis.text.x = element_text(angle = 30))

hist.nl.alt
hist.depth.alt

trace.nl.alt
trace.depth.alt

ggarrange(hist.nl.alt, hist.depth.alt, trace.nl.alt, trace.depth.alt, ncol = 2, nrow = 2)


pred_cond_alt = sapply(1:length(model.list.alt$`LB - alt`$trees), function(i)get_value_tree(model.list.alt[[1]]$trees[[i]],x_new))

est_par_alt = rowMeans(pred_cond_alt)

plot(est_par_alt, BiCopTau2Par(1 , tau_new), ylim = c(0,1), xlim = c(0,1), xlab = "Predicted rho", ylab = "True rho")

simCopula2_alt <- sapply(1:n, function(i)BiCopSim(N=1 , family = 1, par = est_par_alt[i]))

plot(simCopula[1,], simCopula[2,], xlab = "U1", ylab = "U2")
plot(simCopula2_alt[1,], simCopula2_alt[2,], xlab = "Predicted U1", ylab = "Predicted U2")
