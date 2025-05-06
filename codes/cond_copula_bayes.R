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
library(gplots)
library(xtable)


# prior

if(F){
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
    xlab('rho') +
    ylab('p.d.f.') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30))
  
  p_prior_zoom = p_prior +
    # xlim(c(0.95, 0.99)) +
    ylim(c(0, 20))
  
  p_prior_zoom
}

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
tree_top <- assign_term_node_values(tree_top, 0.6, 0.1)
tau_true_1 <- sample_CART(tree_top, X_obs, sigma_ = 0.001) 
tau_true_1 <- matrix(tau_true_1, ncol = 1)
rm(tree_top)
# monotone
tau_true_2 <- 0.3 + 0.2 * sin(3*X_obs) + 0.3*X_obs^2
# convex
tau_true_3 <- 0.5 + 0.3 * sin(3*X_obs)
# concave
tau_true_4 <- 0.8 - 0.3 * sin(3*X_obs)
# non-convex
tau_true_5 <- 0.6 - 0.3 * sin(2*X_obs) + 0.2 * sin(4*X_obs) + 0.3 * X_obs^2


plot(X_obs, tau_true_1, xlab = "Observations", ylab = "tau")
plot(X_obs, tau_true_2, xlab = "Observations", ylab = "tau")
plot(X_obs, tau_true_3, xlab = "Observations", ylab = "tau")
plot(X_obs, tau_true_4, xlab = "Observations", ylab = "tau")
plot(X_obs, tau_true_5, xlab = "Observations", ylab = "tau")

#################################################################
# simulation studies

# gaussian
if(T){
  for (i in 1:5) {
    assign(paste0("copula_uu_gauss_",i), BiCopSim(n, family = 1, par = BiCopTau2Par(family = 1, get(paste0("tau_true_",i)))))
  }
  
  
  # Create 2D histogram (bins only)
  hist_1 <- hist2d(copula_uu_gauss_1[,1], copula_uu_gauss_1[,2], nbins = c(10,10), show = FALSE)
  hist_2 <- hist2d(copula_uu_gauss_2[,1], copula_uu_gauss_2[,2], nbins = c(10,10), show = FALSE)
  hist_3 <- hist2d(copula_uu_gauss_3[,1], copula_uu_gauss_3[,2], nbins = c(10,10), show = FALSE)
  hist_4 <- hist2d(copula_uu_gauss_4[,1], copula_uu_gauss_4[,2], nbins = c(10,10), show = FALSE)
  hist_5 <- hist2d(copula_uu_gauss_5[,1], copula_uu_gauss_5[,2], nbins = c(10,10), show = FALSE)
  
  # Plot 3D histogram
  
  par(mar = c(2,2,1,1))
  
  plot(copula_uu_gauss_1[,1], copula_uu_gauss_1[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_1$x,
    y = hist_1$y,
    z = hist_1$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
  plot(copula_uu_gauss_2[,1], copula_uu_gauss_2[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_2$x,
    y = hist_2$y,
    z = hist_2$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
  plot(copula_uu_gauss_3[,1], copula_uu_gauss_3[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_3$x,
    y = hist_3$y,
    z = hist_3$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
  plot(copula_uu_gauss_4[,1], copula_uu_gauss_4[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_4$x,
    y = hist_4$y,
    z = hist_4$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
  plot(copula_uu_gauss_5[,1], copula_uu_gauss_5[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_5$x,
    y = hist_5$y,
    z = hist_5$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
  ##################################################
  # normalise predictors 
  X_obs.norm <- as.data.frame(apply(X_obs, 2, \(x) (x - min(x))/(max(x) - min(x))))
  X_obs.norm <- as.matrix(X_obs.norm)
  rownames(X_obs.norm) <- 1:nrow(X_obs)
  
  n.chain_par <- 5
  n.iter_par <- 50000
  incl.split_par <- TRUE
  cont.unif_par <- TRUE
  moves.prob_par <- c(0.1, 0.4, 0.25, 0.25)
  
  lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) # c(1.5618883, 0.6293944)
  
  for (i in 1:5) {
    assign(paste0("gauss_mcmc_lb.def_unif_",i), MCMC_copula(#n.chain = n.chain_par,
                                                                 n.iter = n.iter_par,
                                                                 X = X_obs.norm,
                                                                 U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                                 U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                                 mu = 0, cop_type = "Gauss", 
                                                                 sigma = .5, alpha_val = 0, beta_val = 0, 
                                                                 log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
                                                                 prior_list = lb.prior.def, 
                                                                 moves.prob = moves.prob_par, 
                                                                 starting.tree = NULL,
                                                                 cont.unif = cont.unif_par,
                                                                 include.split = incl.split_par))
  }
  
  for (i in 1:5) {
    assign(paste0("gauss_mcmc_lb.def_jeff_",i), MCMC_copula(#n.chain = n.chain_par,
                                                                 n.iter = n.iter_par,
                                                                 X = X_obs.norm,
                                                                 U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                                 U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                                 mu = 0, cop_type = "Gauss", 
                                                                 sigma = .1, alpha_val = 0, beta_val = 0, 
                                                                 log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
                                                                 prior_list = lb.prior.def, 
                                                                 moves.prob = moves.prob_par, 
                                                                 starting.tree = NULL,
                                                                 cont.unif = cont.unif_par,
                                                                 include.split = incl.split_par))
  }
  
  
  for (i in 1:5) {
    assign(paste0("gauss_mcmc_lb.def_half_",i), MCMC_copula(#n.chain = n.chain_par,
                                                                 n.iter = n.iter_par,
                                                                 X = X_obs.norm,
                                                                 U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                                 U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                                 mu = 0, cop_type = "Gauss", 
                                                                 sigma = .1, alpha_val = .5, beta_val = .5, 
                                                                 log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
                                                                 prior_list = lb.prior.def, 
                                                                 moves.prob = moves.prob_par, 
                                                                 starting.tree = NULL,
                                                                 cont.unif = cont.unif_par,
                                                                 include.split = incl.split_par))
  }
  
  for (i in 1:5) {
    assign(paste0("gauss_mcmc_lb.def_two_",i), MCMC_copula(#n.chain = n.chain_par,
                                                                n.iter = n.iter_par,
                                                                X = X_obs.norm,
                                                                U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                                U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                                mu = 0, cop_type = "Gauss", 
                                                                sigma = .1, alpha_val = 2, beta_val = 2, 
                                                                log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
                                                                prior_list = lb.prior.def, 
                                                                moves.prob = moves.prob_par, 
                                                                starting.tree = NULL,
                                                                cont.unif = cont.unif_par,
                                                                include.split = incl.split_par))
  }
  
  for (i in 1:5) {
    assign(paste0("gauss_mcmc_lb.def_LN0.8_",i), MCMC_copula(#n.chain = n.chain_par,
                                                                  n.iter = n.iter_par,
                                                                  X = X_obs.norm,
                                                                  U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                                  U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                                  mu = 0, cop_type = "Gauss", 
                                                                  sigma = .1, alpha_val = 2, beta_val = 2, 
                                                                  log_nor_mu = 0, log_nor_sigma = 0.8, prior_type = "LN",
                                                                  prior_list = lb.prior.def, 
                                                                  moves.prob = moves.prob_par, 
                                                                  starting.tree = NULL,
                                                                  cont.unif = cont.unif_par,
                                                                  include.split = incl.split_par))
  }
  
  for (i in 1:5) {
    assign(paste0("gauss_mcmc_lb.def_LN1_",i), MCMC_copula(#n.chain = n.chain_par,
                                                                n.iter = n.iter_par,
                                                                X = X_obs.norm,
                                                                U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                                U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                                mu = 0, cop_type = "Gauss", 
                                                                sigma = .1, alpha_val = 2, beta_val = 2, 
                                                                log_nor_mu = 0, log_nor_sigma = 1, prior_type = "LN",
                                                                prior_list = lb.prior.def, 
                                                                moves.prob = moves.prob_par, 
                                                                starting.tree = NULL,
                                                                cont.unif = cont.unif_par,
                                                                include.split = incl.split_par))
  }
  
  for (i in 1:5) {
    assign(paste0("gauss_mcmc_lb.def_IG11_",i), MCMC_copula(#n.chain = n.chain_par,
                                                                 n.iter = n.iter_par,
                                                                 X = X_obs.norm,
                                                                 U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                                 U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                                 mu = 0, cop_type = "Gauss", 
                                                                 sigma = .1, alpha_val = 1, beta_val = 1, 
                                                                 log_nor_mu = 0, log_nor_sigma = 0.8, prior_type = "IG",
                                                                 prior_list = lb.prior.def, 
                                                                 moves.prob = moves.prob_par, 
                                                                 starting.tree = NULL,
                                                                 cont.unif = cont.unif_par,
                                                                 include.split = incl.split_par))
  }
  
  for (i in 1:5) {
    assign(paste0("gauss_mcmc_lb.def_IG22_",i), MCMC_copula(#n.chain = n.chain_par,
                                                                 n.iter = n.iter_par,
                                                                 X = X_obs.norm,
                                                                 U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                                 U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                                 mu = 0, cop_type = "Gauss", 
                                                                 sigma = .1, alpha_val = 2, beta_val = 2, 
                                                                 log_nor_mu = 0, log_nor_sigma = 1, prior_type = "IG",
                                                                 prior_list = lb.prior.def, 
                                                                 moves.prob = moves.prob_par, 
                                                                 starting.tree = NULL,
                                                                 cont.unif = cont.unif_par,
                                                                 include.split = incl.split_par))
  }
  
}

# frank
if(T){
  for (i in 1:5) {
    assign(paste0("copula_uu_frank_",i), BiCopSim(n, family = 5, par = BiCopTau2Par(family = 5, get(paste0("tau_true_",i)))))
  }
  
  
  # Create 2D histogram (bins only)
  hist_1 <- hist2d(copula_uu_frank_1[,1], copula_uu_frank_1[,2], nbins = c(10,10), show = FALSE)
  hist_2 <- hist2d(copula_uu_frank_2[,1], copula_uu_frank_2[,2], nbins = c(10,10), show = FALSE)
  hist_3 <- hist2d(copula_uu_frank_3[,1], copula_uu_frank_3[,2], nbins = c(10,10), show = FALSE)
  hist_4 <- hist2d(copula_uu_frank_4[,1], copula_uu_frank_4[,2], nbins = c(10,10), show = FALSE)
  hist_5 <- hist2d(copula_uu_frank_5[,1], copula_uu_frank_5[,2], nbins = c(10,10), show = FALSE)
  
  # Plot 3D histogram
  
  par(mar = c(2,2,1,1))
  
  plot(copula_uu_frank_1[,1], copula_uu_frank_1[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_1$x,
    y = hist_1$y,
    z = hist_1$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
  plot(copula_uu_frank_2[,1], copula_uu_frank_2[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_2$x,
    y = hist_2$y,
    z = hist_2$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
  plot(copula_uu_frank_3[,1], copula_uu_frank_3[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_3$x,
    y = hist_3$y,
    z = hist_3$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
  plot(copula_uu_frank_4[,1], copula_uu_frank_4[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_4$x,
    y = hist_4$y,
    z = hist_4$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
  plot(copula_uu_frank_5[,1], copula_uu_frank_5[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_5$x,
    y = hist_5$y,
    z = hist_5$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
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
  
  lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) 
  
  for (i in 1:5) {
    assign(paste0("frank_mcmc_lb.def_unif_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                 n.iter = n.iter_par,
                                                                 X = X_obs.norm,
                                                                 U1 = get(paste0("copula_uu_frank_",i))[,1],
                                                                 U2 = get(paste0("copula_uu_frank_",i))[,2],
                                                                 mu = 0, cop_type = "Frank", 
                                                                 sigma = .5, alpha_val = 1, beta_val = 1, 
                                                                 log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
                                                                 prior_list = lb.prior.def, 
                                                                 moves.prob = moves.prob_par, 
                                                                 starting.tree = NULL,
                                                                 cont.unif = cont.unif_par,
                                                                 include.split = incl.split_par))
  }
  
  for (i in 1:5) {
    assign(paste0("frank_mcmc_lb.def_jeff_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                 n.iter = n.iter_par,
                                                                 X = X_obs.norm,
                                                                 U1 = get(paste0("copula_uu_frank_",i))[,1],
                                                                 U2 = get(paste0("copula_uu_frank_",i))[,2],
                                                                 mu = 0, cop_type = "Frank", 
                                                                 sigma = .1, alpha_val = 0, beta_val = 0, 
                                                                 log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
                                                                 prior_list = lb.prior.def, 
                                                                 moves.prob = moves.prob_par, 
                                                                 starting.tree = NULL,
                                                                 cont.unif = cont.unif_par,
                                                                 include.split = incl.split_par))
  }
  
  
  for (i in 1:5) {
    assign(paste0("frank_mcmc_lb.def_normal_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                   n.iter = n.iter_par,
                                                                   X = X_obs.norm,
                                                                   U1 = get(paste0("copula_uu_frank_",i))[,1],
                                                                   U2 = get(paste0("copula_uu_frank_",i))[,2],
                                                                   mu = 0, cop_type = "Frank", 
                                                                   sigma = .5, alpha_val = .5, beta_val = .5, 
                                                                   log_nor_mu = 0, log_nor_sigma = 20, prior_type = "Normal",
                                                                   prior_list = lb.prior.def, 
                                                                   moves.prob = moves.prob_par, 
                                                                   starting.tree = NULL,
                                                                   cont.unif = cont.unif_par,
                                                                   include.split = incl.split_par))
  }
  
  for (i in 1:5) {
    assign(paste0("frank_mcmc_lb.def_t_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                              n.iter = n.iter_par,
                                                              X = X_obs.norm,
                                                              U1 = get(paste0("copula_uu_frank_",i))[,1],
                                                              U2 = get(paste0("copula_uu_frank_",i))[,2],
                                                              mu = 0, cop_type = "Frank", 
                                                              sigma = .5, alpha_val = 2, beta_val = 2, 
                                                              log_nor_mu = 0, log_nor_sigma = 20, prior_type = "T",
                                                              prior_list = lb.prior.def, 
                                                              moves.prob = moves.prob_par, 
                                                              starting.tree = NULL,
                                                              cont.unif = cont.unif_par,
                                                              include.split = incl.split_par))
  }
  
  pp.low.omeg.high.gam <- c(0.5, 1.5)
  lossb.prior.list.high.gam <- list(fun = joint.prior.new.tree, param = pp.low.omeg.high.gam)
  
  for (i in 1:5) {
    assign(paste0("frank_mcmc_lossb.prior.high.gam_unif_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                               n.iter = n.iter_par,
                                                                               X = X_obs.norm,
                                                                               U1 = get(paste0("copula_uu_frank_",i))[,1],
                                                                               U2 = get(paste0("copula_uu_frank_",i))[,2],
                                                                               mu = 0, cop_type = "Frank", 
                                                                               sigma = .5, alpha_val = 1, beta_val = 1, 
                                                                               log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
                                                                               prior_list = lossb.prior.list.high.gam, 
                                                                               moves.prob = moves.prob_par, 
                                                                               starting.tree = NULL,
                                                                               cont.unif = cont.unif_par,
                                                                               include.split = incl.split_par))
  }
  
  for (i in 1:5) {
    assign(paste0("frank_mcmc_lossb.prior.high.gam_jeff_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                               n.iter = n.iter_par,
                                                                               X = X_obs.norm,
                                                                               U1 = get(paste0("copula_uu_frank_",i))[,1],
                                                                               U2 = get(paste0("copula_uu_frank_",i))[,2],
                                                                               mu = 0, cop_type = "Frank", 
                                                                               sigma = .1, alpha_val = 0, beta_val = 0, 
                                                                               log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
                                                                               prior_list = lossb.prior.list.high.gam, 
                                                                               moves.prob = moves.prob_par, 
                                                                               starting.tree = NULL,
                                                                               cont.unif = cont.unif_par,
                                                                               include.split = incl.split_par))
  }
  
  
  for (i in 1:5) {
    assign(paste0("frank_mcmc_lossb.prior.high.gam_normal_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                                 n.iter = n.iter_par,
                                                                                 X = X_obs.norm,
                                                                                 U1 = get(paste0("copula_uu_frank_",i))[,1],
                                                                                 U2 = get(paste0("copula_uu_frank_",i))[,2],
                                                                                 mu = 0, cop_type = "Frank", 
                                                                                 sigma = .5, alpha_val = .5, beta_val = .5, 
                                                                                 log_nor_mu = 0, log_nor_sigma = 20, prior_type = "Normal",
                                                                                 prior_list = lossb.prior.list.high.gam, 
                                                                                 moves.prob = moves.prob_par, 
                                                                                 starting.tree = NULL,
                                                                                 cont.unif = cont.unif_par,
                                                                                 include.split = incl.split_par))
  }
  
  for (i in 1:5) {
    assign(paste0("frank_mcmc_lossb.prior.high.gam_t_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                            n.iter = n.iter_par,
                                                                            X = X_obs.norm,
                                                                            U1 = get(paste0("copula_uu_frank_",i))[,1],
                                                                            U2 = get(paste0("copula_uu_frank_",i))[,2],
                                                                            mu = 0, cop_type = "Frank", 
                                                                            sigma = .5, alpha_val = 2, beta_val = 2, 
                                                                            log_nor_mu = 0, log_nor_sigma = 20, prior_type = "T",
                                                                            prior_list = lossb.prior.list.high.gam, 
                                                                            moves.prob = moves.prob_par, 
                                                                            starting.tree = NULL,
                                                                            cont.unif = cont.unif_par,
                                                                            include.split = incl.split_par))
  }
  
}

# clayton
if(T){
  for (i in 1:5) {
    assign(paste0("copula_uu_clayton_",i), BiCopSim(n, family = 3, par = BiCopTau2Par(family = 3, get(paste0("tau_true_",i)))))
  }
  
  
  # Create 2D histogram (bins only)
  hist_1 <- hist2d(copula_uu_clayton_1[,1], copula_uu_clayton_1[,2], nbins = c(10,10), show = FALSE)
  hist_2 <- hist2d(copula_uu_clayton_2[,1], copula_uu_clayton_2[,2], nbins = c(10,10), show = FALSE)
  hist_3 <- hist2d(copula_uu_clayton_3[,1], copula_uu_clayton_3[,2], nbins = c(10,10), show = FALSE)
  hist_4 <- hist2d(copula_uu_clayton_4[,1], copula_uu_clayton_4[,2], nbins = c(10,10), show = FALSE)
  hist_5 <- hist2d(copula_uu_clayton_5[,1], copula_uu_clayton_5[,2], nbins = c(10,10), show = FALSE)
  
  # Plot 3D histogram
  
  par(mar = c(2,2,1,1))
  
  plot(copula_uu_clayton_1[,1], copula_uu_clayton_1[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_1$x,
    y = hist_1$y,
    z = hist_1$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
  plot(copula_uu_clayton_2[,1], copula_uu_clayton_2[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_2$x,
    y = hist_2$y,
    z = hist_2$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
  plot(copula_uu_clayton_3[,1], copula_uu_clayton_3[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_3$x,
    y = hist_3$y,
    z = hist_3$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
  plot(copula_uu_clayton_4[,1], copula_uu_clayton_4[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_4$x,
    y = hist_4$y,
    z = hist_4$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
  plot(copula_uu_clayton_5[,1], copula_uu_clayton_5[,2], xlab = "U1", ylab = "U2")
  
  hist3D(
    x = hist_5$x,
    y = hist_5$y,
    z = hist_5$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = ""
  )
  
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
  
  lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) 
  
  for (i in 1:5) {
    assign(paste0("clayton_mcmc_lb.def_IG22_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                         n.iter = n.iter_par,
                                                                         X = X_obs.norm,
                                                                         U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                         U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                         mu = 0, cop_type = "Clayton", 
                                                                         sigma = .5, alpha_val = 2, beta_val = 2, 
                                                                         log_nor_mu = 0, log_nor_sigma = 1, prior_type = "IG",
                                                                         prior_list = lb.prior.def, 
                                                                         moves.prob = moves.prob_par, 
                                                                         starting.tree = NULL,
                                                                         cont.unif = cont.unif_par,
                                                                         include.split = incl.split_par))
  }
  
  for (i in 1:5) {
    assign(paste0("clayton_mcmc_lb.def_IG11_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                         n.iter = n.iter_par,
                                                                         X = X_obs.norm,
                                                                         U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                         U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                         mu = 0, cop_type = "Clayton", 
                                                                         sigma = .1, alpha_val = 1, beta_val = 1, 
                                                                         log_nor_mu = 0, log_nor_sigma = 1, prior_type = "IG",
                                                                         prior_list = lb.prior.def, 
                                                                         moves.prob = moves.prob_par, 
                                                                         starting.tree = NULL,
                                                                         cont.unif = cont.unif_par,
                                                                         include.split = incl.split_par))
  }
  
  
  for (i in 1:5) {
    assign(paste0("clayton_mcmc_lb.def_LN01_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                         n.iter = n.iter_par,
                                                                         X = X_obs.norm,
                                                                         U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                         U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                         mu = 0, cop_type = "Clayton", 
                                                                         sigma = .5, alpha_val = .5, beta_val = .5, 
                                                                         log_nor_mu = 0, log_nor_sigma = 1, prior_type = "LN",
                                                                         prior_list = lb.prior.def, 
                                                                         moves.prob = moves.prob_par, 
                                                                         starting.tree = NULL,
                                                                         cont.unif = cont.unif_par,
                                                                         include.split = incl.split_par))
  }
  
  for (i in 1:5) {
    assign(paste0("clayton_mcmc_lb.def_LN05_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                         n.iter = n.iter_par,
                                                                         X = X_obs.norm,
                                                                         U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                         U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                         mu = 0, cop_type = "Clayton", 
                                                                         sigma = .5, alpha_val = 2, beta_val = 2, 
                                                                         log_nor_mu = 0, log_nor_sigma = 5, prior_type = "LN",
                                                                         prior_list = lb.prior.def, 
                                                                         moves.prob = moves.prob_par, 
                                                                         starting.tree = NULL,
                                                                         cont.unif = cont.unif_par,
                                                                         include.split = incl.split_par))
  }
  
  lossb.prior.list.high.gam <- list(fun = joint.prior.new.tree, param = pp.low.omeg.high.gam)
  
  for (i in 1:5) {
    assign(paste0("clayton_mcmc_lossb.prior.high.gam_IG11_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                                       n.iter = n.iter_par,
                                                                                       X = X_obs.norm,
                                                                                       U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                                       U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                                       mu = 0, cop_type = "Clayton", 
                                                                                       sigma = .5, alpha_val = 1, beta_val = 1, 
                                                                                       log_nor_mu = 0, log_nor_sigma = 1, prior_type = "IG",
                                                                                       prior_list = lossb.prior.list.high.gam, 
                                                                                       moves.prob = moves.prob_par, 
                                                                                       starting.tree = NULL,
                                                                                       cont.unif = cont.unif_par,
                                                                                       include.split = incl.split_par))
  }
  
  for (i in 1:5) {
    assign(paste0("clayton_mcmc_lossb.prior.high.gam_IG22_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                                       n.iter = n.iter_par,
                                                                                       X = X_obs.norm,
                                                                                       U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                                       U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                                       mu = 0, cop_type = "Clayton", 
                                                                                       sigma = .1, alpha_val = 2, beta_val = 2, 
                                                                                       log_nor_mu = 0, log_nor_sigma = 1, prior_type = "IG",
                                                                                       prior_list = lossb.prior.list.high.gam, 
                                                                                       moves.prob = moves.prob_par, 
                                                                                       starting.tree = NULL,
                                                                                       cont.unif = cont.unif_par,
                                                                                       include.split = incl.split_par))
  }
  
  
  for (i in 1:5) {
    assign(paste0("clayton_mcmc_lossb.prior.high.gam_LN01_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                                       n.iter = n.iter_par,
                                                                                       X = X_obs.norm,
                                                                                       U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                                       U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                                       mu = 0, cop_type = "Clayton", 
                                                                                       sigma = .5, alpha_val = .5, beta_val = .5, 
                                                                                       log_nor_mu = 0, log_nor_sigma = 1, prior_type = "LN",
                                                                                       prior_list = lossb.prior.list.high.gam, 
                                                                                       moves.prob = moves.prob_par, 
                                                                                       starting.tree = NULL,
                                                                                       cont.unif = cont.unif_par,
                                                                                       include.split = incl.split_par))
  }
  
  for (i in 1:5) {
    assign(paste0("clayton_mcmc_lossb.prior.high.gam_LN05_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                                       n.iter = n.iter_par,
                                                                                       X = X_obs.norm,
                                                                                       U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                                       U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                                       mu = 0, cop_type = "Clayton", 
                                                                                       sigma = .5, alpha_val = 2, beta_val = 2, 
                                                                                       log_nor_mu = 0, log_nor_sigma = 5, prior_type = "LN",
                                                                                       prior_list = lossb.prior.list.high.gam, 
                                                                                       moves.prob = moves.prob_par, 
                                                                                       starting.tree = NULL,
                                                                                       cont.unif = cont.unif_par,
                                                                                       include.split = incl.split_par))
  }
  
}
###################################################
# results

# gauss

if(F){
  test_case = 1
  
  n.born.out.par = 250
  
  model.list.def <- list(
    get(paste0("gauss_mcmc_lb.def_unif_",test_case)),
    get(paste0("gauss_mcmc_lb.def_half_",test_case)),
    get(paste0("gauss_mcmc_lb.def_jeff_",test_case)),
    get(paste0("gauss_mcmc_lb.def_two_",test_case)),
    get(paste0("gauss_mcmc_lb.def_LN0.8_",test_case)),
    get(paste0("gauss_mcmc_lb.def_LN1_",test_case)),
    get(paste0("gauss_mcmc_lb.def_IG22_",test_case)),
    get(paste0("gauss_mcmc_lb.def_IG11_",test_case)))
  
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
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
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
    geom_vline(xintercept = seq(0,n.chain_par * (n.iter_par - n.born.out.par),by = 50),
               color = 'grey', size = 0.2, alpha=0.75)+
    geom_line(aes(x, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab(~n[L]) +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,n.chain_par * (n.iter_par - n.born.out.par),by = 50)) +
    theme(axis.text.x = element_text(angle = 30))
  
  trace.depth <- ggplot(depth.df) +
    geom_vline(xintercept = seq(0,n.chain_par * (n.iter_par - n.born.out.par),by = 50),
               color = 'grey', size = 0.2, alpha=0.75)+
    geom_line(aes(x, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab('Depth') +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,n.chain_par * (n.iter_par - n.born.out.par),by = 50)) +
    theme(axis.text.x = element_text(angle = 30))
  
  hist.nl
  hist.depth
  
  trace.nl
  trace.depth
  
  pred_cond = lapply(1:nrow(X_obs.norm), function(i)apply_fun_models(fun_ = function(x)get_value_tree(x, X_obs.norm[i,,drop = FALSE]),
                                                                     mcmc.list = model.list.def,
                                                                     born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par))
  
  
  pred_cond = do.call(rbind,pred_cond)
  
  pred_cond$obs = as.vector(apply(X_obs.norm, 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  pred_cond$theta_true = as.vector(apply(BiCopTau2Par(family = 1, get(paste0("tau_true_",test_case))), 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  
  ggplot(pred_cond) +
    geom_line(aes(obs, y)) +
    geom_line(aes(obs, theta_true), col = 2) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated theta') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30))
  
  
  pred_cond_mod = pred_cond %>%
    group_by(panel.name, obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .95), theta_q025 = quantile(y, .05)) 
  
  ggplot(pred_cond_mod) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated theta') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30))
  
  copula_uu_gauss_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 1, par = pred_cond_mod$theta_mean)
  
  pred_cond_mod$U1 = copula_uu_gauss_pred[,1]
  pred_cond_mod$U2 = copula_uu_gauss_pred[,2]
  
  pred_cond_mod_p1 = pred_cond_mod %>% filter(panel.name == "LB - default - unif")
  pred_cond_mod_p2 = pred_cond_mod %>% filter(panel.name == "LB - default - half")
  pred_cond_mod_p3 = pred_cond_mod %>% filter(panel.name == "LB - default - jeff")
  pred_cond_mod_p4 = pred_cond_mod %>% filter(panel.name == "LB - default - two")
  pred_cond_mod_p5 = pred_cond_mod %>% filter(panel.name == "LB - default - LN0.8")
  pred_cond_mod_p6 = pred_cond_mod %>% filter(panel.name == "LB - default - LN1")
  pred_cond_mod_p7 = pred_cond_mod %>% filter(panel.name == "LB - default - IG11")
  pred_cond_mod_p8 = pred_cond_mod %>% filter(panel.name == "LB - default - IG22")
  
  hist_pred1 <- hist2d(pred_cond_mod_p1$U1, pred_cond_mod_p1$U2, nbins = c(10,10), show = FALSE)
  hist_pred2 <- hist2d(pred_cond_mod_p2$U1, pred_cond_mod_p2$U2, nbins = c(10,10), show = FALSE)
  hist_pred3 <- hist2d(pred_cond_mod_p3$U1, pred_cond_mod_p3$U2, nbins = c(10,10), show = FALSE)
  hist_pred4 <- hist2d(pred_cond_mod_p4$U1, pred_cond_mod_p4$U2, nbins = c(10,10), show = FALSE)
  hist_pred5 <- hist2d(pred_cond_mod_p5$U1, pred_cond_mod_p5$U2, nbins = c(10,10), show = FALSE)
  hist_pred6 <- hist2d(pred_cond_mod_p6$U1, pred_cond_mod_p6$U2, nbins = c(10,10), show = FALSE)
  hist_pred7 <- hist2d(pred_cond_mod_p7$U1, pred_cond_mod_p7$U2, nbins = c(10,10), show = FALSE)
  hist_pred8 <- hist2d(pred_cond_mod_p8$U1, pred_cond_mod_p8$U2, nbins = c(10,10), show = FALSE)
  
  hist3D(
    x = hist_pred1$x,
    y = hist_pred1$y,
    z = hist_pred1$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = unique(pred_cond_mod_p1$panel.name)
  )
  
  hist3D(
    x = hist_pred2$x,
    y = hist_pred2$y,
    z = hist_pred2$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = unique(pred_cond_mod_p2$panel.name)
  )
  
  hist3D(
    x = hist_pred3$x,
    y = hist_pred3$y,
    z = hist_pred3$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = unique(pred_cond_mod_p3$panel.name)
  )
  
  hist3D(
    x = hist_pred4$x,
    y = hist_pred4$y,
    z = hist_pred4$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = unique(pred_cond_mod_p4$panel.name)
  )
  
  hist3D(
    x = hist_pred5$x,
    y = hist_pred5$y,
    z = hist_pred5$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = unique(pred_cond_mod_p5$panel.name)
  )
  
  hist3D(
    x = hist_pred6$x,
    y = hist_pred6$y,
    z = hist_pred6$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = unique(pred_cond_mod_p6$panel.name)
  )
  
  hist3D(
    x = hist_pred7$x,
    y = hist_pred7$y,
    z = hist_pred7$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = unique(pred_cond_mod_p7$panel.name)
  )
  
  hist3D(
    x = hist_pred8$x,
    y = hist_pred8$y,
    z = hist_pred8$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = unique(pred_cond_mod_p8$panel.name)
  )
  
  pred_cond_stat = pred_cond_mod %>%
    group_by(panel.name)%>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = pred_cond_stat %>%
    group_by(panel.name)%>%
    summarise_all(mean)
  
  xtable(pred_cond_summary , digits = 4)
}

# frank

if(F){
  test_case = 5
  
  n.born.out.par = 250
  
  model.list.def <- list(
    get(paste0("frank_mcmc_lb.def_unif_",test_case)),
    get(paste0("frank_mcmc_lb.def_jeff_",test_case)),
    get(paste0("frank_mcmc_lb.def_normal_",test_case)),
    get(paste0("frank_mcmc_lb.def_t_",test_case)),
    get(paste0("frank_mcmc_lossb.prior.high.gam_unif_",test_case)),
    get(paste0("frank_mcmc_lossb.prior.high.gam_jeff_",test_case)),
    get(paste0("frank_mcmc_lossb.prior.high.gam_normal_",test_case)),
    get(paste0("frank_mcmc_lossb.prior.high.gam_t_",test_case)))
  
  names(model.list.def) <- c(
    'LB - default - unif',
    'LB - default - jeff',
    'LB - default - normal',
    'LB - default - t',
    'LB - high.gam - unif',
    'LB - high.gam - jeff',
    'LB - high.gam - normal',
    'LB - high.gam - t')
  
  
  # extract depth, number of terminal nodes, missing rate and loglik of all the trees
  depth.df <- apply_fun_models(fun_ = get_depth,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
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
    geom_vline(xintercept = seq(0,n.chain_par * (n.iter_par - n.born.out.par),by = 50),
               color = 'grey', size = 0.2, alpha=0.75)+
    geom_line(aes(x, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab(~n[L]) +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,n.chain_par * (n.iter_par - n.born.out.par),by = 50)) +
    theme(axis.text.x = element_text(angle = 30))
  
  trace.depth <- ggplot(depth.df) +
    geom_vline(xintercept = seq(0,n.chain_par * (n.iter_par - n.born.out.par),by = 50),
               color = 'grey', size = 0.2, alpha=0.75)+
    geom_line(aes(x, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab('Depth') +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,n.chain_par * (n.iter_par - n.born.out.par),by = 50)) +
    theme(axis.text.x = element_text(angle = 30))
  
  hist.nl
  hist.depth
  
  trace.nl
  trace.depth
  
  pred_cond = lapply(1:nrow(X_obs.norm), function(i)apply_fun_models(fun_ = function(x)get_value_tree(x, X_obs.norm[i,,drop = FALSE]),
                                                                     mcmc.list = model.list.def,
                                                                     born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par))
  
  
  pred_cond = do.call(rbind,pred_cond)
  
  pred_cond$obs = as.vector(apply(X_obs.norm, 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  pred_cond$theta_true = as.vector(apply(BiCopTau2Par(family = 5, get(paste0("tau_true_",test_case))), 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  
  ggplot(pred_cond) +
    geom_line(aes(obs, y)) +
    geom_line(aes(obs, theta_true), col = 2) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated theta') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30))
  
  
  pred_cond_mod = pred_cond %>%
    group_by(panel.name, obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .95), theta_q025 = quantile(y, .05)) 
  
  ggplot(pred_cond_mod) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated theta') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30))
  
  copula_uu_frank_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 5, par = pred_cond_mod$theta_mean)
  
  pred_cond_mod$U1 = copula_uu_frank_pred[,1]
  pred_cond_mod$U2 = copula_uu_frank_pred[,2]
  
  ggplot(pred_cond_mod) +
    geom_point(aes(U1, U2), size = 0.7) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('U1') +
    ylab('U2') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30))
  
  pred_cond_stat = pred_cond_mod %>%
    group_by(panel.name)%>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = pred_cond_stat %>%
    group_by(panel.name)%>%
    summarise_all(mean)
  
  xtable(pred_cond_summary , digits = 4)
}

# clayton

if(F){
  test_case = 5
  
  n.born.out.par = 250
  
  model.list.def <- list(
    get(paste0("clayton_mcmc_lb.def_IG22_",test_case)),
    get(paste0("clayton_mcmc_lb.def_IG11_",test_case)),
    get(paste0("clayton_mcmc_lb.def_LN01_",test_case)),
    get(paste0("clayton_mcmc_lb.def_LN05_",test_case)),
    get(paste0("clayton_mcmc_lossb.prior.high.gam_IG11_",test_case)),
    get(paste0("clayton_mcmc_lossb.prior.high.gam_IG22_",test_case)),
    get(paste0("clayton_mcmc_lossb.prior.high.gam_LN01_",test_case)),
    get(paste0("clayton_mcmc_lossb.prior.high.gam_LN05_",test_case)))
  
  names(model.list.def) <- c(
    'LB - default - IG22',
    'LB - default - IG11',
    'LB - default - LN01',
    'LB - default - LN05',
    'LB - high.gam - IG11',
    'LB - high.gam - IG22',
    'LB - high.gam - LN01',
    'LB - high.gam - LN05')
  
  
  # extract depth, number of terminal nodes, missing rate and loglik of all the trees
  depth.df <- apply_fun_models(fun_ = get_depth,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
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
    geom_vline(xintercept = seq(0,n.chain_par * (n.iter_par - n.born.out.par),by = 50),
               color = 'grey', size = 0.2, alpha=0.75)+
    geom_line(aes(x, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab(~n[L]) +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,n.chain_par * (n.iter_par - n.born.out.par),by = 50)) +
    theme(axis.text.x = element_text(angle = 30))
  
  trace.depth <- ggplot(depth.df) +
    geom_vline(xintercept = seq(0,n.chain_par * (n.iter_par - n.born.out.par),by = 50),
               color = 'grey', size = 0.2, alpha=0.75)+
    geom_line(aes(x, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab('Depth') +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,n.chain_par * (n.iter_par - n.born.out.par),by = 50)) +
    theme(axis.text.x = element_text(angle = 30))
  
  hist.nl
  hist.depth
  
  trace.nl
  trace.depth
  
  pred_cond = lapply(1:nrow(X_obs.norm), function(i)apply_fun_models(fun_ = function(x)get_value_tree(x, X_obs.norm[i,,drop = FALSE]),
                                                                     mcmc.list = model.list.def,
                                                                     born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par))
  
  
  pred_cond = do.call(rbind,pred_cond)
  
  pred_cond$obs = as.vector(apply(X_obs.norm, 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  pred_cond$theta_true = as.vector(apply(BiCopTau2Par(family = 3, get(paste0("tau_true_",test_case))), 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  
  ggplot(pred_cond) +
    geom_line(aes(obs, y)) +
    geom_line(aes(obs, theta_true), col = 2) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated theta') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30))
  
  
  pred_cond_mod = pred_cond %>%
    group_by(panel.name, obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .95), theta_q025 = quantile(y, .05)) 
  
  ggplot(pred_cond_mod) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated theta') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30))
  
  copula_uu_clayton_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 3, par = pred_cond_mod$theta_mean)
  
  pred_cond_mod$U1 = copula_uu_clayton_pred[,1]
  pred_cond_mod$U2 = copula_uu_clayton_pred[,2]
  
  ggplot(pred_cond_mod) +
    geom_point(aes(U1, U2), size = 0.7) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('U1') +
    ylab('U2') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30))
  
  pred_cond_stat = pred_cond_mod %>%
    group_by(panel.name)%>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = pred_cond_stat %>%
    group_by(panel.name)%>%
    summarise_all(mean)
  
  xtable(pred_cond_summary , digits = 4)
}
############################################################
# sanity check

if(F){
  library(GA)
  
  fun_opt = function(rho)loglik_frank(rho, copula_uu_gauss_2[,1], copula_uu_gauss_2[,2]) + log(dnorm(rho, sd = 10))
  
  mle_est = ga(type = "real-valued", fitness = fun_opt, lower = -35, upper = 35)
  
}