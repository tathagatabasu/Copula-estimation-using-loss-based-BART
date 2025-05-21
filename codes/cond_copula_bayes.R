# codes and packages
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
require(foreach)
require(parallel)
require(doParallel)

################################################################################
# data generation
################################################################################

if(T){
  set.seed(123)
  n <- 500
  X_obs <- matrix(runif(n), ncol = 1)
  
  # normalise predictors 
  X_obs.norm <- as.data.frame(apply(X_obs, 2, \(x) (x - min(x))/(max(x) - min(x))))
  X_obs.norm <- as.matrix(X_obs.norm)
  rownames(X_obs.norm) <- 1:nrow(X_obs)
  
  # tau with tree structure
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
  # non-convex
  tau_true_4 <- 0.6 - 0.3 * sin(2*X_obs) + 0.2 * sin(4*X_obs) + 0.3 * X_obs^2
  
  plot(X_obs, tau_true_1, xlab = "Observations", ylab = "tau")
  plot(X_obs, tau_true_2, xlab = "Observations", ylab = "tau")
  plot(X_obs, tau_true_3, xlab = "Observations", ylab = "tau")
  plot(X_obs, tau_true_4, xlab = "Observations", ylab = "tau")
  
  X_obs_pred = matrix(runif(n), ncol = 1)
  
  # normalise predictors 
  X_obs_pred.norm <- as.data.frame(apply(X_obs_pred, 2, \(x) (x - min(x))/(max(x) - min(x))))
  X_obs_pred.norm <- as.matrix(X_obs_pred.norm)
  rownames(X_obs_pred.norm) <- 1:nrow(X_obs_pred)
  
  # tau with tree structure
  tree_top <- generate_random_binary_tree_n_delta(8,4)
  tree_top <- assign_node_idx(tree_top)
  tree_top <- assign_split_rules(tree_top, X_obs_pred)
  tree_top <- assign_term_node_values(tree_top, 0.6, 0.1)
  tau_true_pred_1 <- sample_CART(tree_top, X_obs_pred, sigma_ = 0.001) 
  tau_true_pred_1 <- matrix(tau_true_1, ncol = 1)
  rm(tree_top)
  
  # monotone
  tau_true_pred_2 <- 0.3 + 0.2 * sin(3*X_obs_pred) + 0.3*X_obs_pred^2
  # convex
  tau_true_pred_3 <- 0.5 + 0.3 * sin(3*X_obs_pred)
  # non-convex
  tau_true_pred_4 <- 0.6 - 0.3 * sin(2*X_obs_pred) + 0.2 * sin(4*X_obs_pred) + 0.3 * X_obs_pred^2
  
  # mcmc params
  n.chain_par <- 1
  n.iter_par <- 15000
  n.born.out.par = 5000
  incl.split_par <- TRUE
  cont.unif_par <- TRUE
  moves.prob_par <- c(0.1, 0.4, 0.25, 0.25)
}

################################################################################
# gaussian
################################################################################
if(F){
  for (i in 1:4) {
    assign(paste0("copula_uu_gauss_",i), BiCopSim(n, family = 1, par = param_gauss(get(paste0("tau_true_",i)))))
  }
  
  lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) 
  
  for (i in 1:4) {
    assign(paste0("gauss_mcmc_lb.def_unif_",i), multichain_MCMC_copula(n.chain = n.chain_par,
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
  
  for (i in 1:4) {
    assign(paste0("gauss_mcmc_lb.def_jeff_",i), multichain_MCMC_copula(n.chain = n.chain_par,
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
  
  
  for (i in 1:4) {
    assign(paste0("gauss_mcmc_lb.def_half_",i), multichain_MCMC_copula(n.chain = n.chain_par,
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
  
  for (i in 1:4) {
    assign(paste0("gauss_mcmc_lb.def_two_",i), multichain_MCMC_copula(n.chain = n.chain_par,
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
  
}

# results

if(F){
  test_case = 1
  
  model.list.def <- list(
    get(paste0("gauss_mcmc_lb.def_unif_",test_case)),
    get(paste0("gauss_mcmc_lb.def_half_",test_case)),
    get(paste0("gauss_mcmc_lb.def_jeff_",test_case)),
    get(paste0("gauss_mcmc_lb.def_two_",test_case))
  )
  
  names(model.list.def) <- c(
    'LB - default - unif',
    'LB - default - half',
    'LB - default - jeff',
    'LB - default - two' 
  )
  
  # extract depth, number of terminal nodes and loglik of all the trees
  depth.df <- apply_fun_models(fun_ = get_depth,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  like.df <- apply_fun_models(fun_ = \(x) cart_log_lik_copula(tree_top = x, U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                              U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                              X = X_obs.norm, cop_type = "Gauss"), 
                              mcmc.list = model.list.def, 
                              born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  acc.df <- apply_fun_models_res(fun_ = \(x) (as.factor(x) == TRUE),
                                 mcmc.list = model.list.def,
                                 born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  depth.df_thin <- na.omit(depth.df[c(rep(NA,9), TRUE),])
  nterm.df_thin <- na.omit(nterm.df[c(rep(NA,9), TRUE),])
  like.df_thin <- na.omit(like.df[c(rep(NA,9), TRUE),])
  acc.df_thin <- na.omit(acc.df[c(rep(NA,9), TRUE),])
  
  df.sum.def <- data.frame(tree = nterm.df$x,
                           panel.name = nterm.df$panel.name,
                           loglik = like.df$y,
                           nterm = nterm.df$y)
  
  df.sum.def_thin <- data.frame(tree = nterm.df_thin$x,
                                panel.name = nterm.df_thin$panel.name,
                                loglik = like.df_thin$y,
                                nterm = nterm.df_thin$y,
                                depth = depth.df_thin$y,
                                acc_rate = acc.df_thin$y)
  
  hist.nl <- ggplot(nterm.df_thin) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab(~n[L]) +
    ylab('PMF') +
    theme_classic()
  
  hist.depth <- ggplot(depth.df_thin) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Depth') +
    ylab('PMF') +
    theme_classic()
  
  trace.nl <- ggplot(nterm.df_thin) +
    geom_line(aes(x/10, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab(~n[L]) +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,1000,by = 200))
  
  trace.depth <- ggplot(depth.df_thin) +
    geom_line(aes(x/10, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab('Depth') +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,1000,by = 200))
  
  trace.loglik <- ggplot(df.sum.def, aes(tree/10, loglik)) + 
    geom_line() + 
    facet_wrap(facets = ~panel.name) + 
    theme_classic() + 
    xlab('Iteration') + 
    ylab('Log likelihood') + 
    scale_x_continuous(breaks = seq(0,1000,by = 200))
  
  # plots
  hist.nl
  hist.depth
  trace.nl
  trace.depth
  trace.loglik
  
  # prediction
  
  pred_cond = lapply(1:nrow(X_obs_pred.norm), function(i)apply_fun_models(fun_ = function(x)get_value_tree(x, X_obs_pred.norm[i,,drop = FALSE]),
                                                                     mcmc.list = model.list.def,
                                                                     born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par))
  
  
  pred_cond = do.call(rbind,pred_cond)
  
  pred_cond$obs = as.vector(apply(X_obs_pred.norm, 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  pred_cond$theta_true = as.vector(apply(param_gauss(get(paste0("tau_true_pred_",test_case))), 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,9), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(panel.name, obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  copula_uu_gauss_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 1, par = pred_cond_mod$theta_mean)
  
  pred_cond_mod$U1 = copula_uu_gauss_pred[,1]
  pred_cond_mod$U2 = copula_uu_gauss_pred[,2]
  
  pred_cond_mod_p1 = pred_cond_mod %>% filter(panel.name == "LB - default - unif")
  pred_cond_mod_p2 = pred_cond_mod %>% filter(panel.name == "LB - default - half")
  pred_cond_mod_p3 = pred_cond_mod %>% filter(panel.name == "LB - default - jeff")
  pred_cond_mod_p4 = pred_cond_mod %>% filter(panel.name == "LB - default - two")
  
  # Plot 3D histogram
  
  hist_true <- hist2d(get(paste0("copula_uu_gauss_",i))[,1], get(paste0("copula_uu_gauss_",i))[,2], nbins = c(10,10), show = FALSE)
  
  hist_pred1 <- hist2d(pred_cond_mod_p1$U1, pred_cond_mod_p1$U2, nbins = c(10,10), show = FALSE)
  hist_pred2 <- hist2d(pred_cond_mod_p2$U1, pred_cond_mod_p2$U2, nbins = c(10,10), show = FALSE)
  hist_pred3 <- hist2d(pred_cond_mod_p3$U1, pred_cond_mod_p3$U2, nbins = c(10,10), show = FALSE)
  hist_pred4 <- hist2d(pred_cond_mod_p4$U1, pred_cond_mod_p4$U2, nbins = c(10,10), show = FALSE)
  
  par(mar = c(1,1,1,1), mfrow = c(2,3))
  
  hist3D(
    x = hist_true$x,
    y = hist_true$y,
    z = hist_true$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = "Observed copula"
  )
  
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
  
  pred_cond_stat = pred_cond_mod %>%
    group_by(panel.name)%>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = pred_cond_stat %>%
    group_by(panel.name)%>%
    summarise_all(mean)
  
  tree_stat_summary = df.sum.def_thin %>%
    group_by(panel.name)%>%
    summarise_at(c("nterm","depth","acc_rate"),mean)
  
  xtable(cbind(tree_stat_summary, pred_cond_summary[,-1]) , digits = 4)
}

################################################################################
# t
################################################################################
if(F){
  for (i in 1:4) {
    assign(paste0("copula_uu_t_",i), BiCopSim(n, family = 2, par = param_gauss(get(paste0("tau_true_",i))), par2 = 3))
  }
  
  lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) # c(1.5618883, 0.6293944)
  
  for (i in 1:4) {
    assign(paste0("t_mcmc_lb.def_unif_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                   n.iter = n.iter_par,
                                                                   X = X_obs.norm,
                                                                   U1 = get(paste0("copula_uu_t_",i))[,1],
                                                                   U2 = get(paste0("copula_uu_t_",i))[,2],
                                                                   mu = 0, cop_type = "t", 
                                                                   sigma = .5, alpha_val = 0, beta_val = 0, 
                                                                   log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
                                                                   prior_list = lb.prior.def, 
                                                                   moves.prob = moves.prob_par, 
                                                                   starting.tree = NULL,
                                                                   cont.unif = cont.unif_par,
                                                                   include.split = incl.split_par))
  }
  
  for (i in 1:4) {
    assign(paste0("t_mcmc_lb.def_jeff_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                   n.iter = n.iter_par,
                                                                   X = X_obs.norm,
                                                                   U1 = get(paste0("copula_uu_t_",i))[,1],
                                                                   U2 = get(paste0("copula_uu_t_",i))[,2],
                                                                   mu = 0, cop_type = "t", 
                                                                   sigma = .1, alpha_val = 0, beta_val = 0, 
                                                                   log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
                                                                   prior_list = lb.prior.def, 
                                                                   moves.prob = moves.prob_par, 
                                                                   starting.tree = NULL,
                                                                   cont.unif = cont.unif_par,
                                                                   include.split = incl.split_par))
  }
  
  
  for (i in 1:4) {
    assign(paste0("t_mcmc_lb.def_half_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                   n.iter = n.iter_par,
                                                                   X = X_obs.norm,
                                                                   U1 = get(paste0("copula_uu_t_",i))[,1],
                                                                   U2 = get(paste0("copula_uu_t_",i))[,2],
                                                                   mu = 0, cop_type = "t", 
                                                                   sigma = .1, alpha_val = .5, beta_val = .5, 
                                                                   log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
                                                                   prior_list = lb.prior.def, 
                                                                   moves.prob = moves.prob_par, 
                                                                   starting.tree = NULL,
                                                                   cont.unif = cont.unif_par,
                                                                   include.split = incl.split_par))
  }
  
  for (i in 1:4) {
    assign(paste0("t_mcmc_lb.def_two_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                  n.iter = n.iter_par,
                                                                  X = X_obs.norm,
                                                                  U1 = get(paste0("copula_uu_t_",i))[,1],
                                                                  U2 = get(paste0("copula_uu_t_",i))[,2],
                                                                  mu = 0, cop_type = "t", 
                                                                  sigma = .1, alpha_val = 2, beta_val = 2, 
                                                                  log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
                                                                  prior_list = lb.prior.def, 
                                                                  moves.prob = moves.prob_par, 
                                                                  starting.tree = NULL,
                                                                  cont.unif = cont.unif_par,
                                                                  include.split = incl.split_par))
  }
  
}

# results

if(F){
  test_case = 1
  
  model.list.def <- list(
    get(paste0("t_mcmc_lb.def_unif_",test_case)),
    get(paste0("t_mcmc_lb.def_half_",test_case)),
    get(paste0("t_mcmc_lb.def_jeff_",test_case)),
    get(paste0("t_mcmc_lb.def_two_",test_case))
  )
  
  names(model.list.def) <- c(
    'LB - default - unif',
    'LB - default - half',
    'LB - default - jeff',
    'LB - default - two'
  )
  
  # extract depth, number of terminal nodes and loglik of all the trees
  depth.df <- apply_fun_models(fun_ = get_depth,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  like.df <- apply_fun_models(fun_ = \(x) cart_log_lik_copula(tree_top = x, U1 = get(paste0("copula_uu_t_",i))[,1],
                                                              U2 = get(paste0("copula_uu_t_",i))[,2],
                                                              X = X_obs.norm, cop_type = "t"), 
                              mcmc.list = model.list.def, 
                              born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  acc.df <- apply_fun_models_res(fun_ = \(x) (as.factor(x) == TRUE),
                                 mcmc.list = model.list.def,
                                 born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  
  depth.df_thin <- na.omit(depth.df[c(rep(NA,9), TRUE),])
  nterm.df_thin <- na.omit(nterm.df[c(rep(NA,9), TRUE),])
  like.df_thin <- na.omit(like.df[c(rep(NA,9), TRUE),])
  acc.df_thin <- na.omit(acc.df[c(rep(NA,9), TRUE),])
  
  df.sum.def <- data.frame(tree = nterm.df$x,
                           panel.name = nterm.df$panel.name,
                           loglik = like.df$y,
                           nterm = nterm.df$y)
  
  df.sum.def_thin <- data.frame(tree = nterm.df_thin$x,
                                panel.name = nterm.df_thin$panel.name,
                                loglik = like.df_thin$y,
                                nterm = nterm.df_thin$y,
                                depth = depth.df_thin$y,
                                acc_rate = acc.df_thin$y)
  
  hist.nl <- ggplot(nterm.df_thin) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab(~n[L]) +
    ylab('PMF') +
    theme_classic()
  
  hist.depth <- ggplot(depth.df_thin) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Depth') +
    ylab('PMF') +
    theme_classic()
  
  trace.nl <- ggplot(nterm.df_thin) +
    geom_line(aes(x/10, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab(~n[L]) +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,1000,by = 200))
  
  trace.depth <- ggplot(depth.df_thin) +
    geom_line(aes(x/10, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab('Depth') +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,1000,by = 200))
  
  trace.loglik <- ggplot(df.sum.def, aes(tree/10, loglik)) + 
    geom_line() + 
    facet_wrap(facets = ~panel.name) + 
    theme_classic() + 
    xlab('Iteration') + 
    ylab('Log likelihood') + 
    scale_x_continuous(breaks = seq(0,1000,by = 200))
  
  # plots
  
  hist.nl
  hist.depth
  trace.nl
  trace.depth
  trace.loglik
  
  # prediction
  
  X_obs_pred = matrix(runif(n), ncol = 1)
  
  # normalise predictors 
  X_obs_pred.norm <- as.data.frame(apply(X_obs_pred, 2, \(x) (x - min(x))/(max(x) - min(x))))
  X_obs_pred.norm <- as.matrix(X_obs_pred.norm)
  rownames(X_obs_pred.norm) <- 1:nrow(X_obs_pred)
  
  pred_cond = lapply(1:nrow(X_obs_pred.norm), function(i)apply_fun_models(fun_ = function(x)get_value_tree(x, X_obs_pred.norm[i,,drop = FALSE]),
                                                                     mcmc.list = model.list.def,
                                                                     born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par))
  
  
  pred_cond = do.call(rbind,pred_cond)
  
  pred_cond$obs = as.vector(apply(X_obs_pred.norm, 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  pred_cond$theta_true = as.vector(apply(param_gauss(get(paste0("tau_true_",test_case))), 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,9), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(panel.name, obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  copula_uu_t_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 2, par = pred_cond_mod$theta_mean, par2 = 3)
  
  pred_cond_mod$U1 = copula_uu_t_pred[,1]
  pred_cond_mod$U2 = copula_uu_t_pred[,2]
  
  pred_cond_mod_p1 = pred_cond_mod %>% filter(panel.name == "LB - default - unif")
  pred_cond_mod_p2 = pred_cond_mod %>% filter(panel.name == "LB - default - half")
  pred_cond_mod_p3 = pred_cond_mod %>% filter(panel.name == "LB - default - jeff")
  pred_cond_mod_p4 = pred_cond_mod %>% filter(panel.name == "LB - default - two")
  
  # 3d histogram
  
  hist_true <- hist2d(get(paste0("copula_uu_t_",i))[,1], get(paste0("copula_uu_t_",i))[,2], nbins = c(10,10), show = FALSE)
  
  hist_pred1 <- hist2d(pred_cond_mod_p1$U1, pred_cond_mod_p1$U2, nbins = c(10,10), show = FALSE)
  hist_pred2 <- hist2d(pred_cond_mod_p2$U1, pred_cond_mod_p2$U2, nbins = c(10,10), show = FALSE)
  hist_pred3 <- hist2d(pred_cond_mod_p3$U1, pred_cond_mod_p3$U2, nbins = c(10,10), show = FALSE)
  hist_pred4 <- hist2d(pred_cond_mod_p4$U1, pred_cond_mod_p4$U2, nbins = c(10,10), show = FALSE)
  
  par(mar = c(1,1,1,1), mfrow = c(2,3))
  
  hist3D(
    x = hist_true$x,
    y = hist_true$y,
    z = hist_true$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = "Observed copula"
  )
  
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
  
  pred_cond_stat = pred_cond_mod %>%
    group_by(panel.name)%>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = pred_cond_stat %>%
    group_by(panel.name)%>%
    summarise_all(mean)
  
  tree_stat_summary = df.sum.def_thin %>%
    group_by(panel.name)%>%
    summarise_at(c("nterm","depth","acc_rate"),mean)
  
  xtable(cbind(tree_stat_summary, pred_cond_summary[,-1]) , digits = 4)
}

################################################################################
# Gumbel
################################################################################
if(F){
  for (i in 1:4) {
    assign(paste0("copula_uu_gumbel_",i), BiCopSim(n, family = 4, par = param_gumbel(get(paste0("tau_true_",i)))))
  }
  
  lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) 
  
  for (i in 1:4) {
    assign(paste0("gumbel_mcmc_lb.def_IG22_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                        n.iter = n.iter_par,
                                                                        X = X_obs.norm,
                                                                        U1 = get(paste0("copula_uu_gumbel_",i))[,1],
                                                                        U2 = get(paste0("copula_uu_gumbel_",i))[,2],
                                                                        mu = 0, cop_type = "gumbel", 
                                                                        sigma = .5, alpha_val = 2, beta_val = 2, 
                                                                        log_nor_mu = 0, log_nor_sigma = 1, prior_type = "IG",
                                                                        prior_list = lb.prior.def, 
                                                                        moves.prob = moves.prob_par, 
                                                                        starting.tree = NULL,
                                                                        cont.unif = cont.unif_par,
                                                                        include.split = incl.split_par))
  }
  
  for (i in 1:4) {
    assign(paste0("gumbel_mcmc_lb.def_IG11_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                        n.iter = n.iter_par,
                                                                        X = X_obs.norm,
                                                                        U1 = get(paste0("copula_uu_gumbel_",i))[,1],
                                                                        U2 = get(paste0("copula_uu_gumbel_",i))[,2],
                                                                        mu = 0, cop_type = "gumbel", 
                                                                        sigma = .5, alpha_val = 1, beta_val = 1, 
                                                                        log_nor_mu = 0, log_nor_sigma = 1, prior_type = "IG",
                                                                        prior_list = lb.prior.def, 
                                                                        moves.prob = moves.prob_par, 
                                                                        starting.tree = NULL,
                                                                        cont.unif = cont.unif_par,
                                                                        include.split = incl.split_par))
  }
  
  
  for (i in 1:4) {
    assign(paste0("gumbel_mcmc_lb.def_LN01_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                        n.iter = n.iter_par,
                                                                        X = X_obs.norm,
                                                                        U1 = get(paste0("copula_uu_gumbel_",i))[,1],
                                                                        U2 = get(paste0("copula_uu_gumbel_",i))[,2],
                                                                        mu = 0, cop_type = "gumbel", 
                                                                        sigma = .5, alpha_val = .5, beta_val = .5, 
                                                                        log_nor_mu = 0, log_nor_sigma = 1, prior_type = "LN",
                                                                        prior_list = lb.prior.def, 
                                                                        moves.prob = moves.prob_par, 
                                                                        starting.tree = NULL,
                                                                        cont.unif = cont.unif_par,
                                                                        include.split = incl.split_par))
  }
  
  for (i in 1:4) {
    assign(paste0("gumbel_mcmc_lb.def_LN05_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                        n.iter = n.iter_par,
                                                                        X = X_obs.norm,
                                                                        U1 = get(paste0("copula_uu_gumbel_",i))[,1],
                                                                        U2 = get(paste0("copula_uu_gumbel_",i))[,2],
                                                                        mu = 0, cop_type = "gumbel", 
                                                                        sigma = .5, alpha_val = 2, beta_val = 2, 
                                                                        log_nor_mu = 0, log_nor_sigma = 5, prior_type = "LN",
                                                                        prior_list = lb.prior.def, 
                                                                        moves.prob = moves.prob_par, 
                                                                        starting.tree = NULL,
                                                                        cont.unif = cont.unif_par,
                                                                        include.split = incl.split_par))
  }
}

# results

if(F){
  test_case = 1
  
  model.list.def <- list(
    get(paste0("gumbel_mcmc_lb.def_IG22_",test_case)),
    get(paste0("gumbel_mcmc_lb.def_IG11_",test_case)),
    get(paste0("gumbel_mcmc_lb.def_LN01_",test_case)),
    get(paste0("gumbel_mcmc_lb.def_LN05_",test_case))
  )
  
  names(model.list.def) <- c(
    'LB - default - IG22',
    'LB - default - IG11',
    'LB - default - LN01',
    'LB - default - LN05'
  )
  
  # extract depth, number of terminal nodes, missing rate and loglik of all the trees
  depth.df <- apply_fun_models(fun_ = get_depth,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  like.df <- apply_fun_models(fun_ = \(x) cart_log_lik_copula(tree_top = x, U1 = get(paste0("copula_uu_gumbel_",i))[,1],
                                                              U2 = get(paste0("copula_uu_gumbel_",i))[,2],
                                                              X = X_obs.norm, cop_type = "gumbel"), 
                              mcmc.list = model.list.def, 
                              born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  acc.df <- apply_fun_models_res(fun_ = \(x) (as.factor(x) == TRUE),
                                 mcmc.list = model.list.def,
                                 born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  
  depth.df_thin <- na.omit(depth.df[c(rep(NA,9), TRUE),])
  nterm.df_thin <- na.omit(nterm.df[c(rep(NA,9), TRUE),])
  like.df_thin <- na.omit(like.df[c(rep(NA,9), TRUE),])
  acc.df_thin <- na.omit(acc.df[c(rep(NA,9), TRUE),])
  
  df.sum.def <- data.frame(tree = nterm.df$x,
                           panel.name = nterm.df$panel.name,
                           loglik = like.df$y,
                           nterm = nterm.df$y)
  
  df.sum.def_thin <- data.frame(tree = nterm.df_thin$x,
                                panel.name = nterm.df_thin$panel.name,
                                loglik = like.df_thin$y,
                                nterm = nterm.df_thin$y,
                                depth = depth.df_thin$y,
                                acc_rate = acc.df_thin$y)
  
  hist.nl <- ggplot(nterm.df_thin) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab(~n[L]) +
    ylab('PMF') +
    theme_classic()
  
  hist.depth <- ggplot(depth.df_thin) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Depth') +
    ylab('PMF') +
    theme_classic()
  
  trace.nl <- ggplot(nterm.df_thin) +
    geom_line(aes(x/10, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab(~n[L]) +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,1000,by = 200))
  
  trace.depth <- ggplot(depth.df_thin) +
    geom_line(aes(x/10, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab('Depth') +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,1000,by = 200))
  
  trace.loglik <- ggplot(df.sum.def, aes(tree/10, loglik)) + 
    geom_line() + 
    facet_wrap(facets = ~panel.name) + 
    theme_classic() + 
    xlab('Iteration') + 
    ylab('Log likelihood') + 
    scale_x_continuous(breaks = seq(0,1000,by = 200))
  
  # plots
  
  hist.nl
  hist.depth
  trace.nl
  trace.depth
  trace.loglik
  
  # prediction
  
  pred_cond = lapply(1:nrow(X_obs_pred.norm), function(i)apply_fun_models(fun_ = function(x)get_value_tree(x, X_obs_pred.norm[i,,drop = FALSE]),
                                                                     mcmc.list = model.list.def,
                                                                     born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par))
  
  pred_cond = do.call(rbind,pred_cond)
  
  pred_cond$obs = as.vector(apply(X_obs_pred.norm, 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  pred_cond$theta_true = as.vector(apply(param_gumbel(get(paste0("tau_true_",test_case))), 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,9), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(panel.name, obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  copula_uu_gumbel_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 4, par = pmax(pred_cond_mod$theta_mean,1))
  
  pred_cond_mod$U1 = copula_uu_gumbel_pred[,1]
  pred_cond_mod$U2 = copula_uu_gumbel_pred[,2]
  
  # Plot 3D histogram
  hist_true <- hist2d(get(paste0("copula_uu_gumbel_",i))[,1], get(paste0("copula_uu_gumbel_",i))[,2], nbins = c(10,10), show = FALSE)
  
  hist_pred1 <- hist2d(pred_cond_mod_p1$U1, pred_cond_mod_p1$U2, nbins = c(10,10), show = FALSE)
  hist_pred2 <- hist2d(pred_cond_mod_p2$U1, pred_cond_mod_p2$U2, nbins = c(10,10), show = FALSE)
  hist_pred3 <- hist2d(pred_cond_mod_p3$U1, pred_cond_mod_p3$U2, nbins = c(10,10), show = FALSE)
  hist_pred4 <- hist2d(pred_cond_mod_p4$U1, pred_cond_mod_p4$U2, nbins = c(10,10), show = FALSE)
  
  par(mar = c(1,1,1,1), mfrow = c(2,3))
  
  hist3D(
    x = hist_true$x,
    y = hist_true$y,
    z = hist_true$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = "Observed copula"
  )
  
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
  
  pred_cond_stat = pred_cond_mod %>%
    group_by(panel.name)%>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = pred_cond_stat %>%
    group_by(panel.name)%>%
    summarise_all(mean)
  
  tree_stat_summary = df.sum.def_thin %>%
    group_by(panel.name)%>%
    summarise_at(c("nterm","depth","acc_rate"),mean)
  
  xtable(cbind(tree_stat_summary, pred_cond_summary[,-1]) , digits = 4)
}

################################################################################
# clayton
################################################################################
if(F){
  for (i in 1:4) {
    assign(paste0("copula_uu_clayton_",i), BiCopSim(n, family = 3, par = param_clayton(get(paste0("tau_true_",i)))))
  }
  
  lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) 
  
  for (i in 1:4) {
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
  
  for (i in 1:4) {
    assign(paste0("clayton_mcmc_lb.def_IG11_",i), multichain_MCMC_copula(n.chain = n.chain_par,
                                                                         n.iter = n.iter_par,
                                                                         X = X_obs.norm,
                                                                         U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                         U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                         mu = 0, cop_type = "Clayton", 
                                                                         sigma = .5, alpha_val = 1, beta_val = 1, 
                                                                         log_nor_mu = 0, log_nor_sigma = 1, prior_type = "IG",
                                                                         prior_list = lb.prior.def, 
                                                                         moves.prob = moves.prob_par, 
                                                                         starting.tree = NULL,
                                                                         cont.unif = cont.unif_par,
                                                                         include.split = incl.split_par))
  }
  
  
  for (i in 1:4) {
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
  
  for (i in 1:4) {
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
  
}

# results

if(F){
  test_case = 1
  
  model.list.def <- list(
    get(paste0("clayton_mcmc_lb.def_IG22_",test_case)),
    get(paste0("clayton_mcmc_lb.def_IG11_",test_case)),
    get(paste0("clayton_mcmc_lb.def_LN01_",test_case)),
    get(paste0("clayton_mcmc_lb.def_LN05_",test_case))
    )
  
  names(model.list.def) <- c(
    'LB - default - IG22',
    'LB - default - IG11',
    'LB - default - LN01',
    'LB - default - LN05'
    )
  
  # extract depth, number of terminal nodes, missing rate and loglik of all the trees
  depth.df <- apply_fun_models(fun_ = get_depth,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  like.df <- apply_fun_models(fun_ = \(x) cart_log_lik_copula(tree_top = x, U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                              U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                              X = X_obs.norm, cop_type = "Clayton"), 
                              mcmc.list = model.list.def, 
                              born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  acc.df <- apply_fun_models_res(fun_ = \(x) (as.factor(x) == TRUE),
                                 mcmc.list = model.list.def,
                                 born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  
  depth.df_thin <- na.omit(depth.df[c(rep(NA,9), TRUE),])
  nterm.df_thin <- na.omit(nterm.df[c(rep(NA,9), TRUE),])
  like.df_thin <- na.omit(like.df[c(rep(NA,9), TRUE),])
  acc.df_thin <- na.omit(acc.df[c(rep(NA,9), TRUE),])
  
  df.sum.def <- data.frame(tree = nterm.df$x,
                           panel.name = nterm.df$panel.name,
                           loglik = like.df$y,
                           nterm = nterm.df$y)
  
  df.sum.def_thin <- data.frame(tree = nterm.df_thin$x,
                                panel.name = nterm.df_thin$panel.name,
                                loglik = like.df_thin$y,
                                nterm = nterm.df_thin$y,
                                depth = depth.df_thin$y,
                                acc_rate = acc.df_thin$y)
  
  hist.nl <- ggplot(nterm.df_thin) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab(~n[L]) +
    ylab('PMF') +
    theme_classic()
  
  hist.depth <- ggplot(depth.df_thin) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Depth') +
    ylab('PMF') +
    theme_classic()
  
  trace.nl <- ggplot(nterm.df_thin) +
    geom_line(aes(x/10, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab(~n[L]) +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,1000,by = 200))
  
  trace.depth <- ggplot(depth.df_thin) +
    geom_line(aes(x/10, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab('Depth') +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,1000,by = 200))
  
  trace.loglik <- ggplot(df.sum.def, aes(tree/10, loglik)) + 
    geom_line() + 
    facet_wrap(facets = ~panel.name) + 
    theme_classic() + 
    xlab('Iteration') + 
    ylab('Log likelihood') + 
    scale_x_continuous(breaks = seq(0,1000,by = 200))
  
  # plots
  
  hist.nl
  hist.depth
  trace.nl
  trace.depth
  trace.loglik
  
  # prediction
  
  X_obs_pred = matrix(runif(n), ncol = 1)
  
  # normalise predictors 
  X_obs_pred.norm <- as.data.frame(apply(X_obs_pred, 2, \(x) (x - min(x))/(max(x) - min(x))))
  X_obs_pred.norm <- as.matrix(X_obs_pred.norm)
  rownames(X_obs_pred.norm) <- 1:nrow(X_obs_pred)
  
  pred_cond = lapply(1:nrow(X_obs_pred.norm), function(i)apply_fun_models(fun_ = function(x)get_value_tree(x, X_obs_pred.norm[i,,drop = FALSE]),
                                                                     mcmc.list = model.list.def,
                                                                     born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par))
  
  pred_cond = do.call(rbind,pred_cond)
  
  pred_cond$obs = as.vector(apply(X_obs_pred.norm, 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  pred_cond$theta_true = as.vector(apply(param_clayton(get(paste0("tau_true_",test_case))), 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,9), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(panel.name, obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  copula_uu_clayton_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 3, par = pred_cond_mod$theta_mean)
  
  pred_cond_mod$U1 = copula_uu_clayton_pred[,1]
  pred_cond_mod$U2 = copula_uu_clayton_pred[,2]
  
  # Plot 3D histogram
  hist_true <- hist2d(copula_uu_clayton_1[,1], copula_uu_clayton_1[,2], nbins = c(10,10), show = FALSE)
  
  hist_pred1 <- hist2d(pred_cond_mod_p1$U1, pred_cond_mod_p1$U2, nbins = c(10,10), show = FALSE)
  hist_pred2 <- hist2d(pred_cond_mod_p2$U1, pred_cond_mod_p2$U2, nbins = c(10,10), show = FALSE)
  hist_pred3 <- hist2d(pred_cond_mod_p3$U1, pred_cond_mod_p3$U2, nbins = c(10,10), show = FALSE)
  hist_pred4 <- hist2d(pred_cond_mod_p4$U1, pred_cond_mod_p4$U2, nbins = c(10,10), show = FALSE)
  
  par(mar = c(1,1,1,1), mfrow = c(2,3))
  
  hist3D(
    x = hist_true$x,
    y = hist_true$y,
    z = hist_true$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = "Observed copula"
  )
  
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
  
  pred_cond_stat = pred_cond_mod %>%
    group_by(panel.name)%>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = pred_cond_stat %>%
    group_by(panel.name)%>%
    summarise_all(mean)
  
  tree_stat_summary = df.sum.def_thin %>%
    group_by(panel.name)%>%
    summarise_at(c("nterm","depth","acc_rate"),mean)
  
  xtable(cbind(tree_stat_summary, pred_cond_summary[,-1]) , digits = 4)
}
