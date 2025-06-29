# codes and packages
source('code/import_functions.R')
source('mclapply.R')
# source('MCMC_BART_copula.R')
source('test_MCMC_copula_mult_tree.R')
library(data.tree)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(VineCopula)
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
set.seed(1e5)

# load("analysis_sim_dat_new.RData")

if(F){
  n <- 500
  X_obs <- matrix(runif(n), ncol = 1)
  
  # normalise predictors 
  X_obs.norm <- as.data.frame(apply(X_obs, 2, \(x) (x - min(x))/(max(x) - min(x))))
  X_obs.norm <- as.matrix(X_obs.norm)
  rownames(X_obs.norm) <- 1:nrow(X_obs)
  
  # tau with tree structure
  # initialize list
  tree_ex <- list()
  # one split on the left branch
  tree_ex$left <- list(left = list(left = NULL,
                                   right = NULL),
                       right = list(left = NULL,
                                    right = NULL))
  # no split on the right branch
  tree_ex$right <-  list(left = NULL,
                         right = NULL)
  
  # assign node index
  tree_ex <- assign_node_idx(tree_ex)
  # plot the tree
  get_tree_plot.idx(tree_ex)
  
  # assign first splitting rule
  tree_ex$cond <- list(x.idx = 1, x.val = 0.6)
  # assign second splitting rule
  tree_ex$left$cond <- list(x.idx = 1, x.val = 0.25)
  # set mu_3
  tree_ex$left$left$cond <- NULL
  tree_ex$left$left$value <- .5
  # set mu_4
  tree_ex$left$right$cond <- NULL
  tree_ex$left$right$value <- .7
  # set mu_5
  tree_ex$right$cond <- NULL
  tree_ex$right$value <- .3
  # plot the tree
  get_tree_plot(tree_ex)
  tau_true_1 <- sample_CART(tree_ex, X_obs, sigma_ = 0.001) 
  tau_true_1 <- matrix(tau_true_1, ncol = 1)
  
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
  
  # gauss
  
  for (i in 1:4) {
    assign(paste0("copula_uu_gauss_",i), BiCopSim(n, family = 1, par = param_gauss(get(paste0("tau_true_",i)))))
  }
  
  # t
  
  for (i in 1:4) {
    assign(paste0("copula_uu_t_",i), BiCopSim(n, family = 2, par = param_gauss(get(paste0("tau_true_",i))), par2 = 3))
  }
  
  # gumbel 
  
  for (i in 1:4) {
    assign(paste0("copula_uu_gumbel_",i), BiCopSim(n, family = 4, par = param_gumbel(get(paste0("tau_true_",i)))))
  }
  
  # clayton
  
  for (i in 1:4) {
    assign(paste0("copula_uu_clayton_",i), BiCopSim(n, family = 3, par = param_clayton(get(paste0("tau_true_",i)))))
  }
  
  # dataset for prediction
  
  X_obs_pred = matrix(runif(n), ncol = 1)
  
  # normalise predictors 
  X_obs_pred.norm <- as.data.frame(apply(X_obs_pred, 2, \(x) (x - min(x))/(max(x) - min(x))))
  X_obs_pred.norm <- as.matrix(X_obs_pred.norm)
  rownames(X_obs_pred.norm) <- 1:nrow(X_obs_pred)
  
  # tau with tree structure
  tau_true_pred_1 <- sample_CART(tree_ex, X_obs_pred, sigma_ = 0.001) 
  tau_true_pred_1 <- matrix(tau_true_pred_1, ncol = 1)
  
  # monotone
  tau_true_pred_2 <- 0.3 + 0.2 * sin(3*X_obs_pred) + 0.3*X_obs_pred^2
  # convex
  tau_true_pred_3 <- 0.5 + 0.3 * sin(3*X_obs_pred)
  # non-convex
  tau_true_pred_4 <- 0.6 - 0.3 * sin(2*X_obs_pred) + 0.2 * sin(4*X_obs_pred) + 0.3 * X_obs_pred^2
  
  # mcmc params
  n.chain_par <- 1
  n.iter_par <- 6000
  n.born.out.par <- 1000
  n.thin <- 1
  incl.split_par <- TRUE
  cont.unif_par <- TRUE
  moves.prob_par <- c(0.1, 0.3, 0.3, 0.3)
  lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) 
  
}

################################################################################
# gaussian
################################################################################
if(T){
  for (i in 3:4) {
    assign(paste0("gauss_mcmc_lb.def_single_",i), MCMC_copula(n.iter = n.iter_par,
                                                                  n.tree = 1,
                                                                  X = X_obs.norm,
                                                                  U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                                  U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                                  prior_list = lb.prior.def, 
                                                                  moves.prob = moves.prob_par, 
                                                                  starting.tree = NULL,
                                                                  cont.unif = cont.unif_par,
                                                                  include.split = incl.split_par,
                                                                  prop_mu = 0, prop_sigma = 1,
                                                                  theta_param_1 = 0, theta_param_2 = 1,
                                                                  prior_type = "N",
                                                                  cop_type = "gauss"))
    
    cat('done case', i, '\n')
  }
  
}

# save(gauss_mcmc_lb.def_single_4, file = "gauss_mcmc_lb.def_single_4.Rdata")
# rm(gauss_mcmc_lb.def_single_1,gauss_mcmc_lb.def_single_2,gauss_mcmc_lb.def_single_3,gauss_mcmc_lb.def_single_4)

# results

if(F){
  test_case = 1
  
  load(paste0("gauss_mcmc_lb.def_single_",test_case,".Rdata"))
  
  model <- get(paste0("gauss_mcmc_lb.def_single_",test_case))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs_pred.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  n.thin <- 10
  
  pred_val_vec = as.vector(pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  pred_obs = rep(X_obs_pred.norm, each = (n.chain_par * (n.iter_par - n.born.out.par)))
  
  theta_true = rep(param_gauss(get(paste0("tau_true_pred_",test_case))), each = (n.chain_par * (n.iter_par - n.born.out.par)))
  
  pred_cond <- data.frame("obs" = pred_obs)
  pred_cond$obs = pred_obs
  pred_cond$theta_true = theta_true
  pred_cond$y = link_gauss(pred_val_vec)
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,(n.thin-1)), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  ggplot(pred_cond_mod) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    # facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated rho') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod %>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_val <- apply(pred_val, 1, function(x)loglik_gauss(link_gauss(x), get(paste0("copula_uu_gauss_",i))[,1], get(paste0("copula_uu_gauss_",i))[,2]))
  
  like_df <-data.frame("nn" = like_val)
  like_df$idx <- 1:(n.chain_par*n.iter_par)
  
  pl_like <- ggplot(like_df, aes(idx, nn)) + 
    geom_line() + 
    ylab('log-likelihood') +
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  pl_like
  
  # nterm
  
  nt_lb.df <- nterm_BART(model)
  
  pl_nl <- ggplot(nt_lb.df, aes(idx, nn)) + 
    geom_line() + 
    ylab('nterm') + 
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  pl_nl
  
  # depth
  
  depth_lb.df <- depth_BART(model)
  
  pl_dp <- ggplot(depth_lb.df, aes(idx, nn)) + 
    geom_line() + 
    ylab('depth') +
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  pl_dp
  
  # acceptance
  
  acc_lb.df <- acc_BART(model)
  
  tree_stat_summary = c(mean(nt_lb.df$nn), mean(depth_lb.df$nn), mean(acc_lb.df$nn==TRUE))
  
  names(tree_stat_summary) <- c("nterm", "depth", "acc")
  
  xtable(t(as.matrix(c(tree_stat_summary, pred_cond_summary))), digits = 4)
  
  conv_diag_sum <- cbind(conv_diag(depth_lb.df,n.born.out.par,n.thin), conv_diag(nt_lb.df, n.born.out.par,n.thin), conv_diag(like_df,n.born.out.par,n.thin))
  xtable(conv_diag_sum)
  
  
  if(F){
    copula_uu_gauss_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 1, par = pred_cond_mod$theta_mean)
    
    pred_cond_mod$U1 = copula_uu_gauss_pred[,1]
    pred_cond_mod$U2 = copula_uu_gauss_pred[,2]
    
    # Plot 3D histogram
    
    hist_true <- hist2d(get(paste0("copula_uu_gauss_",i))[,1], get(paste0("copula_uu_gauss_",i))[,2], nbins = c(10,10), show = FALSE)
    
    hist_pred <- hist2d(pred_cond_mod$U1, pred_cond_mod$U2, nbins = c(10,10), show = FALSE)
    
    par(mar = c(1,1,1,1), mfrow = c(1,2))
    
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
      x = hist_pred$x,
      y = hist_pred$y,
      z = hist_pred$counts,
      colvar = NULL,  # disables color mapping
      col = "lightblue",  # solid monochrome color
      border = "black",
      theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
      lighting = TRUE,
      ltheta = 120, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "",
      main = unique(pred_cond_mod$panel.name)
    )
  }
}

################################################################################
# t
################################################################################
if(T){
  
  for (i in 2:3) {
    assign(paste0("t_mcmc_lb.def_single_",i), MCMC_copula(n.iter = n.iter_par,
                                                              n.tree = 1,
                                                              X = X_obs.norm,
                                                              U1 = get(paste0("copula_uu_t_",i))[,1],
                                                              U2 = get(paste0("copula_uu_t_",i))[,2],
                                                          prior_list = lb.prior.def, 
                                                          moves.prob = moves.prob_par, 
                                                          starting.tree = NULL,
                                                          cont.unif = cont.unif_par,
                                                          include.split = incl.split_par,
                                                          prop_mu = 0, prop_sigma = 1,
                                                          theta_param_1 = 0, theta_param_2 = 1,
                                                          prior_type = "N",
                                                          cop_type = "t"))
  }
  
}

# save(t_mcmc_lb.def_single_4, file = "t_mcmc_lb.def_single_4.Rdata")
# rm(t_mcmc_lb.def_single_1,t_mcmc_lb.def_single_2,t_mcmc_lb.def_single_3,t_mcmc_lb.def_single_4)

# results

if(F){
  test_case = 4
  
  load(paste0("t_mcmc_lb.def_single_",test_case,".Rdata"))
  
  model <- get(paste0("t_mcmc_lb.def_single_",test_case))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs_pred.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  n.thin <- 10
  
  pred_val_vec = as.vector(pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  pred_obs = rep(X_obs_pred.norm, each = (n.chain_par * (n.iter_par - n.born.out.par)))
  
  theta_true = rep(param_t(get(paste0("tau_true_pred_",test_case))), each = (n.chain_par * (n.iter_par - n.born.out.par)))
  
  pred_cond <- data.frame("obs" = pred_obs)
  pred_cond$obs = pred_obs
  pred_cond$theta_true = theta_true
  pred_cond$y = link_t(pred_val_vec)
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,(n.thin-1)), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  ggplot(pred_cond_mod) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    # facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated rho') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod %>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_val <- apply(pred_val, 1, function(x)loglik_t(link_t(x), get(paste0("copula_uu_t_",i))[,1], get(paste0("copula_uu_t_",i))[,2]))
  
  like_df <-data.frame("nn" = like_val)
  like_df$idx <- 1:(n.chain_par*n.iter_par)
  
  pl_like <- ggplot(like_df, aes(idx, nn)) + 
    geom_line() + 
    ylab('log-likelihood') +
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  pl_like
  
  # nterm
  
  nt_lb.df <- nterm_BART(model)
  
  pl_nl <- ggplot(nt_lb.df, aes(idx, nn)) + 
    geom_line() + 
    ylab('nterm') + 
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  pl_nl
  
  # depth
  
  depth_lb.df <- depth_BART(model)
  
  pl_dp <- ggplot(depth_lb.df, aes(idx, nn)) + 
    geom_line() + 
    ylab('depth') +
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  pl_dp
  
  # acceptance
  
  acc_lb.df <- acc_BART(model)
  
  tree_stat_summary = c(mean(nt_lb.df$nn), mean(depth_lb.df$nn), mean(acc_lb.df$nn==TRUE))
  
  names(tree_stat_summary) <- c("nterm", "depth", "acc")
  
  xtable(t(as.matrix(c(tree_stat_summary, pred_cond_summary))), digits = 4)
  
  conv_diag_sum <- cbind(conv_diag(depth_lb.df,n.born.out.par,n.thin), conv_diag(nt_lb.df, n.born.out.par,n.thin), conv_diag(like_df,n.born.out.par,n.thin))
  xtable(conv_diag_sum)
  
  
  
  if(F){
    copula_uu_t_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 2, par = pred_cond_mod$theta_mean)
    
    pred_cond_mod$U1 = copula_uu_t_pred[,1]
    pred_cond_mod$U2 = copula_uu_t_pred[,2]
    
    # Plot 3D histogram
    
    hist_true <- hist2d(get(paste0("copula_uu_t_",i))[,1], get(paste0("copula_uu_t_",i))[,2], nbins = c(10,10), show = FALSE)
    
    hist_pred <- hist2d(pred_cond_mod$U1, pred_cond_mod$U2, nbins = c(10,10), show = FALSE)
    
    par(mar = c(1,1,1,1), mfrow = c(1,2))
    
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
      x = hist_pred$x,
      y = hist_pred$y,
      z = hist_pred$counts,
      colvar = NULL,  # disables color mapping
      col = "lightblue",  # solid monochrome color
      border = "black",
      theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
      lighting = TRUE,
      ltheta = 120, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "",
      main = unique(pred_cond_mod$panel.name)
    )
  }
}

################################################################################
# Gumbel
################################################################################
if(T){
  
  for (i in c(1,3:4)) {
    assign(paste0("gumbel_mcmc_lb.def_single_",i), MCMC_copula(n.iter = n.iter_par,
                                                                   n.tree = 1,
                                                                   X = X_obs.norm,
                                                                   U1 = get(paste0("copula_uu_gumbel_",i))[,1],
                                                                   U2 = get(paste0("copula_uu_gumbel_",i))[,2],
                                                               prior_list = lb.prior.def, 
                                                               moves.prob = moves.prob_par, 
                                                               starting.tree = NULL,
                                                               cont.unif = cont.unif_par,
                                                               include.split = incl.split_par,
                                                               prop_mu = 0, prop_sigma = 1,
                                                               theta_param_1 = 0, theta_param_2 = 1,
                                                               prior_type = "N",
                                                               cop_type = "gumbel"))
  }
  
}

# save(gumbel_mcmc_lb.def_single_4, file = "gumbel_mcmc_lb.def_single_4.Rdata")
# rm(gumbel_mcmc_lb.def_single_1,gumbel_mcmc_lb.def_single_2,gumbel_mcmc_lb.def_single_3,gumbel_mcmc_lb.def_single_4)

# results

if(F){
  test_case = 1
  
  load(paste0("gumbel_mcmc_lb.def_single_",test_case,".Rdata"))
  
  model <- get(paste0("gumbel_mcmc_lb.def_single_",test_case))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs_pred.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  n.thin <- 10
  
  pred_val_vec = as.vector(pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  pred_obs = rep(X_obs_pred.norm, each = (n.chain_par * (n.iter_par - n.born.out.par)))
  
  theta_true = rep(param_gumbel(get(paste0("tau_true_pred_",test_case))), each = (n.chain_par * (n.iter_par - n.born.out.par)))
  
  pred_cond <- data.frame("obs" = pred_obs)
  pred_cond$obs = pred_obs
  pred_cond$theta_true = theta_true
  pred_cond$y = link_gumbel(pred_val_vec)
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,(n.thin-1)), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  ggplot(pred_cond_mod) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    # facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated rho') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod %>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_val <- apply(pred_val, 1, function(x)loglik_gumbel(link_gumbel(x), get(paste0("copula_uu_gumbel_",i))[,1], get(paste0("copula_uu_gumbel_",i))[,2]))
  
  like_df <-data.frame("nn" = like_val)
  like_df$idx <- 1:(n.chain_par*n.iter_par)
  
  pl_like <- ggplot(like_df, aes(idx, nn)) + 
    geom_line() + 
    ylab('log-likelihood') +
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  pl_like
  
  # nterm
  
  nt_lb.df <- nterm_BART(model)
  
  pl_nl <- ggplot(nt_lb.df, aes(idx, nn)) + 
    geom_line() + 
    ylab('nterm') + 
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  pl_nl
  
  # depth
  
  depth_lb.df <- depth_BART(model)
  
  pl_dp <- ggplot(depth_lb.df, aes(idx, nn)) + 
    geom_line() + 
    ylab('depth') +
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  pl_dp
  
  # acceptance
  
  acc_lb.df <- acc_BART(model)
  
  tree_stat_summary = c(mean(nt_lb.df$nn), mean(depth_lb.df$nn), mean(acc_lb.df$nn==TRUE))
  
  names(tree_stat_summary) <- c("nterm", "depth", "acc")
  
  xtable(t(as.matrix(c(tree_stat_summary, pred_cond_summary))), digits = 4)
  
  conv_diag_sum <- cbind(conv_diag(depth_lb.df,n.born.out.par,n.thin), conv_diag(nt_lb.df, n.born.out.par,n.thin), conv_diag(like_df,n.born.out.par,n.thin))
  xtable(conv_diag_sum)
  
  
  
  if(F){
    copula_uu_gumbel_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 1, par = pred_cond_mod$theta_mean)
    
    pred_cond_mod$U1 = copula_uu_gumbel_pred[,1]
    pred_cond_mod$U2 = copula_uu_gumbel_pred[,2]
    
    # Plot 3D histogram
    
    hist_true <- hist2d(get(paste0("copula_uu_gumbel_",i))[,1], get(paste0("copula_uu_gumbel_",i))[,2], nbins = c(10,10), show = FALSE)
    
    hist_pred <- hist2d(pred_cond_mod$U1, pred_cond_mod$U2, nbins = c(10,10), show = FALSE)
    
    par(mar = c(1,1,1,1), mfrow = c(1,2))
    
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
      x = hist_pred$x,
      y = hist_pred$y,
      z = hist_pred$counts,
      colvar = NULL,  # disables color mapping
      col = "lightblue",  # solid monochrome color
      border = "black",
      theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
      lighting = TRUE,
      ltheta = 120, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "",
      main = unique(pred_cond_mod$panel.name)
    )
  }
}

################################################################################
# clayton
################################################################################
if(T){
  
  for (i in 2:4) {
    assign(paste0("clayton_mcmc_lb.def_single_",i), MCMC_copula(n.iter = n.iter_par,
                                                                    n.tree = 1,
                                                                    X = X_obs.norm,
                                                                    U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                    U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                prior_list = lb.prior.def, 
                                                                moves.prob = moves.prob_par, 
                                                                starting.tree = NULL,
                                                                cont.unif = cont.unif_par,
                                                                include.split = incl.split_par,
                                                                prop_mu = 0, prop_sigma = 1,
                                                                theta_param_1 = 0, theta_param_2 = 1,
                                                                prior_type = "N",
                                                                cop_type = "clayton"))
  }
  
  # all prop_mu = 0, prop_sigma = .2,
  # 3 prop_mu = 0, prop_sigma = 3,
  
}
# save(clayton_mcmc_lb.def_single_4, file = "clayton_mcmc_lb.def_single_4.Rdata")
# rm(clayton_mcmc_lb.def_single_1,clayton_mcmc_lb.def_single_2,clayton_mcmc_lb.def_single_3,clayton_mcmc_lb.def_single_4)

# results

if(F){
  test_case = 4
  
  load(paste0("clayton_mcmc_lb.def_single_",test_case,".Rdata"))
  
  model <- get(paste0("clayton_mcmc_lb.def_single_",test_case))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs_pred.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  n.thin <- 10
  
  pred_val_vec = as.vector(pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  pred_obs = rep(X_obs_pred.norm, each = (n.chain_par * (n.iter_par - n.born.out.par)))
  
  theta_true = rep(param_clayton(get(paste0("tau_true_pred_",test_case))), each = (n.chain_par * (n.iter_par - n.born.out.par)))
  
  pred_cond <- data.frame("obs" = pred_obs)
  pred_cond$obs = pred_obs
  pred_cond$theta_true = theta_true
  pred_cond$y = link_clayton(pred_val_vec)
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,(n.thin-1)), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  ggplot(pred_cond_mod) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    # facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated rho') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod %>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_val <- apply(pred_val, 1, function(x)loglik_clayton(link_clayton(x), get(paste0("copula_uu_clayton_",i))[,1], get(paste0("copula_uu_clayton_",i))[,2]))
  
  like_df <-data.frame("nn" = like_val)
  like_df$idx <- 1:(n.chain_par*n.iter_par)
  
  pl_like <- ggplot(like_df, aes(idx, nn)) + 
    geom_line() + 
    ylab('log-likelihood') +
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  pl_like
  
  # nterm
  
  nt_lb.df <- nterm_BART(model)
  
  pl_nl <- ggplot(nt_lb.df, aes(idx, nn)) + 
    geom_line() + 
    ylab('nterm') + 
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  pl_nl
  
  # depth
  
  depth_lb.df <- depth_BART(model)
  
  pl_dp <- ggplot(depth_lb.df, aes(idx, nn)) + 
    geom_line() + 
    ylab('depth') +
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  pl_dp
  
  # acceptance
  
  acc_lb.df <- acc_BART(model)
  
  tree_stat_summary = c(mean(nt_lb.df$nn), mean(depth_lb.df$nn), mean(acc_lb.df$nn==TRUE))
  
  names(tree_stat_summary) <- c("nterm", "depth", "acc")
  
  xtable(t(as.matrix(c(tree_stat_summary, pred_cond_summary))), digits = 4)
  
  conv_diag_sum <- cbind(conv_diag(depth_lb.df,n.born.out.par,n.thin), conv_diag(nt_lb.df, n.born.out.par,n.thin), conv_diag(like_df,n.born.out.par,n.thin))
  xtable(conv_diag_sum)
  
  
  if(F){
    copula_uu_clayton_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 1, par = pred_cond_mod$theta_mean)
    
    pred_cond_mod$U1 = copula_uu_clayton_pred[,1]
    pred_cond_mod$U2 = copula_uu_clayton_pred[,2]
    
    # Plot 3D histogram
    
    hist_true <- hist2d(get(paste0("copula_uu_clayton_",i))[,1], get(paste0("copula_uu_clayton_",i))[,2], nbins = c(10,10), show = FALSE)
    
    hist_pred <- hist2d(pred_cond_mod$U1, pred_cond_mod$U2, nbins = c(10,10), show = FALSE)
    
    par(mar = c(1,1,1,1), mfrow = c(1,2))
    
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
      x = hist_pred$x,
      y = hist_pred$y,
      z = hist_pred$counts,
      colvar = NULL,  # disables color mapping
      col = "lightblue",  # solid monochrome color
      border = "black",
      theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
      lighting = TRUE,
      ltheta = 120, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "",
      main = unique(pred_cond_mod$panel.name)
    )
  }
}
