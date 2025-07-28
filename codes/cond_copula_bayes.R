# packages
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
library(calculus)

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
  tau_true_tree <- rep(0,nrow(X_obs))
  tau_true_tree[X_obs<0.33] <- 0.1
  tau_true_tree[(X_obs>=0.33)&(X_obs<0.66)] <- 0.7
  tau_true_tree[(X_obs>=0.66)] <- 0.3
  
  tau_true_1 <- tau_true_tree #+ rnorm(length(tau_true_tree), sd = 0.01)
  tau_true_1 <- matrix(tau_true_1, ncol = 1)
  
  # monotone
  tau_true_2 <- 0.3 + 0.2 * sin(3*X_obs) + 0.3*X_obs^2 #+ rnorm(length(tau_true_1), sd = 0.01)
  # convex
  tau_true_3 <- 0.5 + 0.3 * sin(3*X_obs) #+ rnorm(length(tau_true_1), sd = 0.01)
  # non-convex
  tau_true_4 <- 0.6 - 0.3 * sin(2*X_obs) + 0.2 * sin(4*X_obs) + 0.3 * X_obs^2 #+ rnorm(length(tau_true_1), sd = 0.01)
  
  plot(X_obs, tau_true_1, xlab = "Observations", ylab = "tau")
  plot(X_obs, tau_true_2, xlab = "Observations", ylab = "tau")
  plot(X_obs, tau_true_3, xlab = "Observations", ylab = "tau")
  plot(X_obs, tau_true_4, xlab = "Observations", ylab = "tau")
  
  # gauss # sin(tau*pi/2)
  
  for (i in 1:4) {
    assign(paste0("copula_uu_gauss_",i), BiCopSim(n, family = 1, par = sin(get(paste0("tau_true_",i)) * pi/2)))
  }
  
  # t # sin(tau*pi/2)
  
  for (i in 1:4) {
    assign(paste0("copula_uu_t_",i), BiCopSim(n, family = 2, par = sin(get(paste0("tau_true_",i)) * pi/2), par2 = 3))
  }
  
  # gumbel # 1/(1-tau)
  
  for (i in 1:4) {
    assign(paste0("copula_uu_gumbel_",i), BiCopSim(n, family = 4, par = 1/(1-get(paste0("tau_true_",i)))))
  }
  
  # clayton # (2*tau)/(1-tau)
  
  for (i in 1:4) {
    assign(paste0("copula_uu_clayton_",i), BiCopSim(n, family = 3, par = (2*get(paste0("tau_true_",i)))/(1-get(paste0("tau_true_",i)))))
  }
  
  # frank # numerical
  
  for (i in 1:4) {
    assign(paste0("copula_uu_frank_",i), BiCopSim(n, family = 5, par = param_frank(get(paste0("tau_true_",i)))))
  }
  
  # dataset for prediction
  
  X_obs_pred = matrix(runif(n), ncol = 1)
  
  # normalise predictors 
  X_obs_pred.norm <- as.data.frame(apply(X_obs_pred, 2, \(x) (x - min(x))/(max(x) - min(x))))
  X_obs_pred.norm <- as.matrix(X_obs_pred.norm)
  rownames(X_obs_pred.norm) <- 1:nrow(X_obs_pred)
  
  # tau with tree structure
  tau_true_tree_pred <- rep(0,nrow(X_obs))
  tau_true_tree_pred[X_obs_pred<0.33] <- 0.1
  tau_true_tree_pred[(X_obs_pred>=0.33)&(X_obs_pred<0.66)] <- 0.7
  tau_true_tree_pred[(X_obs_pred>=0.66)] <- 0.3
  
  tau_true_pred_1 <- tau_true_tree_pred #+ rnorm(length(tau_true_tree_pred), sd = 0.01)
  tau_true_pred_1 <- matrix(tau_true_pred_1, ncol = 1)
  
  # monotone
  tau_true_pred_2 <- 0.3 + 0.2 * sin(3*X_obs_pred) + 0.3*X_obs_pred^2 #+ rnorm(length(tau_true_1), sd = 0.01)
  # convex
  tau_true_pred_3 <- 0.5 + 0.3 * sin(3*X_obs_pred) #+ rnorm(length(tau_true_1), sd = 0.01)
  # non-convex
  tau_true_pred_4 <- 0.6 - 0.3 * sin(2*X_obs_pred) + 0.2 * sin(4*X_obs_pred) + 0.3 * X_obs_pred^2 #+ rnorm(length(tau_true_1), sd = 0.01)
  
  # mcmc params
  n.chain_par <- 1
  n.iter_par <- 6000
  n.born.out.par <- 1000
  n.thin <- 1
  incl.split_par <- TRUE
  cont.unif_par <- TRUE
  moves.prob_par <- c(0.4, 0.4, 0.1, 0.1)
}

################################################################################
# source files
################################################################################

# source('mclapply.R') # if run on windows uncomment it

# source('MCMC_BART_copula.R')
source('import_functions.R')
source('test_MCMC_copula_mult_tree.R')


lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) 
################################################################################
# gaussian
################################################################################
if(F){
  n.tree <- 1
  
  for (i in 1) {
    assign(paste0("gauss_mcmc_",i,"_tree_",n.tree), MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                  n.tree = n.tree,
                                                                  X = X_obs.norm,
                                                                  U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                                  U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                                  prior_list = lb.prior.def, 
                                                                  moves.prob = moves.prob_par, 
                                                                  starting.tree = NULL,
                                                                  cont.unif = cont.unif_par,
                                                                  include.split = incl.split_par,
                                                                  prop_mu = 0, prop_sigma = .2,
                                                                  theta_param_1 = 0, theta_param_2 = 1,
                                                                  prior_type = "N",
                                                                  cop_type = "gauss"))
    
    cat('done case', i, '\n')
  }
  
}

save("gauss_mcmc_3_tree_10", file = "gauss_mcmc_3_tree_10.Rdata")
# rm(gauss_mcmc_1_tree_10,gauss_mcmc_2_tree_10,gauss_mcmc_3_tree_10,gauss_mcmc_4_tree_10)

# results

if(F){
  test_case = 1
  
  # load(paste0("gauss_mcmc_",test_case,"_tree_5.Rdata"))
  
  model <- get(paste0("gauss_mcmc_",test_case,"_tree_1"))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs_pred.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  n.thin <- 1
  
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
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
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
if(F){
  n.tree <- 1
  
  for (i in 3) {
    assign(paste0("t_mcmc_",i,"_tree_1"), MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                              n.tree = n.tree,
                                                              X = X_obs.norm,
                                                              U1 = get(paste0("copula_uu_t_",i))[,1],
                                                              U2 = get(paste0("copula_uu_t_",i))[,2],
                                                          prior_list = lb.prior.def, 
                                                          moves.prob = moves.prob_par, 
                                                          starting.tree = NULL,
                                                          cont.unif = cont.unif_par,
                                                          include.split = incl.split_par,
                                                          prop_mu = 0, prop_sigma = .2,
                                                          theta_param_1 = 0, theta_param_2 = 1,
                                                          prior_type = "N",
                                                          cop_type = "t"))
    
    cat('done case', i, '\n')
  }
  
}

save("t_mcmc_4_tree_1", file = "t_mcmc_4_tree_1.Rdata")
# rm(t_mcmc_1_tree_1,t_mcmc_2_tree_1,t_mcmc_3_tree_1,t_mcmc_4_tree_1)

# results

if(F){
  test_case = 3
  
  load(paste0("t_mcmc_",test_case,"_tree_1.Rdata"))
  
  model <- get(paste0("t_mcmc_",test_case,"_tree_1"))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs_pred.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  n.thin <- 1
  
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
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
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
if(F){
  n.tree <- 15
  
  for (i in 3) {
    assign(paste0("gumbel_mcmc_",i,"_tree_5"), MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                   n.tree = n.tree, #n.chain = 10, n.cores = 10,
                                                                   X = X_obs.norm,
                                                                   U1 = get(paste0("copula_uu_gumbel_",i))[,1],
                                                                   U2 = get(paste0("copula_uu_gumbel_",i))[,2],
                                                               prior_list = lb.prior.def, 
                                                               moves.prob = moves.prob_par, 
                                                               starting.tree = NULL,
                                                               cont.unif = cont.unif_par,
                                                               include.split = incl.split_par,
                                                               prop_mu = 0, prop_sigma = rep(.1/n.tree,n.tree),
                                                               theta_param_1 = 0, theta_param_2 = 1,
                                                               var_param_1 = 1, var_param_2 = 2,
                                                               prior_type = "N",
                                                               cop_type = "gumbel"))
    
    cat('done case', i, '\n')
  }
  
}

save("gumbel_mcmc_3_tree_5", file = "gumbel_mcmc_3_tree_15.Rdata")
# rm(gumbel_mcmc_1_tree_1,gumbel_mcmc_2_tree_1,gumbel_mcmc_3_tree_1,gumbel_mcmc_4_tree_1)

# results

if(F){
  test_case = 3
  
  load(paste0("gumbel_mcmc_",test_case,"_tree_5.Rdata"))
  
  model <- get(paste0("gumbel_mcmc_",test_case,"_tree_5"))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs_pred.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  n.thin <- 1
  
  n.born.out.par <- 2500
  
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
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
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
  
  conv_diag(like_df,2500,50)
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
# frank
################################################################################
if(F){
  n.tree <- 1
  
  for (i in 1) {
    assign(paste0("frank_mcmc_",i,"_tree_5"), MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                          n.tree = n.tree, #n.chain = 10, n.cores = 10,
                                                          X = X_obs.norm,
                                                          U1 = get(paste0("copula_uu_frank_",i))[,1],
                                                          U2 = get(paste0("copula_uu_frank_",i))[,2],
                                                          prior_list = lb.prior.def, 
                                                          moves.prob = moves.prob_par, 
                                                          starting.tree = NULL,
                                                          cont.unif = cont.unif_par,
                                                          include.split = incl.split_par,
                                                          prop_mu = 0, prop_sigma = rep(.1/n.tree,n.tree),
                                                          theta_param_1 = 0, theta_param_2 = 1,
                                                          var_param_1 = 2, var_param_2 = 2,
                                                          prior_type = "N",
                                                          cop_type = "frank"))
    
    cat('done case', i, '\n')
  }
  
}

save("frank_mcmc_3_tree_5", file = "frank_mcmc_3_tree_10.Rdata")
# rm(frank_mcmc_1_tree_1,frank_mcmc_2_tree_1,frank_mcmc_3_tree_1,frank_mcmc_4_tree_1)

# results

if(F){
  test_case = 1
  
  load(paste0("frank_mcmc_",test_case,"_tree_5.Rdata"))
  
  model <- get(paste0("frank_mcmc_",test_case,"_tree_5"))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs_pred.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  n.thin <- 1
  
  n.born.out.par <- 2500
  
  pred_val_vec = as.vector(pred_val[(1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))],])
  
  pred_obs = rep(X_obs_pred.norm, each = (n.chain_par * (n.iter_par - n.born.out.par)))
  
  theta_true = rep(param_frank(get(paste0("tau_true_pred_",test_case))), each = (n.chain_par * (n.iter_par - n.born.out.par)))
  
  pred_cond <- data.frame("obs" = pred_obs)
  pred_cond$obs = pred_obs
  pred_cond$theta_true = theta_true
  pred_cond$y = link_frank(pred_val_vec)
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,(n.thin-1)), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  ggplot(pred_cond_mod) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
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
  
  like_val <- apply(pred_val, 1, function(x)loglik_frank(link_frank(x), get(paste0("copula_uu_frank_",i))[,1], get(paste0("copula_uu_frank_",i))[,2]))
  
  like_df <-data.frame("nn" = like_val)
  like_df$idx <- 1:(n.chain_par*n.iter_par)
  
  pl_like <- ggplot(like_df, aes(idx, nn)) + 
    geom_line() + 
    ylab('log-likelihood') +
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  pl_like
  
  conv_diag(like_df,2500,50)
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
    copula_uu_frank_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 5, par = pred_cond_mod$theta_mean)
    
    pred_cond_mod$U1 = copula_uu_frank_pred[,1]
    pred_cond_mod$U2 = copula_uu_frank_pred[,2]
    
    # Plot 3D histogram
    
    hist_true <- hist2d(get(paste0("copula_uu_frank_",i))[,1], get(paste0("copula_uu_frank_",i))[,2], nbins = c(10,10), show = FALSE)
    
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
if(F){
  n.tree <- 1
  
  for (i in 1) {
    assign(paste0("clayton_mcmc_",i,"_tree_1"), MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                    n.tree = n.tree, #n.cores = 10,
                                                                    X = X_obs.norm,
                                                                    U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                    U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                prior_list = lb.prior.def, 
                                                                moves.prob = moves.prob_par, 
                                                                starting.tree = NULL,
                                                                cont.unif = cont.unif_par,
                                                                include.split = incl.split_par,
                                                                prop_mu = 0, prop_sigma = .2,
                                                                theta_param_1 = 0, theta_param_2 = 1,
                                                                prior_type = "N",
                                                                cop_type = "clayton"))
    
    cat('done case', i, '\n')
  }
  
  # all prop_mu = 0, prop_sigma = .2,
  # 3 prop_mu = 0, prop_sigma = 3,
  
}

save("clayton_mcmc_4_tree_1", file = "clayton_mcmc_4_tree_1.Rdata")
# rm(clayton_mcmc_1_tree_1,clayton_mcmc_2_tree_1,clayton_mcmc_3_tree_1,clayton_mcmc_4_tree_1)

# results

if(F){
  test_case = 2
  
  load(paste0("clayton_mcmc_",test_case,"_tree_1.Rdata"))
  
  model <- get(paste0("clayton_mcmc_",test_case,"_tree_1"))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs_pred.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  n.thin <- 1
  
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
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
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
  
  like_df_thin <- like_df[c((1:(n.chain_par * n.iter_par))[rep((n.born.out.par+1):n.iter_par, n.chain_par) + rep(n.iter_par * (0:(n.chain_par-1)), each = (n.iter_par - n.born.out.par))]),]
  
  like_df_thin = na.omit(like_df_thin[c(rep(NA,(n.thin-1)), TRUE),])
  
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
