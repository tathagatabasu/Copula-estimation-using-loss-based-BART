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
library(patchwork)

################################################################################
# data generation
################################################################################
set.seed(1e3)

# load("analysis_sim_dat_red.RData")

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
  
  tau_true_1 <- tau_true_tree + rnorm(length(tau_true_tree), sd = 0.01)
  tau_true_1 <- matrix(tau_true_1, ncol = 1)
  
  # periodoic
  tau_true_2 <- 0.4*sin(2*pi*X_obs) + 0.5 + rnorm(length(tau_true_1), sd = 0.01)
  
  plot(X_obs, tau_true_1, xlab = "Observations", ylab = "tau")
  plot(X_obs, tau_true_2, xlab = "Observations", ylab = "tau")
  
  # gauss # sin(tau*pi/2)
  
  for (i in 1:2) {
    assign(paste0("copula_uu_gauss_",i), BiCopSim(n, family = 1, par = sin(get(paste0("tau_true_",i)) * pi/2)))
  }
  
  # t # sin(tau*pi/2)
  
  for (i in 1:2) {
    assign(paste0("copula_uu_t_",i), BiCopSim(n, family = 2, par = sin(get(paste0("tau_true_",i)) * pi/2), par2 = 3))
  }
  
  # gumbel # 1/(1-tau)
  
  for (i in 1:2) {
    assign(paste0("copula_uu_gumbel_",i), BiCopSim(n, family = 4, par = 1/(1-get(paste0("tau_true_",i)))))
  }
  
  # clayton # (2*tau)/(1-tau)
  
  for (i in 1:2) {
    assign(paste0("copula_uu_clayton_",i), BiCopSim(n, family = 3, par = (2*get(paste0("tau_true_",i)))/(1-get(paste0("tau_true_",i)))))
  }
  
  # frank # numerical
  
  for (i in 1:2) {
    assign(paste0("copula_uu_frank_",i), BiCopSim(n, family = 5, par = BiCopTau2Par(5, get(paste0("tau_true_",i)))))
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
  
  tau_true_pred_1 <- tau_true_tree_pred + rnorm(length(tau_true_tree_pred), sd = 0.01)
  tau_true_pred_1 <- matrix(tau_true_pred_1, ncol = 1)
  
  # periodic
  tau_true_pred_2 <- 0.4*sin(2*pi*X_obs_pred) + 0.5 + rnorm(length(tau_true_1), sd = 0.01)
  
  # mcmc params
  n.chain_par <- 20
  n.iter_par <- 6000
  n.born.out.par <- 1000
  n.thin <- 1
  incl.split_par <- TRUE
  cont.unif_par <- TRUE
  moves.prob_par <- c(0.1, 0.4, 0.25, 0.25)
}

################################################################################
# source files
################################################################################

if(.Platform$OS.type == "windows"){
  source('mclapply.R')
}

# source('MCMC_BART_copula.R')
source('import_functions.R')
source('test_MCMC_copula_mult_tree.R')


lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) 
################################################################################
# gaussian
################################################################################
if(T){
  n.tree <- 10
  
  for (i in 2) {
    assign(paste0("gauss_mcmc_",i,"_tree_",n.tree,"_adapt"), multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                           n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
                                                                           X = X_obs.norm,
                                                                           U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                                           U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                                           prior_list = lb.prior.def, 
                                                                           moves.prob = moves.prob_par, 
                                                                           starting.tree = NULL,
                                                                           cont.unif = cont.unif_par,
                                                                           include.split = incl.split_par,
                                                                           prop_mu = 0, prop_sigma = .4,
                                                                           theta_param_1 = 0, theta_param_2 = 1,
                                                                           var_param_1 = 1, var_param_2 = 2,
                                                                           prior_type = "N",
                                                                           cop_type = "gauss",
                                                                           adapt = T))
    
    cat('done case', i, '\n')
    
    save(list = paste0("gauss_mcmc_",i,"_tree_",n.tree,"_adapt"), file = paste0("gauss_mcmc_",i,"_tree_",n.tree,"_adapt", ".Rdata"))
    rm(list = paste0("gauss_mcmc_",i,"_tree_",n.tree,"_adapt"))
  }
  
  for (i in 2) {
    assign(paste0("gauss_mcmc_",i,"_tree_",n.tree), multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                           n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
                                                                           X = X_obs.norm,
                                                                           U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                                           U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                                           prior_list = lb.prior.def, 
                                                                           moves.prob = moves.prob_par, 
                                                                           starting.tree = NULL,
                                                                           cont.unif = cont.unif_par,
                                                                           include.split = incl.split_par,
                                                                           prop_mu = 0, prop_sigma = .4,
                                                                           theta_param_1 = 0, theta_param_2 = 1,
                                                                           var_param_1 = 1, var_param_2 = 2,
                                                                           prior_type = "N",
                                                                           cop_type = "gauss",
                                                                           adapt = F))
    
    cat('done case', i, '\n')
    
    save(list = paste0("gauss_mcmc_",i,"_tree_",n.tree), file = paste0("gauss_mcmc_",i,"_tree_",n.tree, ".Rdata"))
    rm(list = paste0("gauss_mcmc_",i,"_tree_",n.tree))
  }
  
}

# results

if(F){
  test_case = 2
  
  # model <- get(load(paste0("gauss_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata")))
  # rm(list = paste0("gauss_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  
  model <- get(load(paste0("gauss_mcmc_",test_case,"_tree_",n.tree, ".Rdata")))
  rm(list = paste0("gauss_mcmc_",test_case,"_tree_",n.tree))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  rm(list_pred_lb)
  
  pred_cond <- data.frame("obs" = rep(X_obs.norm, each = (n.chain_par * n.iter_par)))
  pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))), each = (n.chain_par * (n.iter_par)))
  pred_cond$y = BiCopPar2Tau(1,link_gauss(as.vector(pred_val)))
  pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_mod = pred_cond %>%
    filter(idx > n.born.out.par) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df <-data.frame("nn" = apply(pred_val, 1, function(x)loglik_gauss(link_gauss(x), get(paste0("copula_uu_gauss_",test_case))[,1], get(paste0("copula_uu_gauss_",test_case))[,2])))
  like_df$idx <- rep(1:n.iter_par, n.chain_par)
  like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df %>%
    # filter(idx > n.born.out.par) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  library(patchwork)
  pl_like +pl_pred
  
  # nterm
  
  # nt_lb.df <- nterm_BART(model)
  # nt_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # nt_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # 
  # pl_nl <- nt_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "nterm"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_nl
  # 
  # # depth
  # 
  # depth_lb.df <- depth_BART(model)
  # depth_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # depth_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # pl_dp <- depth_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "depth"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_dp
  # 
  # pl_nl + pl_dp
  # 
  # # acceptance
  # 
  # acc_lb.df <- acc_BART(model)
  # acc_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # acc_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.born.out.par]), mean(depth_lb.df$nn[depth_lb.df$idx > n.born.out.par]), mean(acc_lb.df$nn[acc_lb.df$idx > n.born.out.par]==TRUE))
  # 
  # names(tree_stat_summary) <- c("nterm", "depth", "acc")
  # 
  # xtable(t(as.matrix(c(tree_stat_summary, pred_cond_summary))), digits = 4)
  
  # conv_diag_sum <- cbind(conv_diag(depth_lb.df,n.born.out.par,10), conv_diag(nt_lb.df, n.born.out.par,10), conv_diag(like_df,n.born.out.par,10))
  # xtable(conv_diag_sum)
  
  # stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
  # plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  # save(stat_list_adapt, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  # save(plot_list_adapt, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
  
  stat_list <- list("pred" = pred_cond_summary)#, "tree" = tree_stat_summary)
  plot_list <- list("like" = pl_like, "pred" = pl_pred)#, "depth" = pl_dp, "nterm" = pl_nl)
  save(stat_list, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  save(plot_list, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

if(F){
  test_case = 2
  
  model <- get(load(paste0("gauss_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata")))
  rm(list = paste0("gauss_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  
  # model <- get(load(paste0("gauss_mcmc_",test_case,"_tree_",n.tree, ".Rdata")))
  # rm(list = paste0("gauss_mcmc_",test_case,"_tree_",n.tree))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  rm(list_pred_lb)
  
  pred_cond <- data.frame("obs" = rep(X_obs.norm, each = (n.chain_par * n.iter_par)))
  pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))), each = (n.chain_par * (n.iter_par)))
  pred_cond$y = BiCopPar2Tau(1,link_gauss(as.vector(pred_val)))
  pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_mod = pred_cond %>%
    filter(idx > n.born.out.par) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df <-data.frame("nn" = apply(pred_val, 1, function(x)loglik_gauss(link_gauss(x), get(paste0("copula_uu_gauss_",test_case))[,1], get(paste0("copula_uu_gauss_",test_case))[,2])))
  like_df$idx <- rep(1:n.iter_par, n.chain_par)
  like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df %>%
    # filter(idx > n.born.out.par) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  library(patchwork)
  pl_like +pl_pred
  
  # nterm
  
  # nt_lb.df <- nterm_BART(model)
  # nt_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # nt_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # 
  # pl_nl <- nt_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "nterm"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_nl
  # 
  # # depth
  # 
  # depth_lb.df <- depth_BART(model)
  # depth_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # depth_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # pl_dp <- depth_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "depth"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_dp
  # 
  # pl_nl + pl_dp
  # 
  # # acceptance
  # 
  # acc_lb.df <- acc_BART(model)
  # acc_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # acc_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.born.out.par]), mean(depth_lb.df$nn[depth_lb.df$idx > n.born.out.par]), mean(acc_lb.df$nn[acc_lb.df$idx > n.born.out.par]==TRUE))
  # 
  # names(tree_stat_summary) <- c("nterm", "depth", "acc")
  # 
  # xtable(t(as.matrix(c(tree_stat_summary, pred_cond_summary))), digits = 4)
  
  # conv_diag_sum <- cbind(conv_diag(depth_lb.df,n.born.out.par,10), conv_diag(nt_lb.df, n.born.out.par,10), conv_diag(like_df,n.born.out.par,10))
  # xtable(conv_diag_sum)
  
  stat_list_adapt <- list("pred" = pred_cond_summary)#, "tree" = tree_stat_summary)
  plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)#, "depth" = pl_dp, "nterm" = pl_nl)
  save(stat_list_adapt, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  save(plot_list_adapt, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
  
  # stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
  # plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  # save(stat_list, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  # save(plot_list, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

################################################################################
# t
################################################################################
if(T){
  n.tree <- 10
  
  for (i in 2) {
    assign(paste0("t_mcmc_",i,"_tree_",n.tree,"_adapt"), multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                       n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
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
                                                                       var_param_1 = 1, var_param_2 = 2,
                                                                       prior_type = "N",
                                                                       cop_type = "t",
                                                                       adapt = T))
    
    cat('done case', i, '\n')
    
    save(list = paste0("t_mcmc_",i,"_tree_",n.tree,"_adapt"), file = paste0("t_mcmc_",i,"_tree_",n.tree,"_adapt", ".Rdata"))
    rm(list = paste0("t_mcmc_",i,"_tree_",n.tree,"_adapt"))
  }
  
  for (i in 2) {
    assign(paste0("t_mcmc_",i,"_tree_",n.tree), multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                       n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
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
                                                                       var_param_1 = 1, var_param_2 = 2,
                                                                       prior_type = "N",
                                                                       cop_type = "t",
                                                                       adapt = F))
    
    cat('done case', i, '\n')
    
    save(list = paste0("t_mcmc_",i,"_tree_",n.tree), file = paste0("t_mcmc_",i,"_tree_",n.tree, ".Rdata"))
    rm(list = paste0("t_mcmc_",i,"_tree_",n.tree))
  }
  
}

# results

if(F){
  test_case = 2
  
  # model <- get(load(paste0("t_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata")))
  # rm(list = paste0("t_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  
  model <- get(load(paste0("t_mcmc_",test_case,"_tree_",n.tree, ".Rdata")))
  rm(list = paste0("t_mcmc_",test_case,"_tree_",n.tree))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  rm(list_pred_lb)
  
  pred_cond <- data.frame("obs" = rep(X_obs.norm, each = (n.chain_par * n.iter_par)))
  pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))), each = (n.chain_par * (n.iter_par)))
  pred_cond$y = BiCopPar2Tau(2,link_t(as.vector(pred_val)))
  pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_mod = pred_cond %>%
    filter(idx > n.born.out.par) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df <-data.frame("nn" = apply(pred_val, 1, function(x)loglik_t(link_t(x), get(paste0("copula_uu_t_",test_case))[,1], get(paste0("copula_uu_t_",test_case))[,2])))
  like_df$idx <- rep(1:n.iter_par, n.chain_par)
  like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df %>%
    # filter(idx > n.born.out.par) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  library(patchwork)
  pl_like +pl_pred
  
  # nterm
  
  # nt_lb.df <- nterm_BART(model)
  # nt_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # nt_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # 
  # pl_nl <- nt_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "nterm"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_nl
  # 
  # # depth
  # 
  # depth_lb.df <- depth_BART(model)
  # depth_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # depth_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # pl_dp <- depth_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "depth"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_dp
  # 
  # pl_nl + pl_dp
  # 
  # # acceptance
  # 
  # acc_lb.df <- acc_BART(model)
  # acc_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # acc_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.born.out.par]), mean(depth_lb.df$nn[depth_lb.df$idx > n.born.out.par]), mean(acc_lb.df$nn[acc_lb.df$idx > n.born.out.par]==TRUE))
  # 
  # names(tree_stat_summary) <- c("nterm", "depth", "acc")
  # 
  # xtable(t(as.matrix(c(tree_stat_summary, pred_cond_summary))), digits = 4)
  
  # conv_diag_sum <- cbind(conv_diag(depth_lb.df,n.born.out.par,10), conv_diag(nt_lb.df, n.born.out.par,10), conv_diag(like_df,n.born.out.par,10))
  # xtable(conv_diag_sum)
  
  # stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
  # plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  # save(stat_list_adapt, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  # save(plot_list_adapt, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
  
  stat_list <- list("pred" = pred_cond_summary)#, "tree" = tree_stat_summary)
  plot_list <- list("like" = pl_like, "pred" = pl_pred)#, "depth" = pl_dp, "nterm" = pl_nl)
  save(stat_list, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  save(plot_list, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

if(F){
  test_case = 2
  
  model <- get(load(paste0("t_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata")))
  rm(list = paste0("t_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  
  # model <- get(load(paste0("t_mcmc_",test_case,"_tree_",n.tree, ".Rdata")))
  # rm(list = paste0("t_mcmc_",test_case,"_tree_",n.tree))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  rm(list_pred_lb)
  
  pred_cond <- data.frame("obs" = rep(X_obs.norm, each = (n.chain_par * n.iter_par)))
  pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))), each = (n.chain_par * (n.iter_par)))
  pred_cond$y = BiCopPar2Tau(2,link_t(as.vector(pred_val)))
  pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_mod = pred_cond %>%
    filter(idx > n.born.out.par) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df <-data.frame("nn" = apply(pred_val, 1, function(x)loglik_t(link_t(x), get(paste0("copula_uu_t_",test_case))[,1], get(paste0("copula_uu_t_",test_case))[,2])))
  like_df$idx <- rep(1:n.iter_par, n.chain_par)
  like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df %>%
    # filter(idx > n.born.out.par) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  library(patchwork)
  pl_like +pl_pred
  
  # nterm
  
  # nt_lb.df <- nterm_BART(model)
  # nt_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # nt_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # 
  # pl_nl <- nt_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "nterm"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_nl
  # 
  # # depth
  # 
  # depth_lb.df <- depth_BART(model)
  # depth_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # depth_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # pl_dp <- depth_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "depth"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_dp
  # 
  # pl_nl + pl_dp
  # 
  # # acceptance
  # 
  # acc_lb.df <- acc_BART(model)
  # acc_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # acc_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.born.out.par]), mean(depth_lb.df$nn[depth_lb.df$idx > n.born.out.par]), mean(acc_lb.df$nn[acc_lb.df$idx > n.born.out.par]==TRUE))
  # 
  # names(tree_stat_summary) <- c("nterm", "depth", "acc")
  # 
  # xtable(t(as.matrix(c(tree_stat_summary, pred_cond_summary))), digits = 4)
  
  # conv_diag_sum <- cbind(conv_diag(depth_lb.df,n.born.out.par,10), conv_diag(nt_lb.df, n.born.out.par,10), conv_diag(like_df,n.born.out.par,10))
  # xtable(conv_diag_sum)
  
  stat_list_adapt <- list("pred" = pred_cond_summary)#, "tree" = tree_stat_summary)
  plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)#, "depth" = pl_dp, "nterm" = pl_nl)
  save(stat_list_adapt, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  save(plot_list_adapt, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
  
  # stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
  # plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  # save(stat_list, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  # save(plot_list, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

################################################################################
# Gumbel
################################################################################
if(T){
  n.tree <- 10
  
  for (i in 2) {
    assign(paste0("gumbel_mcmc_",i,"_tree_",n.tree,"_adapt"), multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                            n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
                                                                            X = X_obs.norm,
                                                                            U1 = get(paste0("copula_uu_gumbel_",i))[,1],
                                                                            U2 = get(paste0("copula_uu_gumbel_",i))[,2],
                                                                            prior_list = lb.prior.def, 
                                                                            moves.prob = moves.prob_par, 
                                                                            starting.tree = NULL,
                                                                            cont.unif = cont.unif_par,
                                                                            include.split = incl.split_par,
                                                                            prop_mu = 0, prop_sigma = .2/n.tree,
                                                                            theta_param_1 = 0, theta_param_2 = 1,
                                                                            var_param_1 = 1, var_param_2 = 2,
                                                                            prior_type = "N",
                                                                            cop_type = "gumbel",
                                                                            adapt = T))
    
    cat('done case', i, '\n')
    
    save(list = paste0("gumbel_mcmc_",i,"_tree_",n.tree,"_adapt"), file = paste0("gumbel_mcmc_",i,"_tree_",n.tree,"_adapt", ".Rdata"))
    rm(list = paste0("gumbel_mcmc_",i,"_tree_",n.tree,"_adapt"))
  }
  
  for (i in 2) {
    assign(paste0("gumbel_mcmc_",i,"_tree_",n.tree), multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                            n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
                                                                            X = X_obs.norm,
                                                                            U1 = get(paste0("copula_uu_gumbel_",i))[,1],
                                                                            U2 = get(paste0("copula_uu_gumbel_",i))[,2],
                                                                            prior_list = lb.prior.def, 
                                                                            moves.prob = moves.prob_par, 
                                                                            starting.tree = NULL,
                                                                            cont.unif = cont.unif_par,
                                                                            include.split = incl.split_par,
                                                                            prop_mu = 0, prop_sigma = .2/n.tree,
                                                                            theta_param_1 = 0, theta_param_2 = 1,
                                                                            var_param_1 = 1, var_param_2 = 2,
                                                                            prior_type = "N",
                                                                            cop_type = "gumbel",
                                                                            adapt = F))
    
    cat('done case', i, '\n')
    
    save(list = paste0("gumbel_mcmc_",i,"_tree_",n.tree), file = paste0("gumbel_mcmc_",i,"_tree_",n.tree, ".Rdata"))
    rm(list = paste0("gumbel_mcmc_",i,"_tree_",n.tree))
  }
  
}

# results

if(F){
  test_case = 2
  
  # model <- get(load(paste0("gumbel_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata")))
  # rm(list = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  
  model <- get(load(paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, ".Rdata")))
  rm(list = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  rm(list_pred_lb)
  
  pred_cond <- data.frame("obs" = rep(X_obs.norm, each = (n.chain_par * n.iter_par)))
  pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))), each = (n.chain_par * (n.iter_par)))
  pred_cond$y = BiCopPar2Tau(4,link_gumbel(as.vector(pred_val)))
  pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_mod = pred_cond %>%
    filter(idx > n.born.out.par) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df <-data.frame("nn" = apply(pred_val, 1, function(x)loglik_gumbel(link_gumbel(x), get(paste0("copula_uu_gumbel_",test_case))[,1], get(paste0("copula_uu_gumbel_",test_case))[,2])))
  like_df$idx <- rep(1:n.iter_par, n.chain_par)
  like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df %>%
    # filter(idx > n.born.out.par) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  library(patchwork)
  pl_like +pl_pred
  
  # nterm
  
  # nt_lb.df <- nterm_BART(model)
  # nt_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # nt_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # 
  # pl_nl <- nt_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "nterm"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_nl
  # 
  # # depth
  # 
  # depth_lb.df <- depth_BART(model)
  # depth_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # depth_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # pl_dp <- depth_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "depth"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_dp
  # 
  # pl_nl + pl_dp
  # 
  # # acceptance
  # 
  # acc_lb.df <- acc_BART(model)
  # acc_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # acc_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.born.out.par]), mean(depth_lb.df$nn[depth_lb.df$idx > n.born.out.par]), mean(acc_lb.df$nn[acc_lb.df$idx > n.born.out.par]==TRUE))
  # 
  # names(tree_stat_summary) <- c("nterm", "depth", "acc")
  # 
  # xtable(t(as.matrix(c(tree_stat_summary, pred_cond_summary))), digits = 4)
  
  # conv_diag_sum <- cbind(conv_diag(depth_lb.df,n.born.out.par,10), conv_diag(nt_lb.df, n.born.out.par,10), conv_diag(like_df,n.born.out.par,10))
  # xtable(conv_diag_sum)
  
  # stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
  # plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  # save(stat_list_adapt, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  # save(plot_list_adapt, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
  
  stat_list <- list("pred" = pred_cond_summary)#, "tree" = tree_stat_summary)
  plot_list <- list("like" = pl_like, "pred" = pl_pred)#, "depth" = pl_dp, "nterm" = pl_nl)
  save(stat_list, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  save(plot_list, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

if(F){
  test_case = 2
  
  model <- get(load(paste0("gumbel_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata")))
  rm(list = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  
  # model <- get(load(paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, ".Rdata")))
  # rm(list = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  rm(list_pred_lb)
  
  pred_cond <- data.frame("obs" = rep(X_obs.norm, each = (n.chain_par * n.iter_par)))
  pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))), each = (n.chain_par * (n.iter_par)))
  pred_cond$y = BiCopPar2Tau(4,link_gumbel(as.vector(pred_val)))
  pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_mod = pred_cond %>%
    filter(idx > n.born.out.par) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df <-data.frame("nn" = apply(pred_val, 1, function(x)loglik_gumbel(link_gumbel(x), get(paste0("copula_uu_gumbel_",test_case))[,1], get(paste0("copula_uu_gumbel_",test_case))[,2])))
  like_df$idx <- rep(1:n.iter_par, n.chain_par)
  like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df %>%
    # filter(idx > n.born.out.par) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  library(patchwork)
  pl_like +pl_pred
  
  # nterm
  
  # nt_lb.df <- nterm_BART(model)
  # nt_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # nt_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # 
  # pl_nl <- nt_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "nterm"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_nl
  # 
  # # depth
  # 
  # depth_lb.df <- depth_BART(model)
  # depth_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # depth_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # pl_dp <- depth_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "depth"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_dp
  # 
  # pl_nl + pl_dp
  # 
  # # acceptance
  # 
  # acc_lb.df <- acc_BART(model)
  # acc_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # acc_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.born.out.par]), mean(depth_lb.df$nn[depth_lb.df$idx > n.born.out.par]), mean(acc_lb.df$nn[acc_lb.df$idx > n.born.out.par]==TRUE))
  # 
  # names(tree_stat_summary) <- c("nterm", "depth", "acc")
  # 
  # xtable(t(as.matrix(c(tree_stat_summary, pred_cond_summary))), digits = 4)
  
  # conv_diag_sum <- cbind(conv_diag(depth_lb.df,n.born.out.par,10), conv_diag(nt_lb.df, n.born.out.par,10), conv_diag(like_df,n.born.out.par,10))
  # xtable(conv_diag_sum)
  
  stat_list_adapt <- list("pred" = pred_cond_summary)#, "tree" = tree_stat_summary)
  plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)#, "depth" = pl_dp, "nterm" = pl_nl)
  save(stat_list_adapt, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  save(plot_list_adapt, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
  
  # stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
  # plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  # save(stat_list, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  # save(plot_list, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

################################################################################
# frank
################################################################################
if(T){
  n.tree <- 10
  
  for (i in 2) {
    assign(paste0("frank_mcmc_",i,"_tree_",n.tree,"_adapt"), multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                           n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
                                                                           X = X_obs.norm,
                                                                           U1 = get(paste0("copula_uu_frank_",i))[,1],
                                                                           U2 = get(paste0("copula_uu_frank_",i))[,2],
                                                                           prior_list = lb.prior.def, 
                                                                           moves.prob = moves.prob_par, 
                                                                           starting.tree = NULL,
                                                                           cont.unif = cont.unif_par,
                                                                           include.split = incl.split_par,
                                                                           prop_mu = 0, prop_sigma = rep(1/n.tree,n.tree),
                                                                           theta_param_1 = 0, theta_param_2 = 1,
                                                                           var_param_1 = 2, var_param_2 = 2,
                                                                           prior_type = "N",
                                                                           cop_type = "frank",
                                                                           adapt = T))
    
    cat('done case', i, '\n')
    
    save(list = paste0("frank_mcmc_",i,"_tree_",n.tree,"_adapt"), file = paste0("frank_mcmc_",i,"_tree_",n.tree,"_adapt", ".Rdata"))
    rm(list = paste0("frank_mcmc_",i,"_tree_",n.tree,"_adapt"))
  }
  
  for (i in 2) {
    assign(paste0("frank_mcmc_",i,"_tree_",n.tree), multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                           n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
                                                                           X = X_obs.norm,
                                                                           U1 = get(paste0("copula_uu_frank_",i))[,1],
                                                                           U2 = get(paste0("copula_uu_frank_",i))[,2],
                                                                           prior_list = lb.prior.def, 
                                                                           moves.prob = moves.prob_par, 
                                                                           starting.tree = NULL,
                                                                           cont.unif = cont.unif_par,
                                                                           include.split = incl.split_par,
                                                                           prop_mu = 0, prop_sigma = rep(1/n.tree,n.tree),
                                                                           theta_param_1 = 0, theta_param_2 = 1,
                                                                           var_param_1 = 2, var_param_2 = 2,
                                                                           prior_type = "N",
                                                                           cop_type = "frank",
                                                                           adapt = F))
    
    cat('done case', i, '\n')
    
    save(list = paste0("frank_mcmc_",i,"_tree_",n.tree), file = paste0("frank_mcmc_",i,"_tree_",n.tree, ".Rdata"))
    rm(list = paste0("frank_mcmc_",i,"_tree_",n.tree))
  }
}

# results

if(F){
  test_case = 2
  
  # model <- get(load(paste0("frank_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata")))
  # rm(list = paste0("frank_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  
  model <- get(load(paste0("frank_mcmc_",test_case,"_tree_",n.tree, ".Rdata")))
  rm(list = paste0("frank_mcmc_",test_case,"_tree_",n.tree))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  rm(list_pred_lb)
  
  pred_cond <- data.frame("obs" = rep(X_obs.norm, each = (n.chain_par * n.iter_par)))
  pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))), each = (n.chain_par * (n.iter_par)))
  pred_cond$y = BiCopPar2Tau(5,link_frank(as.vector(pred_val)))
  pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_mod = pred_cond %>%
    filter(idx > n.born.out.par) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df <-data.frame("nn" = apply(pred_val, 1, function(x)loglik_frank(link_frank(x), get(paste0("copula_uu_frank_",test_case))[,1], get(paste0("copula_uu_frank_",test_case))[,2])))
  like_df$idx <- rep(1:n.iter_par, n.chain_par)
  like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df %>%
    # filter(idx > n.born.out.par) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  library(patchwork)
  pl_like +pl_pred
  
  # nterm
  
  # nt_lb.df <- nterm_BART(model)
  # nt_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # nt_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # 
  # pl_nl <- nt_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "nterm"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_nl
  # 
  # # depth
  # 
  # depth_lb.df <- depth_BART(model)
  # depth_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # depth_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # pl_dp <- depth_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "depth"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_dp
  # 
  # pl_nl + pl_dp
  # 
  # # acceptance
  # 
  # acc_lb.df <- acc_BART(model)
  # acc_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # acc_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.born.out.par]), mean(depth_lb.df$nn[depth_lb.df$idx > n.born.out.par]), mean(acc_lb.df$nn[acc_lb.df$idx > n.born.out.par]==TRUE))
  # 
  # names(tree_stat_summary) <- c("nterm", "depth", "acc")
  # 
  # xtable(t(as.matrix(c(tree_stat_summary, pred_cond_summary))), digits = 4)
  
  # conv_diag_sum <- cbind(conv_diag(depth_lb.df,n.born.out.par,10), conv_diag(nt_lb.df, n.born.out.par,10), conv_diag(like_df,n.born.out.par,10))
  # xtable(conv_diag_sum)
  
  # stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
  # plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  # save(stat_list_adapt, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  # save(plot_list_adapt, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
  
  stat_list <- list("pred" = pred_cond_summary)#, "tree" = tree_stat_summary)
  plot_list <- list("like" = pl_like, "pred" = pl_pred)#, "depth" = pl_dp, "nterm" = pl_nl)
  save(stat_list, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  save(plot_list, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

if(F){
  test_case = 2
  
  model <- get(load(paste0("frank_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata")))
  rm(list = paste0("frank_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  
  # model <- get(load(paste0("frank_mcmc_",test_case,"_tree_",n.tree, ".Rdata")))
  # rm(list = paste0("frank_mcmc_",test_case,"_tree_",n.tree))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  rm(list_pred_lb)
  
  pred_cond <- data.frame("obs" = rep(X_obs.norm, each = (n.chain_par * n.iter_par)))
  pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))), each = (n.chain_par * (n.iter_par)))
  pred_cond$y = BiCopPar2Tau(5,link_frank(as.vector(pred_val)))
  pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_mod = pred_cond %>%
    filter(idx > n.born.out.par) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df <-data.frame("nn" = apply(pred_val, 1, function(x)loglik_frank(link_frank(x), get(paste0("copula_uu_frank_",test_case))[,1], get(paste0("copula_uu_frank_",test_case))[,2])))
  like_df$idx <- rep(1:n.iter_par, n.chain_par)
  like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df %>%
    # filter(idx > n.born.out.par) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  library(patchwork)
  pl_like +pl_pred
  
  # nterm
  
  # nt_lb.df <- nterm_BART(model)
  # nt_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # nt_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # 
  # pl_nl <- nt_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "nterm"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_nl
  # 
  # # depth
  # 
  # depth_lb.df <- depth_BART(model)
  # depth_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # depth_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # pl_dp <- depth_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "depth"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_dp
  # 
  # pl_nl + pl_dp
  # 
  # # acceptance
  # 
  # acc_lb.df <- acc_BART(model)
  # acc_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # acc_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.born.out.par]), mean(depth_lb.df$nn[depth_lb.df$idx > n.born.out.par]), mean(acc_lb.df$nn[acc_lb.df$idx > n.born.out.par]==TRUE))
  # 
  # names(tree_stat_summary) <- c("nterm", "depth", "acc")
  # 
  # xtable(t(as.matrix(c(tree_stat_summary, pred_cond_summary))), digits = 4)
  
  # conv_diag_sum <- cbind(conv_diag(depth_lb.df,n.born.out.par,10), conv_diag(nt_lb.df, n.born.out.par,10), conv_diag(like_df,n.born.out.par,10))
  # xtable(conv_diag_sum)
  
  stat_list_adapt <- list("pred" = pred_cond_summary)#, "tree" = tree_stat_summary)
  plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)#, "depth" = pl_dp, "nterm" = pl_nl)
  save(stat_list_adapt, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  save(plot_list_adapt, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
  
  # stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
  # plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  # save(stat_list, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  # save(plot_list, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

################################################################################
# clayton
################################################################################
if(T){
  n.tree <- 10
  
  for (i in 2) {
    assign(paste0("clayton_mcmc_",i,"_tree_",n.tree,"_adapt"), multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                             n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
                                                                             X = X_obs.norm,
                                                                             U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                             U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                             prior_list = lb.prior.def, 
                                                                             moves.prob = moves.prob_par, 
                                                                             starting.tree = NULL,
                                                                             cont.unif = cont.unif_par,
                                                                             include.split = incl.split_par,
                                                                             prop_mu = 0, prop_sigma = rep(1/n.tree,n.tree),
                                                                             theta_param_1 = 0, theta_param_2 = 1,
                                                                             var_param_1 = 2, var_param_2 = 2,
                                                                             prior_type = "N",
                                                                             cop_type = "clayton",
                                                                             adapt = T))
    
    cat('done case', i, '\n')
    
    save(list = paste0("clayton_mcmc_",i,"_tree_",n.tree,"_adapt"), file = paste0("clayton_mcmc_",i,"_tree_",n.tree,"_adapt", ".Rdata"))
    rm(list = paste0("clayton_mcmc_",i,"_tree_",n.tree,"_adapt"))
  }
  
  for (i in 2) {
    assign(paste0("clayton_mcmc_",i,"_tree_",n.tree), multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                             n.tree = n.tree, n.chain = n.chain_par, n.cores = 10,
                                                                             X = X_obs.norm,
                                                                             U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                             U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                             prior_list = lb.prior.def, 
                                                                             moves.prob = moves.prob_par, 
                                                                             starting.tree = NULL,
                                                                             cont.unif = cont.unif_par,
                                                                             include.split = incl.split_par,
                                                                             prop_mu = 0, prop_sigma = rep(1/n.tree,n.tree),
                                                                             theta_param_1 = 0, theta_param_2 = 1,
                                                                             var_param_1 = 2, var_param_2 = 2,
                                                                             prior_type = "N",
                                                                             cop_type = "clayton",
                                                                             adapt = F))
    
    cat('done case', i, '\n')
    
    save(list = paste0("clayton_mcmc_",i,"_tree_",n.tree), file = paste0("clayton_mcmc_",i,"_tree_",n.tree, ".Rdata"))
    rm(list = paste0("clayton_mcmc_",i,"_tree_",n.tree))
  }
  
}
# results

if(F){
  test_case = 2
  
  # model <- get(load(paste0("clayton_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata")))
  # rm(list = paste0("clayton_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  
  model <- get(load(paste0("clayton_mcmc_",test_case,"_tree_",n.tree, ".Rdata")))
  rm(list = paste0("clayton_mcmc_",test_case,"_tree_",n.tree))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  rm(list_pred_lb)
  
  pred_cond <- data.frame("obs" = rep(X_obs.norm, each = (n.chain_par * n.iter_par)))
  pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))), each = (n.chain_par * (n.iter_par)))
  pred_cond$y = BiCopPar2Tau(3,link_clayton(as.vector(pred_val)))
  pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_mod = pred_cond %>%
    filter(idx > n.born.out.par) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df <-data.frame("nn" = apply(pred_val, 1, function(x)loglik_clayton(link_clayton(x), get(paste0("copula_uu_clayton_",test_case))[,1], get(paste0("copula_uu_clayton_",test_case))[,2])))
  like_df$idx <- rep(1:n.iter_par, n.chain_par)
  like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df %>%
    # filter(idx > n.born.out.par) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  library(patchwork)
  pl_like +pl_pred
  
  # nterm
  
  # nt_lb.df <- nterm_BART(model)
  # nt_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # nt_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # 
  # pl_nl <- nt_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "nterm"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_nl
  # 
  # # depth
  # 
  # depth_lb.df <- depth_BART(model)
  # depth_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # depth_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # pl_dp <- depth_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "depth"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_dp
  # 
  # pl_nl + pl_dp
  # 
  # # acceptance
  # 
  # acc_lb.df <- acc_BART(model)
  # acc_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # acc_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.born.out.par]), mean(depth_lb.df$nn[depth_lb.df$idx > n.born.out.par]), mean(acc_lb.df$nn[acc_lb.df$idx > n.born.out.par]==TRUE))
  # 
  # names(tree_stat_summary) <- c("nterm", "depth", "acc")
  # 
  # xtable(t(as.matrix(c(tree_stat_summary, pred_cond_summary))), digits = 4)
  
  # conv_diag_sum <- cbind(conv_diag(depth_lb.df,n.born.out.par,10), conv_diag(nt_lb.df, n.born.out.par,10), conv_diag(like_df,n.born.out.par,10))
  # xtable(conv_diag_sum)
  
  # stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
  # plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  # save(stat_list_adapt, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  # save(plot_list_adapt, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
  
  stat_list <- list("pred" = pred_cond_summary)#, "tree" = tree_stat_summary)
  plot_list <- list("like" = pl_like, "pred" = pl_pred)#, "depth" = pl_dp, "nterm" = pl_nl)
  save(stat_list, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  save(plot_list, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

if(F){
  test_case = 2
  
  model <- get(load(paste0("clayton_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata")))
  rm(list = paste0("clayton_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  
  # model <- get(load(paste0("clayton_mcmc_",test_case,"_tree_",n.tree, ".Rdata")))
  # rm(list = paste0("clayton_mcmc_",test_case,"_tree_",n.tree))
  
  list_pred_lb <- lapply(1:length(model$trees), \(idx) BART_calculate_pred(model$trees[[idx]], X_obs.norm))
  
  pred_val = do.call(rbind,list_pred_lb)
  
  rm(list_pred_lb)
  
  pred_cond <- data.frame("obs" = rep(X_obs.norm, each = (n.chain_par * n.iter_par)))
  pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))), each = (n.chain_par * (n.iter_par)))
  pred_cond$y = BiCopPar2Tau(3,link_clayton(as.vector(pred_val)))
  pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_mod = pred_cond %>%
    filter(idx > n.born.out.par) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df <-data.frame("nn" = apply(pred_val, 1, function(x)loglik_clayton(link_clayton(x), get(paste0("copula_uu_clayton_",test_case))[,1], get(paste0("copula_uu_clayton_",test_case))[,2])))
  like_df$idx <- rep(1:n.iter_par, n.chain_par)
  like_df$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df %>%
    # filter(idx > n.born.out.par) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  library(patchwork)
  pl_like +pl_pred
  
  # nterm
  
  # nt_lb.df <- nterm_BART(model)
  # nt_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # nt_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # 
  # pl_nl <- nt_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "nterm"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_nl
  # 
  # # depth
  # 
  # depth_lb.df <- depth_BART(model)
  # depth_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # depth_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # pl_dp <- depth_lb.df %>%
  #   # filter(idx > n.born.out.par) %>%
  #   ggplot(aes(x = idx, y = nn, color = factor(chain))) +
  #   geom_line() +
  #   labs(
  #     x = "Iteration",
  #     y = "depth"
  #   ) +
  #   guides(color = "none") +
  #   theme_minimal()
  # 
  # pl_dp
  # 
  # pl_nl + pl_dp
  # 
  # # acceptance
  # 
  # acc_lb.df <- acc_BART(model)
  # acc_lb.df$idx <- rep(1:n.iter_par, n.chain_par)
  # acc_lb.df$chain <- rep(1:n.chain_par, each = n.iter_par)
  # 
  # tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.born.out.par]), mean(depth_lb.df$nn[depth_lb.df$idx > n.born.out.par]), mean(acc_lb.df$nn[acc_lb.df$idx > n.born.out.par]==TRUE))
  # 
  # names(tree_stat_summary) <- c("nterm", "depth", "acc")
  # 
  # xtable(t(as.matrix(c(tree_stat_summary, pred_cond_summary))), digits = 4)
  
  # conv_diag_sum <- cbind(conv_diag(depth_lb.df,n.born.out.par,10), conv_diag(nt_lb.df, n.born.out.par,10), conv_diag(like_df,n.born.out.par,10))
  # xtable(conv_diag_sum)
  
  stat_list_adapt <- list("pred" = pred_cond_summary)#, "tree" = tree_stat_summary)
  plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)#, "depth" = pl_dp, "nterm" = pl_nl)
  save(stat_list_adapt, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  save(plot_list_adapt, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
  
  # stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
  # plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  # save(stat_list, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  # save(plot_list, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

# rm(acc_lb.df, depth_lb.df, like_df, like_df_thin, list_pred_lb, model, nt_lb.df, pl_dp, pl_like, pl_nl, pred_cond, pred_cond_mod, pred_cond_stat, pred_cond_thin,
#    pred_obs, pred_val, pred_val_vec, theta_true, pred_cond_summary, tree_stat_summary, conv_diag_sum, pl_pred, plot_list, plot_list_adapt,
#    pred_cond_mod_avg,stat_list,stat_list_adapt)
