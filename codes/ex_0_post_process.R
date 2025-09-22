if(F){
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(MASS)   # For multivariate normal functions
  library(coda)   # For MCMC diagnostics
  library(plot3D)
  library(gplots)
  library(xtable)
  library(patchwork)
  library(parallel)
  library(mc2d)
  library(pracma)
  
}

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



gauss <- get(load("gauss_mcmc_1_tree_1_plot.Rdata"))
gauss_stat <- get(load("gauss_mcmc_1_tree_1_stat.Rdata"))

gauss_adapt <- get(load("gauss_mcmc_1_tree_1_plot_adapt.Rdata"))
gauss_stat_adapt <- get(load("gauss_mcmc_1_tree_1_stat_adapt.Rdata"))

gauss_like_true <- loglik_gauss((sin(tau_true_1 * pi/2)), copula_uu_gauss_1[,1] , copula_uu_gauss_1[,2])
# gauss_like_true <- loglik_gauss((sin(tau_true_2 * pi/2)), copula_uu_gauss_2[,1] , copula_uu_gauss_2[,2])


t <- get(load("t_mcmc_1_tree_1_plot.Rdata"))
t_stat <- get(load("t_mcmc_1_tree_1_stat.Rdata"))

t_adapt <- get(load("t_mcmc_1_tree_1_plot_adapt.Rdata"))
t_stat_adapt <- get(load("t_mcmc_1_tree_1_stat_adapt.Rdata"))

t_like_true <- loglik_t((sin(tau_true_1 * pi/2)), copula_uu_t_1[,1] , copula_uu_t_1[,2])
# t_like_true <- loglik_t((sin(tau_true_2 * pi/2)), copula_uu_t_2[,1] , copula_uu_t_2[,2])


clayton <- get(load("clayton_mcmc_1_tree_1_plot.Rdata"))
clayton_stat <- get(load("clayton_mcmc_1_tree_1_stat.Rdata"))

clayton_adapt <- get(load("clayton_mcmc_1_tree_1_plot_adapt.Rdata"))
clayton_stat_adapt <- get(load("clayton_mcmc_1_tree_1_stat_adapt.Rdata"))

clayton_like_true <- loglik_clayton((2*tau_true_1 / (1-tau_true_1)), copula_uu_clayton_1[,1] , copula_uu_clayton_1[,2])
# clayton_like_true <- loglik_clayton((2*tau_true_2 / (1-tau_true_2)), copula_uu_clayton_2[,1] , copula_uu_clayton_2[,2])


gumbel <- get(load("gumbel_mcmc_1_tree_1_plot.Rdata"))
gumbel_stat <- get(load("gumbel_mcmc_1_tree_1_stat.Rdata"))

gumbel_adapt <- get(load("gumbel_mcmc_1_tree_1_plot_adapt.Rdata"))
gumbel_stat_adapt <- get(load("gumbel_mcmc_1_tree_1_stat_adapt.Rdata"))

gumbel_like_true <- loglik_gumbel((1 / (1-tau_true_1)), copula_uu_gumbel_1[,1] , copula_uu_gumbel_1[,2])
# gumbel_like_true <- loglik_gumbel((1 / (1-tau_true_2)), copula_uu_gumbel_2[,1] , copula_uu_gumbel_2[,2])


frank <- get(load("frank_mcmc_1_tree_1_plot.Rdata"))
frank_stat <- get(load("frank_mcmc_1_tree_1_stat.Rdata"))

frank_adapt <- get(load("frank_mcmc_1_tree_1_plot_adapt.Rdata"))
frank_stat_adapt <- get(load("frank_mcmc_1_tree_1_stat_adapt.Rdata"))


frank_like_true <- loglik_frank(BiCopTau2Par(5,tau_true_1), copula_uu_frank_1[,1] , copula_uu_frank_1[,2])
# frank_like_true <- loglik_frank(BiCopTau2Par(5,tau_true_2), copula_uu_frank_2[,1] , copula_uu_frank_2[,2])

c(gauss_like_true, t_like_true, clayton_like_true, gumbel_like_true, frank_like_true)

summ_stat_acc <- rbind(c(as.vector(gauss_stat$pred),as.vector(gauss_stat_adapt$pred)),
                       c(as.vector(t_stat$pred),as.vector(t_stat_adapt$pred)),
                       c(as.vector(clayton_stat$pred),as.vector(clayton_stat_adapt$pred)),
                       c(as.vector(gumbel_stat$pred),as.vector(gumbel_stat_adapt$pred)),
                       c(as.vector(frank_stat$pred),as.vector(frank_stat_adapt$pred)))

xtable(summ_stat_acc, digits = 4)

summ_stat_tree <- rbind(c(as.vector(gauss_stat$tree),as.vector(gauss_stat_adapt$tree)),
                        c(as.vector(t_stat$tree),as.vector(t_stat_adapt$tree)),
                        c(as.vector(clayton_stat$tree),as.vector(clayton_stat_adapt$tree)),
                        c(as.vector(gumbel_stat$tree),as.vector(gumbel_stat_adapt$tree)),
                        c(as.vector(frank_stat$tree),as.vector(frank_stat_adapt$tree)))

xtable(summ_stat_tree, digits = 2)

(gauss$like + labs(title="Gaussian (without adaption)") + geom_hline(yintercept = gauss_like_true, linetype = 2) + gauss_adapt$like + labs(title="Gaussian (with adaption)") + geom_hline(yintercept = gauss_like_true, linetype = 2)) / 
  (t$like + labs(title="Student-t (without adaption)") + geom_hline(yintercept = t_like_true, linetype = 2) + t_adapt$like + labs(title="Student-t (with adaption)") + geom_hline(yintercept = t_like_true, linetype = 2)) /
  (clayton$like + labs(title="Clayton (without adaption)")+ geom_hline(yintercept = clayton_like_true, linetype = 2) + clayton_adapt$like + labs(title="Clayton (with adaption)")+ geom_hline(yintercept = clayton_like_true, linetype = 2)) /
  (gumbel$like + labs(title="Gumbel (without adaption)")+ geom_hline(yintercept = gumbel_like_true, linetype = 2) + gumbel_adapt$like + labs(title="Gumbel (with adaption)")+ geom_hline(yintercept = gumbel_like_true, linetype = 2)) /
  (frank$like + labs(title="Frank (without adaption)")+ geom_hline(yintercept = frank_like_true, linetype = 2) + frank_adapt$like + labs(title="Frank (with adaption)")+ geom_hline(yintercept = frank_like_true, linetype = 2))

(gauss$pred + labs(title="Gaussian (without adaption)") + ylim(-0.2,1) + gauss_adapt$pred + labs(title="Gaussian (with adaption)") + ylim(-0.2,1)) / 
  (t$pred + labs(title="Student-t (without adaption)") + ylim(-0.2,1) + t_adapt$pred + labs(title="Student-t (with adaption)") + ylim(-0.2,1)) /
  (clayton$pred + labs(title="Clayton (without adaption)") + ylim(0,1) + clayton_adapt$pred + labs(title="Clayton (with adaption)") + ylim(0,1)) /
  (gumbel$pred + labs(title="Gumbel (without adaption)") + ylim(0,1) + gumbel_adapt$pred + labs(title="Gumbel (with adaption)") + ylim(0,1)) /
  (frank$pred + labs(title="Frank (without adaption)") + ylim(-0.1,1) + frank_adapt$pred + labs(title="Frank (with adaption)") + ylim(-0.1,1))

# (gauss$nterm + labs(title="Gaussian (without adaption)") + geom_hline(yintercept = 3, linetype = 2) + ylim(0, 8) + gauss_adapt$nterm + labs(title="Gaussian (with adaption)") + geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8)) / 
#   (t$nterm + labs(title="Student-t (without adaption)") + geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8) + t_adapt$nterm + labs(title="Student-t (with adaption)") + geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8)) /
#   (clayton$nterm + labs(title="Clayton (without adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8) + clayton_adapt$nterm + labs(title="Clayton (with adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8)) /
#   (gumbel$nterm + labs(title="Gumbel (without adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8) + gumbel_adapt$nterm + labs(title="Gumbel (with adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8)) /
#   (frank$nterm + labs(title="Frank (without adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8) + frank_adapt$nterm + labs(title="Frank (with adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8))
# 
# (gauss$depth + labs(title="Gaussian (without adaption)") + geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + gauss_adapt$depth + labs(title="Gaussian (with adaption)") + geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6)) / 
#   (t$depth + labs(title="Student-t (without adaption)") + geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + t_adapt$depth + labs(title="Student-t (with adaption)") + geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6)) /
#   (clayton$depth + labs(title="Clayton (without adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + clayton_adapt$depth + labs(title="Clayton (with adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6)) /
#   (gumbel$depth + labs(title="Gumbel (without adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + gumbel_adapt$depth + labs(title="Gumbel (with adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6)) /
#   (frank$depth + labs(title="Frank (without adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + frank_adapt$depth + labs(title="Frank (with adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6))
