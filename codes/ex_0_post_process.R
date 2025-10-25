# packages
library(VineCopula)
library(dplyr)
library(ggplot2)
library(patchwork)

################################################################################
# data generation
################################################################################
set.seed(1e3)

n <- 200
R <- 100
X_obs <- lapply(1:R, function(x)matrix(runif(n), ncol = 1))

# normalise predictors 
data_norm <- function(X_obs){
  X_obs.norm <- as.data.frame(apply(X_obs, 2, \(x) (x - min(x))/(max(x) - min(x))))
  X_obs.norm <- as.matrix(X_obs.norm)
  rownames(X_obs.norm) <- 1:nrow(X_obs)
  return(X_obs.norm)
}

X_obs.norm <- lapply(1:R, function(i)data_norm(X_obs[[i]]))

# tau with tree structure
tree_tau_func <- function(X_obs){
  tau_true_1 <- rep(0,nrow(X_obs))
  tau_true_1[X_obs<0.33] <- 0.3
  tau_true_1[(X_obs>=0.33)&(X_obs<0.66)] <- 0.7
  tau_true_1[(X_obs>=0.66)] <- 0.4
  
  tau_true_1 <- matrix(tau_true_1, ncol = 1)
}

tau_true_1 <- lapply(1:R, function(i){tree_tau_func(X_obs[[i]])})

# periodoic
tau_true_2 <- lapply(1:R, function(i){0.2*sin(2*pi*X_obs[[i]]) + 0.5})

plot(X_obs[[1]], tau_true_1[[1]], xlab = "Observations", ylab = "tau")

plot(X_obs[[1]], tau_true_2[[1]], xlab = "Observations", ylab = "tau")

# gauss # sin(tau*pi/2)

for (i in 1:2) {
  assign(paste0("copula_uu_gauss_",i), lapply(1:R,function(k){BiCopSim(n, family = 1, par = sin(get(paste0("tau_true_",i))[[k]] * pi/2))}))
}

# t # sin(tau*pi/2)

for (i in 1:2) {
  assign(paste0("copula_uu_t_",i), lapply(1:R, function(k){BiCopSim(n, family = 2, par = sin(get(paste0("tau_true_",i))[[k]] * pi/2), par2 = 3)}))
}

# gumbel # 1/(1-tau)

for (i in 1:2) {
  assign(paste0("copula_uu_gumbel_",i), lapply(1:R, function(k){BiCopSim(n, family = 4, par = 1/(1-get(paste0("tau_true_",i))[[k]]))}))
}

# clayton # (2*tau)/(1-tau)

for (i in 1:2) {
  assign(paste0("copula_uu_clayton_",i), lapply(1:R, function(k){BiCopSim(n, family = 3, par = (2*get(paste0("tau_true_",i))[[k]])/(1-get(paste0("tau_true_",i))[[k]]))}))
}

# frank # numerical

for (i in 1:2) {
  assign(paste0("copula_uu_frank_",i), lapply(1:R, function(k){BiCopSim(n, family = 5, par = BiCopTau2Par(5, get(paste0("tau_true_",i))[[k]]))}))
}

# mcmc params

n.chain_par <- 4
n.iter_par <- 3000
n.discard <- 1500
n.thin <- 1

incl.split_par <- TRUE
cont.unif_par <- TRUE

################################################################################
# source files
################################################################################

source('MCMC_BART_copula.R')
source('import_functions.R')

moves.prob_par <- c(0.1, 0.4, 0.25, 0.25)
lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) 

n.tree <- 1
test_case <- 1

data_pred <- function(pred_list, x_list){
  
  pred_cond_mod <- list()
  for (i in c(1:length(x_list))) {
    pred_cond <- data.frame("obs" = rep(x_list[[i]], each = (n.chain_par * n.iter_par)))
    pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))[[i]]), each = (n.chain_par * (n.iter_par)))
    pred_cond$y = BiCopPar2Tau(1,link_gauss(as.vector(pred_list[[i]])))
    pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
    pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
    pred_cond$set = rep(i, (n.iter_par*n.chain_par*n))
    
    pred_cond_mod[[i]] = pred_cond %>%
      filter(idx > n.discard) %>%
      group_by(obs, chain, theta_true) %>%
      summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
    
  }
  return(pred_cond_mod)
}


# gauss

if(T){
  
  load(paste0("gauss_mcmc_",test_case,"_tree_",n.tree, ".Rdata"))
  
  summ_list <- BART_summary(get(paste0("gauss_mcmc_",test_case,"_tree_",n.tree)), X_obs.norm, others = T)
  
  rm(list = paste0("gauss_mcmc_",test_case,"_tree_",n.tree))
  gc()
  
  pred_cond_mod_list <- data_pred(summ_list$pred, X_obs.norm)
  
  pred_cond_case_1 <- data.frame("obs" = rep(X_obs.norm[[1]], each = (n.chain_par * n.iter_par)))
  pred_cond_case_1$theta_true = rep((get(paste0("tau_true_",test_case))[[1]]), each = (n.chain_par * (n.iter_par)))
  pred_cond_case_1$y = BiCopPar2Tau(1,link_gauss(as.vector(summ_list$pred[[1]])))
  pred_cond_case_1$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond_case_1$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_case_1_mod = pred_cond_case_1 %>%
    filter(idx > n.discard) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_case_1_mod_avg = pred_cond_case_1_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_case_1_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df_case_1 <-data.frame("nn" = apply(summ_list$pred[[1]], 1, function(x)loglik_gauss(link_gauss(x), get(paste0("copula_uu_gauss_",test_case))[[1]][,1], get(paste0("copula_uu_gauss_",test_case))[[1]][,2])))
  like_df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df_case_1 %>%
    filter(idx > n.discard) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  pl_like
  
  gauss_like_true <- loglik_gauss((sin(get(paste0("tau_true_",test_case))[[1]] * pi/2)), get(paste0("copula_uu_gauss_",test_case))[[1]][,1] , get(paste0("copula_uu_gauss_",test_case))[[1]][,2])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    nt_lb.df_case_1 <- summ_list$nterm[[1]]
    nt_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    nt_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    
    pl_nl <- nt_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "nterm"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_nl
    # depth
    
    depth_lb.df_case_1 <- summ_list$depth[[1]]
    depth_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    depth_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    pl_dp <- depth_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "depth"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_dp
    # acceptance
    
    acc_lb.df <- do.call(rbind,summ_list$acc)
    nt_lb.df <- do.call(rbind, summ_list$nterm)
    depth_lb.df <- do.call(rbind, summ_list$depth)
    
    
    tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.discard]), mean(depth_lb.df$nn[depth_lb.df$idx > n.discard]), mean(acc_lb.df$nn[acc_lb.df$idx > n.discard]==TRUE))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc")
  }

  
  gauss_dat_list <- list("all" = summ_list, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    gauss_stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    gauss_plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    gauss_stat_list <- list("pred" = pred_cond_summary)
    gauss_plot_list <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(gauss_dat_list, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  save(gauss_stat_list, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  save(gauss_plot_list, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, gauss_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, gauss_dat_list, gauss_stat_list, gauss_plot_list)
gc()

if(T){
  
  load(paste0("gauss_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata"))
  
  summ_list <- BART_summary(get(paste0("gauss_mcmc_",test_case,"_tree_",n.tree,"_adapt")), X_obs.norm, others = T)
  
  rm(list = paste0("gauss_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  gc()
  
  pred_cond_mod_list <- data_pred(summ_list$pred, X_obs.norm)
  
  pred_cond_case_1 <- data.frame("obs" = rep(X_obs.norm[[1]], each = (n.chain_par * n.iter_par)))
  pred_cond_case_1$theta_true = rep((get(paste0("tau_true_",test_case))[[1]]), each = (n.chain_par * (n.iter_par)))
  pred_cond_case_1$y = BiCopPar2Tau(1,link_gauss(as.vector(summ_list$pred[[1]])))
  pred_cond_case_1$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond_case_1$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_case_1_mod = pred_cond_case_1 %>%
    filter(idx > n.discard) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_case_1_mod_avg = pred_cond_case_1_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_case_1_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df_case_1 <-data.frame("nn" = apply(summ_list$pred[[1]], 1, function(x)loglik_gauss(link_gauss(x), get(paste0("copula_uu_gauss_",test_case))[[1]][,1], get(paste0("copula_uu_gauss_",test_case))[[1]][,2])))
  like_df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df_case_1 %>%
    filter(idx > n.discard) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  pl_like
  
  gauss_like_true <- loglik_gauss((sin(get(paste0("tau_true_",test_case))[[1]] * pi/2)), get(paste0("copula_uu_gauss_",test_case))[[1]][,1] , get(paste0("copula_uu_gauss_",test_case))[[1]][,2])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    nt_lb.df_case_1 <- summ_list$nterm[[1]]
    nt_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    nt_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    
    pl_nl <- nt_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "nterm"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_nl
    # depth
    
    depth_lb.df_case_1 <- summ_list$depth[[1]]
    depth_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    depth_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    pl_dp <- depth_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "depth"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_dp
    # acceptance
    
    acc_lb.df <- do.call(rbind,summ_list$acc)
    nt_lb.df <- do.call(rbind, summ_list$nterm)
    depth_lb.df <- do.call(rbind, summ_list$depth)
    
    
    tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.discard]), mean(depth_lb.df$nn[depth_lb.df$idx > n.discard]), mean(acc_lb.df$nn[acc_lb.df$idx > n.discard]==TRUE))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc")
  }
  
  
  gauss_dat_list_adapt <- list("all" = summ_list, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    gauss_stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    gauss_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    gauss_stat_list_adapt <- list("pred" = pred_cond_summary)
    gauss_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(gauss_dat_list_adapt, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  save(gauss_stat_list_adapt, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  save(gauss_plot_list_adapt, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, gauss_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, gauss_dat_list_adapt, gauss_stat_list_adapt, gauss_plot_list_adapt)
gc()

# t

if(T){
  
  load(paste0("t_mcmc_",test_case,"_tree_",n.tree, ".Rdata"))
  
  summ_list <- BART_summary(get(paste0("t_mcmc_",test_case,"_tree_",n.tree)), X_obs.norm, others = T)
  
  rm(list = paste0("t_mcmc_",test_case,"_tree_",n.tree))
  gc()
  
  pred_cond_mod_list <- data_pred(summ_list$pred, X_obs.norm)
  
  pred_cond_case_1 <- data.frame("obs" = rep(X_obs.norm[[1]], each = (n.chain_par * n.iter_par)))
  pred_cond_case_1$theta_true = rep((get(paste0("tau_true_",test_case))[[1]]), each = (n.chain_par * (n.iter_par)))
  pred_cond_case_1$y = BiCopPar2Tau(2,link_t(as.vector(summ_list$pred[[1]])))
  pred_cond_case_1$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond_case_1$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_case_1_mod = pred_cond_case_1 %>%
    filter(idx > n.discard) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_case_1_mod_avg = pred_cond_case_1_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_case_1_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df_case_1 <-data.frame("nn" = apply(summ_list$pred[[1]], 1, function(x)loglik_t(link_t(x), get(paste0("copula_uu_t_",test_case))[[1]][,1], get(paste0("copula_uu_t_",test_case))[[1]][,2])))
  like_df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df_case_1 %>%
    filter(idx > n.discard) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  pl_like
  
  t_like_true <- loglik_t((sin(get(paste0("tau_true_",test_case))[[1]] * pi/2)), get(paste0("copula_uu_t_",test_case))[[1]][,1] , get(paste0("copula_uu_t_",test_case))[[1]][,2])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    nt_lb.df_case_1 <- summ_list$nterm[[1]]
    nt_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    nt_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    
    pl_nl <- nt_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "nterm"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_nl
    # depth
    
    depth_lb.df_case_1 <- summ_list$depth[[1]]
    depth_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    depth_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    pl_dp <- depth_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "depth"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_dp
    # acceptance
    
    acc_lb.df <- do.call(rbind,summ_list$acc)
    nt_lb.df <- do.call(rbind, summ_list$nterm)
    depth_lb.df <- do.call(rbind, summ_list$depth)
    
    
    tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.discard]), mean(depth_lb.df$nn[depth_lb.df$idx > n.discard]), mean(acc_lb.df$nn[acc_lb.df$idx > n.discard]==TRUE))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc")
  }
  
  
  t_dat_list <- list("all" = summ_list, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    t_stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    t_plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    t_stat_list <- list("pred" = pred_cond_summary)
    t_plot_list <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(t_dat_list, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  save(t_stat_list, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  save(t_plot_list, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, t_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, t_dat_list, t_stat_list, t_plot_list)
gc()

if(T){
  
  load(paste0("t_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata"))
  
  summ_list <- BART_summary(get(paste0("t_mcmc_",test_case,"_tree_",n.tree,"_adapt")), X_obs.norm, others = T)
  
  rm(list = paste0("t_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  gc()
  
  pred_cond_mod_list <- data_pred(summ_list$pred, X_obs.norm)
  
  pred_cond_case_1 <- data.frame("obs" = rep(X_obs.norm[[1]], each = (n.chain_par * n.iter_par)))
  pred_cond_case_1$theta_true = rep((get(paste0("tau_true_",test_case))[[1]]), each = (n.chain_par * (n.iter_par)))
  pred_cond_case_1$y = BiCopPar2Tau(2,link_t(as.vector(summ_list$pred[[1]])))
  pred_cond_case_1$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond_case_1$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_case_1_mod = pred_cond_case_1 %>%
    filter(idx > n.discard) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_case_1_mod_avg = pred_cond_case_1_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_case_1_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df_case_1 <-data.frame("nn" = apply(summ_list$pred[[1]], 1, function(x)loglik_t(link_t(x), get(paste0("copula_uu_t_",test_case))[[1]][,1], get(paste0("copula_uu_t_",test_case))[[1]][,2])))
  like_df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df_case_1 %>%
    filter(idx > n.discard) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  pl_like
  
  t_like_true <- loglik_t((sin(get(paste0("tau_true_",test_case))[[1]] * pi/2)), get(paste0("copula_uu_t_",test_case))[[1]][,1] , get(paste0("copula_uu_t_",test_case))[[1]][,2])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    nt_lb.df_case_1 <- summ_list$nterm[[1]]
    nt_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    nt_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    
    pl_nl <- nt_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "nterm"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_nl
    # depth
    
    depth_lb.df_case_1 <- summ_list$depth[[1]]
    depth_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    depth_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    pl_dp <- depth_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "depth"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_dp
    # acceptance
    
    acc_lb.df <- do.call(rbind,summ_list$acc)
    nt_lb.df <- do.call(rbind, summ_list$nterm)
    depth_lb.df <- do.call(rbind, summ_list$depth)
    
    
    tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.discard]), mean(depth_lb.df$nn[depth_lb.df$idx > n.discard]), mean(acc_lb.df$nn[acc_lb.df$idx > n.discard]==TRUE))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc")
  }
  
  
  t_dat_list_adapt <- list("all" = summ_list, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    t_stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    t_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    t_stat_list_adapt <- list("pred" = pred_cond_summary)
    t_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(t_dat_list_adapt, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  save(t_stat_list_adapt, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  save(t_plot_list_adapt, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, t_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, t_dat_list_adapt, t_stat_list_adapt, t_plot_list_adapt)
gc()

# clayton

if(T){
  
  load(paste0("clayton_mcmc_",test_case,"_tree_",n.tree, ".Rdata"))
  
  summ_list <- BART_summary(get(paste0("clayton_mcmc_",test_case,"_tree_",n.tree)), X_obs.norm, others = T)
  
  rm(list = paste0("clayton_mcmc_",test_case,"_tree_",n.tree))
  gc()
  
  pred_cond_mod_list <- data_pred(summ_list$pred, X_obs.norm)
  
  pred_cond_case_1 <- data.frame("obs" = rep(X_obs.norm[[1]], each = (n.chain_par * n.iter_par)))
  pred_cond_case_1$theta_true = rep((get(paste0("tau_true_",test_case))[[1]]), each = (n.chain_par * (n.iter_par)))
  pred_cond_case_1$y = BiCopPar2Tau(3,link_clayton(as.vector(summ_list$pred[[1]])))
  pred_cond_case_1$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond_case_1$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_case_1_mod = pred_cond_case_1 %>%
    filter(idx > n.discard) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_case_1_mod_avg = pred_cond_case_1_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_case_1_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df_case_1 <-data.frame("nn" = apply(summ_list$pred[[1]], 1, function(x)loglik_clayton(link_clayton(x), get(paste0("copula_uu_clayton_",test_case))[[1]][,1], get(paste0("copula_uu_clayton_",test_case))[[1]][,2])))
  like_df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df_case_1 %>%
    filter(idx > n.discard) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  pl_like
  
  clayton_like_true <- loglik_clayton((sin(get(paste0("tau_true_",test_case))[[1]] * pi/2)), get(paste0("copula_uu_clayton_",test_case))[[1]][,1] , get(paste0("copula_uu_clayton_",test_case))[[1]][,2])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    nt_lb.df_case_1 <- summ_list$nterm[[1]]
    nt_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    nt_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    
    pl_nl <- nt_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "nterm"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_nl
    # depth
    
    depth_lb.df_case_1 <- summ_list$depth[[1]]
    depth_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    depth_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    pl_dp <- depth_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "depth"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_dp
    # acceptance
    
    acc_lb.df <- do.call(rbind,summ_list$acc)
    nt_lb.df <- do.call(rbind, summ_list$nterm)
    depth_lb.df <- do.call(rbind, summ_list$depth)
    
    
    tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.discard]), mean(depth_lb.df$nn[depth_lb.df$idx > n.discard]), mean(acc_lb.df$nn[acc_lb.df$idx > n.discard]==TRUE))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc")
  }
  
  
  clayton_dat_list <- list("all" = summ_list, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    clayton_stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    clayton_plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    clayton_stat_list <- list("pred" = pred_cond_summary)
    clayton_plot_list <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(clayton_dat_list, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  save(clayton_stat_list, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  save(clayton_plot_list, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, clayton_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, clayton_dat_list, clayton_stat_list, clayton_plot_list)
gc()

if(T){
  
  load(paste0("clayton_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata"))
  
  summ_list <- BART_summary(get(paste0("clayton_mcmc_",test_case,"_tree_",n.tree,"_adapt")), X_obs.norm, others = T)
  
  rm(list = paste0("clayton_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  gc()
  
  pred_cond_mod_list <- data_pred(summ_list$pred, X_obs.norm)
  
  pred_cond_case_1 <- data.frame("obs" = rep(X_obs.norm[[1]], each = (n.chain_par * n.iter_par)))
  pred_cond_case_1$theta_true = rep((get(paste0("tau_true_",test_case))[[1]]), each = (n.chain_par * (n.iter_par)))
  pred_cond_case_1$y = BiCopPar2Tau(3,link_clayton(as.vector(summ_list$pred[[1]])))
  pred_cond_case_1$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond_case_1$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_case_1_mod = pred_cond_case_1 %>%
    filter(idx > n.discard) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_case_1_mod_avg = pred_cond_case_1_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_case_1_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df_case_1 <-data.frame("nn" = apply(summ_list$pred[[1]], 1, function(x)loglik_clayton(link_clayton(x), get(paste0("copula_uu_clayton_",test_case))[[1]][,1], get(paste0("copula_uu_clayton_",test_case))[[1]][,2])))
  like_df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df_case_1 %>%
    filter(idx > n.discard) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  pl_like
  
  clayton_like_true <- loglik_clayton((sin(get(paste0("tau_true_",test_case))[[1]] * pi/2)), get(paste0("copula_uu_clayton_",test_case))[[1]][,1] , get(paste0("copula_uu_clayton_",test_case))[[1]][,2])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    nt_lb.df_case_1 <- summ_list$nterm[[1]]
    nt_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    nt_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    
    pl_nl <- nt_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "nterm"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_nl
    # depth
    
    depth_lb.df_case_1 <- summ_list$depth[[1]]
    depth_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    depth_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    pl_dp <- depth_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "depth"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_dp
    # acceptance
    
    acc_lb.df <- do.call(rbind,summ_list$acc)
    nt_lb.df <- do.call(rbind, summ_list$nterm)
    depth_lb.df <- do.call(rbind, summ_list$depth)
    
    
    tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.discard]), mean(depth_lb.df$nn[depth_lb.df$idx > n.discard]), mean(acc_lb.df$nn[acc_lb.df$idx > n.discard]==TRUE))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc")
  }
  
  
  clayton_dat_list_adapt <- list("all" = summ_list, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    clayton_stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    clayton_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    clayton_stat_list_adapt <- list("pred" = pred_cond_summary)
    clayton_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(clayton_dat_list_adapt, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  save(clayton_stat_list_adapt, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  save(clayton_plot_list_adapt, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, clayton_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, clayton_dat_list_adapt, clayton_stat_list_adapt, clayton_plot_list_adapt)
gc()

# gumbel

if(T){
  
  load(paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, ".Rdata"))
  
  summ_list <- BART_summary(get(paste0("gumbel_mcmc_",test_case,"_tree_",n.tree)), X_obs.norm, others = T)
  
  rm(list = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree))
  gc()
  
  pred_cond_mod_list <- data_pred(summ_list$pred, X_obs.norm)
  
  pred_cond_case_1 <- data.frame("obs" = rep(X_obs.norm[[1]], each = (n.chain_par * n.iter_par)))
  pred_cond_case_1$theta_true = rep((get(paste0("tau_true_",test_case))[[1]]), each = (n.chain_par * (n.iter_par)))
  pred_cond_case_1$y = BiCopPar2Tau(4,link_gumbel(as.vector(summ_list$pred[[1]])))
  pred_cond_case_1$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond_case_1$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_case_1_mod = pred_cond_case_1 %>%
    filter(idx > n.discard) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_case_1_mod_avg = pred_cond_case_1_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_case_1_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df_case_1 <-data.frame("nn" = apply(summ_list$pred[[1]], 1, function(x)loglik_gumbel(link_gumbel(x), get(paste0("copula_uu_gumbel_",test_case))[[1]][,1], get(paste0("copula_uu_gumbel_",test_case))[[1]][,2])))
  like_df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df_case_1 %>%
    filter(idx > n.discard) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  pl_like
  
  gumbel_like_true <- loglik_gumbel((sin(get(paste0("tau_true_",test_case))[[1]] * pi/2)), get(paste0("copula_uu_gumbel_",test_case))[[1]][,1] , get(paste0("copula_uu_gumbel_",test_case))[[1]][,2])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    nt_lb.df_case_1 <- summ_list$nterm[[1]]
    nt_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    nt_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    
    pl_nl <- nt_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "nterm"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_nl
    # depth
    
    depth_lb.df_case_1 <- summ_list$depth[[1]]
    depth_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    depth_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    pl_dp <- depth_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "depth"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_dp
    # acceptance
    
    acc_lb.df <- do.call(rbind,summ_list$acc)
    nt_lb.df <- do.call(rbind, summ_list$nterm)
    depth_lb.df <- do.call(rbind, summ_list$depth)
    
    
    tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.discard]), mean(depth_lb.df$nn[depth_lb.df$idx > n.discard]), mean(acc_lb.df$nn[acc_lb.df$idx > n.discard]==TRUE))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc")
  }
  
  
  gumbel_dat_list <- list("all" = summ_list, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    gumbel_stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    gumbel_plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    gumbel_stat_list <- list("pred" = pred_cond_summary)
    gumbel_plot_list <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(gumbel_dat_list, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  save(gumbel_stat_list, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  save(gumbel_plot_list, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, gumbel_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, gumbel_dat_list, gumbel_stat_list, gumbel_plot_list)
gc()

if(T){
  
  load(paste0("gumbel_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata"))
  
  summ_list <- BART_summary(get(paste0("gumbel_mcmc_",test_case,"_tree_",n.tree,"_adapt")), X_obs.norm, others = T)
  
  rm(list = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  gc()
  
  pred_cond_mod_list <- data_pred(summ_list$pred, X_obs.norm)
  
  pred_cond_case_1 <- data.frame("obs" = rep(X_obs.norm[[1]], each = (n.chain_par * n.iter_par)))
  pred_cond_case_1$theta_true = rep((get(paste0("tau_true_",test_case))[[1]]), each = (n.chain_par * (n.iter_par)))
  pred_cond_case_1$y = BiCopPar2Tau(4,link_gumbel(as.vector(summ_list$pred[[1]])))
  pred_cond_case_1$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond_case_1$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_case_1_mod = pred_cond_case_1 %>%
    filter(idx > n.discard) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_case_1_mod_avg = pred_cond_case_1_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_case_1_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df_case_1 <-data.frame("nn" = apply(summ_list$pred[[1]], 1, function(x)loglik_gumbel(link_gumbel(x), get(paste0("copula_uu_gumbel_",test_case))[[1]][,1], get(paste0("copula_uu_gumbel_",test_case))[[1]][,2])))
  like_df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df_case_1 %>%
    filter(idx > n.discard) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  pl_like
  
  gumbel_like_true <- loglik_gumbel((sin(get(paste0("tau_true_",test_case))[[1]] * pi/2)), get(paste0("copula_uu_gumbel_",test_case))[[1]][,1] , get(paste0("copula_uu_gumbel_",test_case))[[1]][,2])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    nt_lb.df_case_1 <- summ_list$nterm[[1]]
    nt_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    nt_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    
    pl_nl <- nt_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "nterm"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_nl
    # depth
    
    depth_lb.df_case_1 <- summ_list$depth[[1]]
    depth_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    depth_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    pl_dp <- depth_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "depth"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_dp
    # acceptance
    
    acc_lb.df <- do.call(rbind,summ_list$acc)
    nt_lb.df <- do.call(rbind, summ_list$nterm)
    depth_lb.df <- do.call(rbind, summ_list$depth)
    
    
    tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.discard]), mean(depth_lb.df$nn[depth_lb.df$idx > n.discard]), mean(acc_lb.df$nn[acc_lb.df$idx > n.discard]==TRUE))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc")
  }
  
  
  gumbel_dat_list_adapt <- list("all" = summ_list, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    gumbel_stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    gumbel_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    gumbel_stat_list_adapt <- list("pred" = pred_cond_summary)
    gumbel_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(gumbel_dat_list_adapt, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  save(gumbel_stat_list_adapt, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  save(gumbel_plot_list_adapt, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, gumbel_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, gumbel_dat_list_adapt, gumbel_stat_list_adapt, gumbel_plot_list_adapt)
gc()

# frank

if(T){
  
  load(paste0("frank_mcmc_",test_case,"_tree_",n.tree, ".Rdata"))
  
  summ_list <- BART_summary(get(paste0("frank_mcmc_",test_case,"_tree_",n.tree)), X_obs.norm, others = T)
  
  rm(list = paste0("frank_mcmc_",test_case,"_tree_",n.tree))
  gc()
  
  pred_cond_mod_list <- data_pred(summ_list$pred, X_obs.norm)
  
  pred_cond_case_1 <- data.frame("obs" = rep(X_obs.norm[[1]], each = (n.chain_par * n.iter_par)))
  pred_cond_case_1$theta_true = rep((get(paste0("tau_true_",test_case))[[1]]), each = (n.chain_par * (n.iter_par)))
  pred_cond_case_1$y = BiCopPar2Tau(5,link_frank(as.vector(summ_list$pred[[1]])))
  pred_cond_case_1$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond_case_1$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_case_1_mod = pred_cond_case_1 %>%
    filter(idx > n.discard) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_case_1_mod_avg = pred_cond_case_1_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_case_1_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df_case_1 <-data.frame("nn" = apply(summ_list$pred[[1]], 1, function(x)loglik_frank(link_frank(x), get(paste0("copula_uu_frank_",test_case))[[1]][,1], get(paste0("copula_uu_frank_",test_case))[[1]][,2])))
  like_df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df_case_1 %>%
    filter(idx > n.discard) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  pl_like
  
  frank_like_true <- loglik_frank((sin(get(paste0("tau_true_",test_case))[[1]] * pi/2)), get(paste0("copula_uu_frank_",test_case))[[1]][,1] , get(paste0("copula_uu_frank_",test_case))[[1]][,2])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    nt_lb.df_case_1 <- summ_list$nterm[[1]]
    nt_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    nt_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    
    pl_nl <- nt_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "nterm"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_nl
    # depth
    
    depth_lb.df_case_1 <- summ_list$depth[[1]]
    depth_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    depth_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    pl_dp <- depth_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "depth"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_dp
    # acceptance
    
    acc_lb.df <- do.call(rbind,summ_list$acc)
    nt_lb.df <- do.call(rbind, summ_list$nterm)
    depth_lb.df <- do.call(rbind, summ_list$depth)
    
    
    tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.discard]), mean(depth_lb.df$nn[depth_lb.df$idx > n.discard]), mean(acc_lb.df$nn[acc_lb.df$idx > n.discard]==TRUE))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc")
  }
  
  
  frank_dat_list <- list("all" = summ_list, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    frank_stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    frank_plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    frank_stat_list <- list("pred" = pred_cond_summary)
    frank_plot_list <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(frank_dat_list, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  save(frank_stat_list, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  save(frank_plot_list, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, frank_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, frank_dat_list, frank_stat_list, frank_plot_list)
gc()

if(T){
  
  load(paste0("frank_mcmc_",test_case,"_tree_",n.tree,"_adapt", ".Rdata"))
  
  summ_list <- BART_summary(get(paste0("frank_mcmc_",test_case,"_tree_",n.tree,"_adapt")), X_obs.norm, others = T)
  
  rm(list = paste0("frank_mcmc_",test_case,"_tree_",n.tree,"_adapt"))
  gc()
  
  pred_cond_mod_list <- data_pred(summ_list$pred, X_obs.norm)
  
  pred_cond_case_1 <- data.frame("obs" = rep(X_obs.norm[[1]], each = (n.chain_par * n.iter_par)))
  pred_cond_case_1$theta_true = rep((get(paste0("tau_true_",test_case))[[1]]), each = (n.chain_par * (n.iter_par)))
  pred_cond_case_1$y = BiCopPar2Tau(5,link_frank(as.vector(summ_list$pred[[1]])))
  pred_cond_case_1$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
  pred_cond_case_1$idx = rep(rep(1:n.iter_par, n.chain_par),n)
  
  pred_cond_case_1_mod = pred_cond_case_1 %>%
    filter(idx > n.discard) %>%
    group_by(obs, chain, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  pred_cond_case_1_mod_avg = pred_cond_case_1_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  
  pl_pred <- ggplot(pred_cond_case_1_mod_avg) +
    geom_point(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    xlab('X') +
    ylab('estimated tau') +
    theme_classic()
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(theta_mean), theta_q975 = mean(theta_q975), theta_q025 = mean(theta_q025)) 
  
  pred_cond_stat = pred_cond_mod_avg %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # like
  
  like_df_case_1 <-data.frame("nn" = apply(summ_list$pred[[1]], 1, function(x)loglik_frank(link_frank(x), get(paste0("copula_uu_frank_",test_case))[[1]][,1], get(paste0("copula_uu_frank_",test_case))[[1]][,2])))
  like_df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
  
  
  pl_like <- like_df_case_1 %>%
    filter(idx > n.discard) %>%
    ggplot(aes(x = idx, y = nn, color = factor(chain))) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Log-likelihood"
    ) +
    guides(color = "none") +
    theme_minimal()
  
  pl_like
  
  frank_like_true <- loglik_frank((sin(get(paste0("tau_true_",test_case))[[1]] * pi/2)), get(paste0("copula_uu_frank_",test_case))[[1]][,1] , get(paste0("copula_uu_frank_",test_case))[[1]][,2])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    nt_lb.df_case_1 <- summ_list$nterm[[1]]
    nt_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    nt_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    
    pl_nl <- nt_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "nterm"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_nl
    # depth
    
    depth_lb.df_case_1 <- summ_list$depth[[1]]
    depth_lb.df_case_1$idx <- rep(1:n.iter_par, n.chain_par)
    depth_lb.df_case_1$chain <- rep(1:n.chain_par, each = n.iter_par)
    
    pl_dp <- depth_lb.df_case_1 %>%
      # filter(idx > n.discard) %>%
      ggplot(aes(x = idx, y = nn, color = factor(chain))) +
      geom_line() +
      labs(
        x = "Iteration",
        y = "depth"
      ) +
      guides(color = "none") +
      theme_minimal()
    
    pl_dp
    # acceptance
    
    acc_lb.df <- do.call(rbind,summ_list$acc)
    nt_lb.df <- do.call(rbind, summ_list$nterm)
    depth_lb.df <- do.call(rbind, summ_list$depth)
    
    
    tree_stat_summary = c(mean(nt_lb.df$nn[nt_lb.df$idx > n.discard]), mean(depth_lb.df$nn[depth_lb.df$idx > n.discard]), mean(acc_lb.df$nn[acc_lb.df$idx > n.discard]==TRUE))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc")
  }
  
  
  frank_dat_list_adapt <- list("all" = summ_list, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    frank_stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    frank_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    frank_stat_list_adapt <- list("pred" = pred_cond_summary)
    frank_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(frank_dat_list_adapt, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  save(frank_stat_list_adapt, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  save(frank_plot_list_adapt, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, frank_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, frank_dat_list_adapt, frank_stat_list_adapt, frank_plot_list_adapt)
gc()

if(F){
  gauss <- get(load(paste0("gauss_mcmc_",test_case,"_tree_",n.tree,"_plot", ".Rdata")))
  gauss_stat <- get(load(paste0("gauss_mcmc_",test_case,"_tree_",n.tree,"_stat", ".Rdata")))
  
  gauss_adapt <- get(load(paste0("gauss_mcmc_",test_case,"_tree_",n.tree,"_plot_adapt", ".Rdata")))
  gauss_stat_adapt <- get(load(paste0("gauss_mcmc_",test_case,"_tree_",n.tree,"_stat_adapt", ".Rdata")))
  
  gauss_like_true <- loglik_gauss((sin(get(paste0("tau_true_",test_case)) * pi/2)), get(paste0("copula_uu_gauss_",test_case))[,1] , get(paste0("copula_uu_gauss_",test_case))[,2])
  
  
  t <- get(load(paste0("t_mcmc_",test_case,"_tree_",n.tree,"_plot", ".Rdata")))
  t_stat <- get(load(paste0("t_mcmc_",test_case,"_tree_",n.tree,"_stat", ".Rdata")))
  
  t_adapt <- get(load(paste0("t_mcmc_",test_case,"_tree_",n.tree,"_plot_adapt", ".Rdata")))
  t_stat_adapt <- get(load(paste0("t_mcmc_",test_case,"_tree_",n.tree,"_stat_adapt", ".Rdata")))
  
  t_like_true <- loglik_t((sin(get(paste0("tau_true_",test_case)) * pi/2)), get(paste0("copula_uu_t_",test_case))[,1] , get(paste0("copula_uu_t_",test_case))[,2])
  
  
  clayton <- get(load(paste0("clayton_mcmc_",test_case,"_tree_",n.tree,"_plot", ".Rdata")))
  clayton_stat <- get(load(paste0("clayton_mcmc_",test_case,"_tree_",n.tree,"_stat", ".Rdata")))
  
  clayton_adapt <- get(load(paste0("clayton_mcmc_",test_case,"_tree_",n.tree,"_plot_adapt", ".Rdata")))
  clayton_stat_adapt <- get(load(paste0("clayton_mcmc_",test_case,"_tree_",n.tree,"_stat_adapt", ".Rdata")))
  
  clayton_like_true <- loglik_clayton((2*get(paste0("tau_true_",test_case)) / (1-get(paste0("tau_true_",test_case)))), get(paste0("copula_uu_clayton_",test_case))[,1] , get(paste0("copula_uu_clayton_",test_case))[,2])
  
  
  gumbel <- get(load(paste0("gumbel_mcmc_",test_case,"_tree_",n.tree,"_plot", ".Rdata")))
  gumbel_stat <- get(load(paste0("gumbel_mcmc_",test_case,"_tree_",n.tree,"_stat", ".Rdata")))
  
  gumbel_adapt <- get(load(paste0("gumbel_mcmc_",test_case,"_tree_",n.tree,"_plot_adapt", ".Rdata")))
  gumbel_stat_adapt <- get(load(paste0("gumbel_mcmc_",test_case,"_tree_",n.tree,"_stat_adapt", ".Rdata")))
  
  gumbel_like_true <- loglik_gumbel((1 / (1-get(paste0("tau_true_",test_case)))), get(paste0("copula_uu_gumbel_",test_case))[,1] , get(paste0("copula_uu_gumbel_",test_case))[,2])
  
  
  frank <- get(load(paste0("frank_mcmc_",test_case,"_tree_",n.tree,"_plot", ".Rdata")))
  frank_stat <- get(load(paste0("frank_mcmc_",test_case,"_tree_",n.tree,"_stat", ".Rdata")))
  
  frank_adapt <- get(load(paste0("frank_mcmc_",test_case,"_tree_",n.tree,"_plot_adapt", ".Rdata")))
  frank_stat_adapt <- get(load(paste0("frank_mcmc_",test_case,"_tree_",n.tree,"_stat_adapt", ".Rdata")))
  
  
  frank_like_true <- loglik_frank(BiCopTau2Par(5,get(paste0("tau_true_",test_case))), get(paste0("copula_uu_frank_",test_case))[,1] , get(paste0("copula_uu_frank_",test_case))[,2])
  
  c(gauss_like_true, t_like_true, clayton_like_true, gumbel_like_true, frank_like_true)
  
  summ_stat_acc <- rbind(c(sqrt(gauss_stat$pred[1]),as.vector(gauss_stat$pred[-1]),sqrt(gauss_stat_adapt$pred[1]),as.vector(gauss_stat_adapt$pred[-1])),
                         c(sqrt(t_stat$pred[1]),as.vector(t_stat$pred[-1]),sqrt(t_stat_adapt$pred[1]),as.vector(t_stat_adapt$pred[-1])),
                         c(sqrt(clayton_stat$pred[1]),as.vector(clayton_stat$pred[-1]),sqrt(clayton_stat_adapt$pred[1]),as.vector(clayton_stat_adapt$pred[-1])),
                         c(sqrt(gumbel_stat$pred[1]),as.vector(gumbel_stat$pred[-1]),sqrt(gumbel_stat_adapt$pred[1]),as.vector(gumbel_stat_adapt$pred[-1])),
                         c(sqrt(frank_stat$pred[1]),as.vector(frank_stat$pred[-1]),sqrt(frank_stat_adapt$pred[1]),as.vector(frank_stat_adapt$pred[-1])))
  
  xtable(summ_stat_acc, digits = 3)
  
  if(n.tree ==1 || test_case == 1){
    summ_stat_tree <- rbind(c(as.vector(gauss_stat$tree),as.vector(gauss_stat_adapt$tree)),
                            c(as.vector(t_stat$tree),as.vector(t_stat_adapt$tree)),
                            c(as.vector(clayton_stat$tree),as.vector(clayton_stat_adapt$tree)),
                            c(as.vector(gumbel_stat$tree),as.vector(gumbel_stat_adapt$tree)),
                            c(as.vector(frank_stat$tree),as.vector(frank_stat_adapt$tree)))
    
    xtable(summ_stat_tree, digits = 2)
    
  }
  
  # 6, 7 potrait
  
  (gauss$like + labs(title="Gaussian (without adaption)") + geom_hline(yintercept = gauss_like_true, linetype = 2) + gauss_adapt$like + labs(title="Gaussian (with adaption)") + geom_hline(yintercept = gauss_like_true, linetype = 2)) / 
    (t$like + labs(title="Student-t (without adaption)") + geom_hline(yintercept = t_like_true, linetype = 2) + t_adapt$like + labs(title="Student-t (with adaption)") + geom_hline(yintercept = t_like_true, linetype = 2)) /
    (clayton$like + labs(title="Clayton (without adaption)")+ geom_hline(yintercept = clayton_like_true, linetype = 2) + clayton_adapt$like + labs(title="Clayton (with adaption)")+ geom_hline(yintercept = clayton_like_true, linetype = 2)) /
    (gumbel$like + labs(title="Gumbel (without adaption)")+ geom_hline(yintercept = gumbel_like_true, linetype = 2) + gumbel_adapt$like + labs(title="Gumbel (with adaption)")+ geom_hline(yintercept = gumbel_like_true, linetype = 2)) /
    (frank$like + labs(title="Frank (without adaption)")+ geom_hline(yintercept = frank_like_true, linetype = 2) + frank_adapt$like + labs(title="Frank (with adaption)")+ geom_hline(yintercept = frank_like_true, linetype = 2))
  
  (gauss$pred + labs(title="Gaussian (without adaption)") + ylim(0,1) + gauss_adapt$pred + labs(title="Gaussian (with adaption)") + ylim(0,1)) / 
    (t$pred + labs(title="Student-t (without adaption)") + ylim(0,1) + t_adapt$pred + labs(title="Student-t (with adaption)") + ylim(0,1)) /
    (clayton$pred + labs(title="Clayton (without adaption)") + ylim(0,1) + clayton_adapt$pred + labs(title="Clayton (with adaption)") + ylim(0,1)) /
    (gumbel$pred + labs(title="Gumbel (without adaption)") + ylim(0,1) + gumbel_adapt$pred + labs(title="Gumbel (with adaption)") + ylim(0,1)) /
    (frank$pred + labs(title="Frank (without adaption)") + ylim(0,1) + frank_adapt$pred + labs(title="Frank (with adaption)") + ylim(0,1))
  
  if(n.tree == 1 || test_case == 1){
    (gauss$nterm + labs(title="Gaussian (without adaption)") + geom_hline(yintercept = 3, linetype = 2) + ylim(0, 8) + gauss_adapt$nterm + labs(title="Gaussian (with adaption)") + geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8)) /
      (t$nterm + labs(title="Student-t (without adaption)") + geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8) + t_adapt$nterm + labs(title="Student-t (with adaption)") + geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8)) /
      (clayton$nterm + labs(title="Clayton (without adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8) + clayton_adapt$nterm + labs(title="Clayton (with adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8)) /
      (gumbel$nterm + labs(title="Gumbel (without adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8) + gumbel_adapt$nterm + labs(title="Gumbel (with adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8)) /
      (frank$nterm + labs(title="Frank (without adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8) + frank_adapt$nterm + labs(title="Frank (with adaption)")+ geom_hline(yintercept = 3, linetype = 2)+ ylim(0, 8))
    
    (gauss$depth + labs(title="Gaussian (without adaption)") + geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + gauss_adapt$depth + labs(title="Gaussian (with adaption)") + geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6)) /
      (t$depth + labs(title="Student-t (without adaption)") + geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + t_adapt$depth + labs(title="Student-t (with adaption)") + geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6)) /
      (clayton$depth + labs(title="Clayton (without adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + clayton_adapt$depth + labs(title="Clayton (with adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6)) /
      (gumbel$depth + labs(title="Gumbel (without adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + gumbel_adapt$depth + labs(title="Gumbel (with adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6)) /
      (frank$depth + labs(title="Frank (without adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6) + frank_adapt$depth + labs(title="Frank (with adaption)")+ geom_hline(yintercept = 2, linetype = 2)+ ylim(0, 6))
  }
  
}