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
  tau_true_1[(X_obs>=0.33)&(X_obs<0.66)] <- 0.9
  tau_true_1[(X_obs>=0.66)] <- 0.3
  
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

data_pred <- function(pred_list, x_list, family = "gauss"){
  
  pred_cond_mod <- list()
  if(family == "gauss"){
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
  } else if(family == "t"){
    for (i in c(1:length(x_list))) {
      pred_cond <- data.frame("obs" = rep(x_list[[i]], each = (n.chain_par * n.iter_par)))
      pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))[[i]]), each = (n.chain_par * (n.iter_par)))
      pred_cond$y = BiCopPar2Tau(2,link_t(as.vector(pred_list[[i]])))
      pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
      pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
      pred_cond$set = rep(i, (n.iter_par*n.chain_par*n))
      
      pred_cond_mod[[i]] = pred_cond %>%
        filter(idx > n.discard) %>%
        group_by(obs, chain, theta_true) %>%
        summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
      
    }
  } else if(family == "clayton"){
    for (i in c(1:length(x_list))) {
      pred_cond <- data.frame("obs" = rep(x_list[[i]], each = (n.chain_par * n.iter_par)))
      pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))[[i]]), each = (n.chain_par * (n.iter_par)))
      pred_cond$y = BiCopPar2Tau(3,link_clayton(as.vector(pred_list[[i]])))
      pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
      pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
      pred_cond$set = rep(i, (n.iter_par*n.chain_par*n))
      
      pred_cond_mod[[i]] = pred_cond %>%
        filter(idx > n.discard) %>%
        group_by(obs, chain, theta_true) %>%
        summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
      
    }
  } else if(family == "gumbel"){
    for (i in c(1:length(x_list))) {
      pred_cond <- data.frame("obs" = rep(x_list[[i]], each = (n.chain_par * n.iter_par)))
      pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))[[i]]), each = (n.chain_par * (n.iter_par)))
      pred_cond$y = BiCopPar2Tau(4,link_gumbel(as.vector(pred_list[[i]])))
      pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
      pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
      pred_cond$set = rep(i, (n.iter_par*n.chain_par*n))
      
      pred_cond_mod[[i]] = pred_cond %>%
        filter(idx > n.discard) %>%
        group_by(obs, chain, theta_true) %>%
        summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
      
    }
  } else if(family == "frank"){
    for (i in c(1:length(x_list))) {
      pred_cond <- data.frame("obs" = rep(x_list[[i]], each = (n.chain_par * n.iter_par)))
      pred_cond$theta_true = rep((get(paste0("tau_true_",test_case))[[i]]), each = (n.chain_par * (n.iter_par)))
      pred_cond$y = BiCopPar2Tau(5,link_frank(as.vector(pred_list[[i]])))
      pred_cond$chain = rep(rep(1:n.chain_par, each = n.iter_par),n)
      pred_cond$idx = rep(rep(1:n.iter_par, n.chain_par),n)
      pred_cond$set = rep(i, (n.iter_par*n.chain_par*n))
      
      pred_cond_mod[[i]] = pred_cond %>%
        filter(idx > n.discard) %>%
        group_by(obs, chain, theta_true) %>%
        summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
      
    }
  }
  return(pred_cond_mod)
}


# gauss

if(F){
  
  load(paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  
  pred_cond_mod_list <- data_pred(gauss_dat_list$all$pred, X_obs.norm, family = "gauss")
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  pred_cond_mod$set <- rep(1:R, each = n*n.chain_par)
  
  pred_cond_mod = pred_cond_mod %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025)))
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(set) %>%
    summarise(RMSE = mean(RMSE), mean_CI.length = mean(CI.length), CI.cov = mean(CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_mod_avg[,-1])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    acc_lb.df <- do.call(rbind,gauss_dat_list$all$acc)
    acc_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_acc <- acc_lb.df %>%
      group_by(set) %>%
      summarise(acc = mean(nn==TRUE))
    
    nt_lb.df <- do.call(rbind, gauss_dat_list$all$nterm)
    nt_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_nt <- nt_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(nterm = mean(nn))
    
    depth_lb.df <- do.call(rbind, gauss_dat_list$all$depth)
    depth_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_depth <- depth_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(depth = mean(nn))
    
    tree_stat_df = cbind(mean_nt$nterm, mean_depth$depth, mean_acc$acc)
    
    tree_stat_summary = c(colMeans(tree_stat_df), apply(tree_stat_df,2,sd))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc", "sd.nterm", "sd.depth", "sd.acc")
  }
  
  
  gauss_dat_list <- list("all" = gauss_dat_list$all, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    gauss_stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    # gauss_plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    gauss_stat_list <- list("pred" = pred_cond_summary)
    # gauss_plot_list <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(gauss_dat_list, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  save(gauss_stat_list, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  # save(gauss_plot_list, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, gauss_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, gauss_dat_list, gauss_stat_list, gauss_plot_list,
   mean_acc,mean_depth,mean_nt,tree_stat_df)
gc()

if(F){
  
  load(paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  
  pred_cond_mod_list <- data_pred(gauss_dat_list_adapt$all$pred, X_obs.norm, family = "gauss")
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  pred_cond_mod$set <- rep(1:R, each = n*n.chain_par)
  
  pred_cond_mod = pred_cond_mod %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025)))
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(set) %>%
    summarise(RMSE = mean(RMSE), mean_CI.length = mean(CI.length), CI.cov = mean(CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_mod_avg[,-1])
  
  if((n.tree==1)||(test_case == 1)){
    
    acc_lb.df <- do.call(rbind,gauss_dat_list_adapt$all$acc)
    acc_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_acc <- acc_lb.df %>%
      group_by(set) %>%
      summarise(acc = mean(nn==TRUE))
    
    nt_lb.df <- do.call(rbind, gauss_dat_list_adapt$all$nterm)
    nt_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_nt <- nt_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(nterm = mean(nn))
    
    depth_lb.df <- do.call(rbind, gauss_dat_list_adapt$all$depth)
    depth_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_depth <- depth_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(depth = mean(nn))
    
    tree_stat_df = cbind(mean_nt$nterm, mean_depth$depth, mean_acc$acc)
    
    tree_stat_summary = c(colMeans(tree_stat_df), apply(tree_stat_df,2,sd))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc", "sd.nterm", "sd.depth", "sd.acc")
  }
  
  
  gauss_dat_list_adapt <- list("all" = gauss_dat_list_adapt$all, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    gauss_stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    # gauss_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    gauss_stat_list_adapt <- list("pred" = pred_cond_summary)
    # gauss_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(gauss_dat_list_adapt, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  save(gauss_stat_list_adapt, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  # save(gauss_plot_list_adapt, file = paste0("gauss_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, gauss_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, gauss_dat_list_adapt, gauss_stat_list_adapt, gauss_plot_list_adapt,
   mean_acc,mean_depth,mean_nt,tree_stat_df)
gc()

# t

if(T){
  
  load(paste0("t_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  
  pred_cond_mod_list <- data_pred(t_dat_list$all$pred, X_obs.norm, family = "t")
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  pred_cond_mod$set <- rep(1:R, each = n*n.chain_par)
  
  pred_cond_mod = pred_cond_mod %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025)))
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(set) %>%
    summarise(RMSE = mean(RMSE), mean_CI.length = mean(CI.length), CI.cov = mean(CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_mod_avg[,-1])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    acc_lb.df <- do.call(rbind,t_dat_list$all$acc)
    acc_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_acc <- acc_lb.df %>%
      group_by(set) %>%
      summarise(acc = mean(nn==TRUE))
    
    nt_lb.df <- do.call(rbind, t_dat_list$all$nterm)
    nt_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_nt <- nt_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(nterm = mean(nn))
    
    depth_lb.df <- do.call(rbind, t_dat_list$all$depth)
    depth_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_depth <- depth_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(depth = mean(nn))
    
    tree_stat_df = cbind(mean_nt$nterm, mean_depth$depth, mean_acc$acc)
    
    tree_stat_summary = c(colMeans(tree_stat_df), apply(tree_stat_df,2,sd))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc", "sd.nterm", "sd.depth", "sd.acc")
  }
  
  
  t_dat_list <- list("all" = t_dat_list$all, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    t_stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    # t_plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    t_stat_list <- list("pred" = pred_cond_summary)
    # t_plot_list <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(t_dat_list, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  save(t_stat_list, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  # save(t_plot_list, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, t_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, t_dat_list, t_stat_list, t_plot_list,
   mean_acc,mean_depth,mean_nt,tree_stat_df)
gc()

if(T){
  
  load(paste0("t_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  
  pred_cond_mod_list <- data_pred(t_dat_list_adapt$all$pred, X_obs.norm, family = "t")
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  pred_cond_mod$set <- rep(1:R, each = n*n.chain_par)
  
  pred_cond_mod = pred_cond_mod %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025)))
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(set) %>%
    summarise(RMSE = mean(RMSE), mean_CI.length = mean(CI.length), CI.cov = mean(CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_mod_avg[,-1])
  
  if((n.tree==1)||(test_case == 1)){
    
    acc_lb.df <- do.call(rbind,t_dat_list_adapt$all$acc)
    acc_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_acc <- acc_lb.df %>%
      group_by(set) %>%
      summarise(acc = mean(nn==TRUE))
    
    nt_lb.df <- do.call(rbind, t_dat_list_adapt$all$nterm)
    nt_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_nt <- nt_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(nterm = mean(nn))
    
    depth_lb.df <- do.call(rbind, t_dat_list_adapt$all$depth)
    depth_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_depth <- depth_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(depth = mean(nn))
    
    tree_stat_df = cbind(mean_nt$nterm, mean_depth$depth, mean_acc$acc)
    
    tree_stat_summary = c(colMeans(tree_stat_df), apply(tree_stat_df,2,sd))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc", "sd.nterm", "sd.depth", "sd.acc")
  }
  
  
  t_dat_list_adapt <- list("all" = t_dat_list_adapt$all, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    t_stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    # t_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    t_stat_list_adapt <- list("pred" = pred_cond_summary)
    # t_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(t_dat_list_adapt, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  save(t_stat_list_adapt, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  # save(t_plot_list_adapt, file = paste0("t_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, t_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, t_dat_list_adapt, t_stat_list_adapt, t_plot_list_adapt,
   mean_acc,mean_depth,mean_nt,tree_stat_df)
gc()

# clayton

if(T){
  
  load(paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  
  pred_cond_mod_list <- data_pred(clayton_dat_list$all$pred, X_obs.norm, family = "clayton")
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  pred_cond_mod$set <- rep(1:R, each = n*n.chain_par)
  
  pred_cond_mod = pred_cond_mod %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025)))
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(set) %>%
    summarise(RMSE = mean(RMSE), mean_CI.length = mean(CI.length), CI.cov = mean(CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_mod_avg[,-1])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    acc_lb.df <- do.call(rbind,clayton_dat_list$all$acc)
    acc_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_acc <- acc_lb.df %>%
      group_by(set) %>%
      summarise(acc = mean(nn==TRUE))
    
    nt_lb.df <- do.call(rbind, clayton_dat_list$all$nterm)
    nt_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_nt <- nt_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(nterm = mean(nn))
    
    depth_lb.df <- do.call(rbind, clayton_dat_list$all$depth)
    depth_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_depth <- depth_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(depth = mean(nn))
    
    tree_stat_df = cbind(mean_nt$nterm, mean_depth$depth, mean_acc$acc)
    
    tree_stat_summary = c(colMeans(tree_stat_df), apply(tree_stat_df,2,sd))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc", "sd.nterm", "sd.depth", "sd.acc")
  }
  
  
  clayton_dat_list <- list("all" = clayton_dat_list$all, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    clayton_stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    # clayton_plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    clayton_stat_list <- list("pred" = pred_cond_summary)
    # clayton_plot_list <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(clayton_dat_list, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  save(clayton_stat_list, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  # save(clayton_plot_list, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, clayton_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, clayton_dat_list, clayton_stat_list, clayton_plot_list,
   mean_acc,mean_depth,mean_nt,tree_stat_df)
gc()

if(T){
  
  load(paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  
  pred_cond_mod_list <- data_pred(clayton_dat_list_adapt$all$pred, X_obs.norm, family = "clayton")
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  pred_cond_mod$set <- rep(1:R, each = n*n.chain_par)
  
  pred_cond_mod = pred_cond_mod %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025)))
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(set) %>%
    summarise(RMSE = mean(RMSE), mean_CI.length = mean(CI.length), CI.cov = mean(CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_mod_avg[,-1])
  
  if((n.tree==1)||(test_case == 1)){
    
    acc_lb.df <- do.call(rbind,clayton_dat_list_adapt$all$acc)
    acc_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_acc <- acc_lb.df %>%
      group_by(set) %>%
      summarise(acc = mean(nn==TRUE))
    
    nt_lb.df <- do.call(rbind, clayton_dat_list_adapt$all$nterm)
    nt_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_nt <- nt_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(nterm = mean(nn))
    
    depth_lb.df <- do.call(rbind, clayton_dat_list_adapt$all$depth)
    depth_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_depth <- depth_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(depth = mean(nn))
    
    tree_stat_df = cbind(mean_nt$nterm, mean_depth$depth, mean_acc$acc)
    
    tree_stat_summary = c(colMeans(tree_stat_df), apply(tree_stat_df,2,sd))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc", "sd.nterm", "sd.depth", "sd.acc")
  }
  
  
  clayton_dat_list_adapt <- list("all" = clayton_dat_list_adapt$all, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    clayton_stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    # clayton_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    clayton_stat_list_adapt <- list("pred" = pred_cond_summary)
    # clayton_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(clayton_dat_list_adapt, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  save(clayton_stat_list_adapt, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  # save(clayton_plot_list_adapt, file = paste0("clayton_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, clayton_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, clayton_dat_list_adapt, clayton_stat_list_adapt, clayton_plot_list_adapt,
   mean_acc,mean_depth,mean_nt,tree_stat_df)
gc()

# gumbel

if(T){
  
  load(paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  
  pred_cond_mod_list <- data_pred(gumbel_dat_list$all$pred, X_obs.norm, family = "gumbel")
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  pred_cond_mod$set <- rep(1:R, each = n*n.chain_par)
  
  pred_cond_mod = pred_cond_mod %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025)))
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(set) %>%
    summarise(RMSE = mean(RMSE), mean_CI.length = mean(CI.length), CI.cov = mean(CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_mod_avg[,-1])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    acc_lb.df <- do.call(rbind,gumbel_dat_list$all$acc)
    acc_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_acc <- acc_lb.df %>%
      group_by(set) %>%
      summarise(acc = mean(nn==TRUE))
    
    nt_lb.df <- do.call(rbind, gumbel_dat_list$all$nterm)
    nt_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_nt <- nt_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(nterm = mean(nn))
    
    depth_lb.df <- do.call(rbind, gumbel_dat_list$all$depth)
    depth_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_depth <- depth_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(depth = mean(nn))
    
    tree_stat_df = cbind(mean_nt$nterm, mean_depth$depth, mean_acc$acc)
    
    tree_stat_summary = c(colMeans(tree_stat_df), apply(tree_stat_df,2,sd))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc", "sd.nterm", "sd.depth", "sd.acc")
  }
  
  
  gumbel_dat_list <- list("all" = gumbel_dat_list$all, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    gumbel_stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    # gumbel_plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    gumbel_stat_list <- list("pred" = pred_cond_summary)
    # gumbel_plot_list <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(gumbel_dat_list, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  save(gumbel_stat_list, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  # save(gumbel_plot_list, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, gumbel_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, gumbel_dat_list, gumbel_stat_list, gumbel_plot_list,
   mean_acc,mean_depth,mean_nt,tree_stat_df)
gc()

if(T){
  
  load(paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  
  pred_cond_mod_list <- data_pred(gumbel_dat_list_adapt$all$pred, X_obs.norm, family = "gumbel")
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  pred_cond_mod$set <- rep(1:R, each = n*n.chain_par)
  
  pred_cond_mod = pred_cond_mod %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025)))
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(set) %>%
    summarise(RMSE = mean(RMSE), mean_CI.length = mean(CI.length), CI.cov = mean(CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_mod_avg[,-1])
  
  if((n.tree==1)||(test_case == 1)){
    
    acc_lb.df <- do.call(rbind,gumbel_dat_list_adapt$all$acc)
    acc_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_acc <- acc_lb.df %>%
      group_by(set) %>%
      summarise(acc = mean(nn==TRUE))
    
    nt_lb.df <- do.call(rbind, gumbel_dat_list_adapt$all$nterm)
    nt_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_nt <- nt_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(nterm = mean(nn))
    
    depth_lb.df <- do.call(rbind, gumbel_dat_list_adapt$all$depth)
    depth_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_depth <- depth_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(depth = mean(nn))
    
    tree_stat_df = cbind(mean_nt$nterm, mean_depth$depth, mean_acc$acc)
    
    tree_stat_summary = c(colMeans(tree_stat_df), apply(tree_stat_df,2,sd))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc", "sd.nterm", "sd.depth", "sd.acc")
  }
  
  
  gumbel_dat_list_adapt <- list("all" = gumbel_dat_list_adapt$all, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    gumbel_stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    # gumbel_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    gumbel_stat_list_adapt <- list("pred" = pred_cond_summary)
    # gumbel_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(gumbel_dat_list_adapt, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  save(gumbel_stat_list_adapt, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  # save(gumbel_plot_list_adapt, file = paste0("gumbel_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, gumbel_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, gumbel_dat_list_adapt, gumbel_stat_list_adapt, gumbel_plot_list_adapt,
   mean_acc,mean_depth,mean_nt,tree_stat_df)
gc()

# frank

if(T){
  
  load(paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  
  pred_cond_mod_list <- data_pred(frank_dat_list$all$pred, X_obs.norm, family = "frank")
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  pred_cond_mod$set <- rep(1:R, each = n*n.chain_par)
  
  pred_cond_mod = pred_cond_mod %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025)))
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(set) %>%
    summarise(RMSE = mean(RMSE), mean_CI.length = mean(CI.length), CI.cov = mean(CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_mod_avg[,-1])
  
  if((n.tree==1)||(test_case == 1)){
    # nterm
    
    acc_lb.df <- do.call(rbind,frank_dat_list$all$acc)
    acc_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_acc <- acc_lb.df %>%
      group_by(set) %>%
      summarise(acc = mean(nn==TRUE))
    
    nt_lb.df <- do.call(rbind, frank_dat_list$all$nterm)
    nt_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_nt <- nt_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(nterm = mean(nn))
    
    depth_lb.df <- do.call(rbind, frank_dat_list$all$depth)
    depth_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_depth <- depth_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(depth = mean(nn))
    
    tree_stat_df = cbind(mean_nt$nterm, mean_depth$depth, mean_acc$acc)
    
    tree_stat_summary = c(colMeans(tree_stat_df), apply(tree_stat_df,2,sd))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc", "sd.nterm", "sd.depth", "sd.acc")
  }
  
  
  frank_dat_list <- list("all" = frank_dat_list$all, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    frank_stat_list <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    # frank_plot_list <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    frank_stat_list <- list("pred" = pred_cond_summary)
    # frank_plot_list <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(frank_dat_list, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_dat.Rdata"))
  save(frank_stat_list, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_stat.Rdata"))
  # save(frank_plot_list, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_plot.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, frank_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, frank_dat_list, frank_stat_list, frank_plot_list,
   mean_acc,mean_depth,mean_nt,tree_stat_df)
gc()

if(T){
  
  load(paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  
  pred_cond_mod_list <- data_pred(frank_dat_list_adapt$all$pred, X_obs.norm, family = "frank")
  
  pred_cond_mod <- do.call(rbind, pred_cond_mod_list)
  pred_cond_mod$set <- rep(1:R, each = n*n.chain_par)
  
  pred_cond_mod = pred_cond_mod %>%
    mutate(RMSE = ((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = (theta_q975 - theta_q025)) %>%
    mutate(CI.cov = ((theta_true < theta_q975) & (theta_true > theta_q025)))
  
  pred_cond_mod_avg = pred_cond_mod %>%
    group_by(set) %>%
    summarise(RMSE = mean(RMSE), mean_CI.length = mean(CI.length), CI.cov = mean(CI.cov)) 
  
  pred_cond_summary = colMeans(pred_cond_mod_avg[,-1])
  
  if((n.tree==1)||(test_case == 1)){
    
    acc_lb.df <- do.call(rbind,frank_dat_list_adapt$all$acc)
    acc_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_acc <- acc_lb.df %>%
      group_by(set) %>%
      summarise(acc = mean(nn==TRUE))
    
    nt_lb.df <- do.call(rbind, frank_dat_list_adapt$all$nterm)
    nt_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_nt <- nt_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(nterm = mean(nn))
    
    depth_lb.df <- do.call(rbind, frank_dat_list_adapt$all$depth)
    depth_lb.df$set <- rep(1:R, each = n.chain_par * n.iter_par)
    
    mean_depth <- depth_lb.df %>%
      filter(idx>n.discard) %>%
      group_by(set) %>%
      summarise(depth = mean(nn))
    
    tree_stat_df = cbind(mean_nt$nterm, mean_depth$depth, mean_acc$acc)
    
    tree_stat_summary = c(colMeans(tree_stat_df), apply(tree_stat_df,2,sd))
    
    names(tree_stat_summary) <- c("nterm", "depth", "acc", "sd.nterm", "sd.depth", "sd.acc")
  }
  
  
  frank_dat_list_adapt <- list("all" = frank_dat_list_adapt$all, "pred" = pred_cond_mod_list)
  
  if((n.tree==1)||(test_case == 1)){
    frank_stat_list_adapt <- list("pred" = pred_cond_summary, "tree" = tree_stat_summary)
    # frank_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred, "depth" = pl_dp, "nterm" = pl_nl)
  } else{
    frank_stat_list_adapt <- list("pred" = pred_cond_summary)
    # frank_plot_list_adapt <- list("like" = pl_like, "pred" = pl_pred)
  }
  
  save(frank_dat_list_adapt, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_dat_adapt.Rdata"))
  save(frank_stat_list_adapt, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_stat_adapt.Rdata"))
  # save(frank_plot_list_adapt, file = paste0("frank_mcmc_",test_case,"_tree_",n.tree, "_plot_adapt.Rdata"))
  
}

rm(summ_list, pred_cond_mod_list, pred_cond_case_1, pred_cond_case_1_mod, pred_cond_case_1_mod_avg, pred_cond_mod, pred_cond_mod_avg, pred_cond_stat,
   pred_cond_summary, like_df_case_1, pl_like, pl_pred, frank_like_true, nt_lb.df_case_1, pl_nl, depth_lb.df_case_1, pl_dp,
   acc_lb.df, nt_lb.df, depth_lb.df, tree_stat_summary, frank_dat_list_adapt, frank_stat_list_adapt, frank_plot_list_adapt,
   mean_acc,mean_depth,mean_nt,tree_stat_df)
gc()
