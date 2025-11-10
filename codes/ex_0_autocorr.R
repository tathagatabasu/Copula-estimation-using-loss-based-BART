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
  tau_true_1[(X_obs>=0.33)&(X_obs<0.66)] <- 0.8
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
################################################################################

load("mcmc_1/gauss_mcmc_1_tree_1_dat_adapt.Rdata")

load("mcmc_1/gauss_mcmc_1_tree_1_dat.Rdata")

like_df <- list()

for (i in 1:R) {
  like_df[[i]] <-data.frame("nn" = apply(gauss_dat_list$all$pred[[i]], 1, function(x)loglik_gauss(link_gauss(x), get(paste0("copula_uu_gauss_",test_case))[[i]][,1], get(paste0("copula_uu_gauss_",test_case))[[i]][,2])))
  like_df[[i]]$idx <- rep(1:n.iter_par, n.chain_par)
  like_df[[i]]$chain <- rep(1:n.chain_par, each = n.iter_par)
}

auto_corr_gauss <- sapply(1:R, function(i)autocorr(mcmc(like_df[[i]]$nn[like_df[[i]]$idx > 1500]), lags = seq(0,400,10)))

like_df_adapt <- list()

for (i in 1:R) {
  like_df_adapt[[i]] <-data.frame("nn" = apply(gauss_dat_list_adapt$all$pred[[i]], 1, function(x)loglik_gauss(link_gauss(x), get(paste0("copula_uu_gauss_",test_case))[[i]][,1], get(paste0("copula_uu_gauss_",test_case))[[i]][,2])))
  like_df_adapt[[i]]$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_adapt[[i]]$chain <- rep(1:n.chain_par, each = n.iter_par)
}

auto_corr_gauss_adapt <- sapply(1:R, function(i)autocorr(mcmc(like_df_adapt[[i]]$nn[like_df_adapt[[i]]$idx > 1500]), lags = seq(0,400,10)))

rm(gauss_dat_list, gauss_dat_list_adapt)

load("mcmc_1/t_mcmc_1_tree_1_dat_adapt.Rdata")

load("mcmc_1/t_mcmc_1_tree_1_dat.Rdata")

like_df <- list()

for (i in 1:R) {
  like_df[[i]] <-data.frame("nn" = apply(t_dat_list$all$pred[[i]], 1, function(x)loglik_t(link_t(x), get(paste0("copula_uu_t_",test_case))[[i]][,1], get(paste0("copula_uu_t_",test_case))[[i]][,2])))
  like_df[[i]]$idx <- rep(1:n.iter_par, n.chain_par)
  like_df[[i]]$chain <- rep(1:n.chain_par, each = n.iter_par)
  print(i)
}

auto_corr_t <- sapply(1:R, function(i)autocorr(mcmc(like_df[[i]]$nn[like_df[[i]]$idx > 1500]), lags = seq(0,400,10)))

like_df_adapt <- list()

for (i in 1:R) {
  like_df_adapt[[i]] <-data.frame("nn" = apply(t_dat_list_adapt$all$pred[[i]], 1, function(x)loglik_t(link_t(x), get(paste0("copula_uu_t_",test_case))[[i]][,1], get(paste0("copula_uu_t_",test_case))[[i]][,2])))
  like_df_adapt[[i]]$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_adapt[[i]]$chain <- rep(1:n.chain_par, each = n.iter_par)
  print(i)
}

auto_corr_t_adapt <- sapply(1:R, function(i)autocorr(mcmc(like_df_adapt[[i]]$nn[like_df_adapt[[i]]$idx > 1500]), lags = seq(0,400,10)))

rm(t_dat_list, t_dat_list_adapt)

load("mcmc_1/clayton_mcmc_1_tree_1_dat_adapt.Rdata")

load("mcmc_1/clayton_mcmc_1_tree_1_dat.Rdata")

like_df <- list()

for (i in 1:R) {
  like_df[[i]] <-data.frame("nn" = apply(clayton_dat_list$all$pred[[i]], 1, function(x)loglik_clayton(link_clayton(x), get(paste0("copula_uu_clayton_",test_case))[[i]][,1], get(paste0("copula_uu_clayton_",test_case))[[i]][,2])))
  like_df[[i]]$idx <- rep(1:n.iter_par, n.chain_par)
  like_df[[i]]$chain <- rep(1:n.chain_par, each = n.iter_par)
  print(i)
}

auto_corr_clayton <- sapply(1:R, function(i)autocorr(mcmc(like_df[[i]]$nn[like_df[[i]]$idx > 1500]), lags = seq(0,400,10)))

like_df_adapt <- list()

for (i in 1:R) {
  like_df_adapt[[i]] <-data.frame("nn" = apply(clayton_dat_list_adapt$all$pred[[i]], 1, function(x)loglik_clayton(link_clayton(x), get(paste0("copula_uu_clayton_",test_case))[[i]][,1], get(paste0("copula_uu_clayton_",test_case))[[i]][,2])))
  like_df_adapt[[i]]$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_adapt[[i]]$chain <- rep(1:n.chain_par, each = n.iter_par)
  print(i)
}

auto_corr_clayton_adapt <- sapply(1:R, function(i)autocorr(mcmc(like_df_adapt[[i]]$nn[like_df_adapt[[i]]$idx > 1500]), lags = seq(0,400,10)))

rm(clayton_dat_list, clayton_dat_list_adapt)

load("mcmc_1/gumbel_mcmc_1_tree_1_dat_adapt.Rdata")

load("mcmc_1/gumbel_mcmc_1_tree_1_dat.Rdata")

like_df <- list()

for (i in 1:R) {
  like_df[[i]] <-data.frame("nn" = apply(gumbel_dat_list$all$pred[[i]], 1, function(x)loglik_gumbel(link_gumbel(x), get(paste0("copula_uu_gumbel_",test_case))[[i]][,1], get(paste0("copula_uu_gumbel_",test_case))[[i]][,2])))
  like_df[[i]]$idx <- rep(1:n.iter_par, n.chain_par)
  like_df[[i]]$chain <- rep(1:n.chain_par, each = n.iter_par)
  print(i)
}

auto_corr_gumbel <- sapply(1:R, function(i)autocorr(mcmc(like_df[[i]]$nn[like_df[[i]]$idx > 1500]), lags = seq(0,400,10)))

like_df_adapt <- list()

for (i in 1:R) {
  like_df_adapt[[i]] <-data.frame("nn" = apply(gumbel_dat_list_adapt$all$pred[[i]], 1, function(x)loglik_gumbel(link_gumbel(x), get(paste0("copula_uu_gumbel_",test_case))[[i]][,1], get(paste0("copula_uu_gumbel_",test_case))[[i]][,2])))
  like_df_adapt[[i]]$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_adapt[[i]]$chain <- rep(1:n.chain_par, each = n.iter_par)
  print(i)
}

auto_corr_gumbel_adapt <- sapply(1:R, function(i)autocorr(mcmc(like_df_adapt[[i]]$nn[like_df_adapt[[i]]$idx > 1500]), lags = seq(0,400,10)))

rm(gumbel_dat_list, gumbel_dat_list_adapt)

load("mcmc_1/frank_mcmc_1_tree_1_dat_adapt_new_prop.Rdata")

load("mcmc_1/frank_mcmc_1_tree_1_dat_new_prop.Rdata")

like_df <- list()

for (i in 1:R) {
  like_df[[i]] <-data.frame("nn" = apply(frank_dat_list$all$pred[[i]], 1, function(x)loglik_frank(link_frank(x), get(paste0("copula_uu_frank_",test_case))[[i]][,1], get(paste0("copula_uu_frank_",test_case))[[i]][,2])))
  like_df[[i]]$idx <- rep(1:n.iter_par, n.chain_par)
  like_df[[i]]$chain <- rep(1:n.chain_par, each = n.iter_par)
  print(i)
}

auto_corr_frank <- sapply(1:R, function(i)autocorr(mcmc(like_df[[i]]$nn[like_df[[i]]$idx > 1500]), lags = seq(0,400,10)))

like_df_adapt <- list()

for (i in 1:R) {
  like_df_adapt[[i]] <-data.frame("nn" = apply(frank_dat_list_adapt$all$pred[[i]], 1, function(x)loglik_frank(link_frank(x), get(paste0("copula_uu_frank_",test_case))[[i]][,1], get(paste0("copula_uu_frank_",test_case))[[i]][,2])))
  like_df_adapt[[i]]$idx <- rep(1:n.iter_par, n.chain_par)
  like_df_adapt[[i]]$chain <- rep(1:n.chain_par, each = n.iter_par)
  print(i)
}

auto_corr_frank_adapt <- sapply(1:R, function(i)autocorr(mcmc(like_df_adapt[[i]]$nn[like_df_adapt[[i]]$idx > 1500]), lags = seq(0,400,10)))

rm(frank_dat_list, frank_dat_list_adapt)

save(auto_corr_gauss, auto_corr_gauss_adapt, auto_corr_t, auto_corr_t_adapt, 
     auto_corr_clayton, auto_corr_clayton_adapt, auto_corr_gumbel, auto_corr_gumbel_adapt,
     auto_corr_frank, auto_corr_frank_adapt, file = "autocorr_dat_case_1.Rdata")


if(F){
  library(ggplot2)
  
  load("autocorr_dat_case_2.Rdata")
  
  auto_corr_clayton_med <- apply(auto_corr_clayton,1,median)
  auto_corr_clayton_q025 <- apply(auto_corr_clayton,1,\(x)quantile(x, .025))
  auto_corr_clayton_q975 <- apply(auto_corr_clayton,1,\(x)quantile(x, .975))
  
  auto_corr_clayton_adapt_med <- apply(auto_corr_clayton_adapt,1,median)
  auto_corr_clayton_adapt_q025 <- apply(auto_corr_clayton_adapt,1,\(x)quantile(x, .025))
  auto_corr_clayton_adapt_q975 <- apply(auto_corr_clayton_adapt,1,\(x)quantile(x, .975))
  
  ggplot() +
    geom_line(aes(x = seq(0,400,10), y = auto_corr_clayton_med), col = "red") +
    geom_line(aes(x = seq(0,400,10), y = auto_corr_clayton_adapt_med), col = "green")+
    geom_ribbon(aes(x = seq(0,400,10), ymin = auto_corr_clayton_q025, ymax = auto_corr_clayton_q975), fill = "red", alpha = 0.2)+
    geom_ribbon(aes(x = seq(0,400,10), ymin = auto_corr_clayton_adapt_q025, ymax = auto_corr_clayton_adapt_q975), fill = "green", alpha = 0.2)
  
  auto_corr_gumbel_med <- apply(auto_corr_gumbel,1,median)
  auto_corr_gumbel_q025 <- apply(auto_corr_gumbel,1,\(x)quantile(x, .025))
  auto_corr_gumbel_q975 <- apply(auto_corr_gumbel,1,\(x)quantile(x, .975))
  
  auto_corr_gumbel_adapt_med <- apply(auto_corr_gumbel_adapt,1,median)
  auto_corr_gumbel_adapt_q025 <- apply(auto_corr_gumbel_adapt,1,\(x)quantile(x, .025))
  auto_corr_gumbel_adapt_q975 <- apply(auto_corr_gumbel_adapt,1,\(x)quantile(x, .975))
  
  ggplot() +
    geom_line(aes(x = seq(0,400,10), y = auto_corr_gumbel_med), col = "red") +
    geom_line(aes(x = seq(0,400,10), y = auto_corr_gumbel_adapt_med), col = "green")+
    geom_ribbon(aes(x = seq(0,400,10), ymin = auto_corr_gumbel_q025, ymax = auto_corr_gumbel_q975), fill = "red", alpha = 0.2)+
    geom_ribbon(aes(x = seq(0,400,10), ymin = auto_corr_gumbel_adapt_q025, ymax = auto_corr_gumbel_adapt_q975), fill = "green", alpha = 0.2)
  
  auto_corr_gauss_med <- apply(auto_corr_gauss,1,median)
  auto_corr_gauss_q025 <- apply(auto_corr_gauss,1,\(x)quantile(x, .025))
  auto_corr_gauss_q975 <- apply(auto_corr_gauss,1,\(x)quantile(x, .975))
  
  auto_corr_gauss_adapt_med <- apply(auto_corr_gauss_adapt,1,median)
  auto_corr_gauss_adapt_q025 <- apply(auto_corr_gauss_adapt,1,\(x)quantile(x, .025))
  auto_corr_gauss_adapt_q975 <- apply(auto_corr_gauss_adapt,1,\(x)quantile(x, .975))
  
  ggplot() +
    geom_line(aes(x = seq(0,400,10), y = auto_corr_gauss_med), col = "red") +
    geom_line(aes(x = seq(0,400,10), y = auto_corr_gauss_adapt_med), col = "green")+
    geom_ribbon(aes(x = seq(0,400,10), ymin = auto_corr_gauss_q025, ymax = auto_corr_gauss_q975), fill = "red", alpha = 0.2)+
    geom_ribbon(aes(x = seq(0,400,10), ymin = auto_corr_gauss_adapt_q025, ymax = auto_corr_gauss_adapt_q975), fill = "green", alpha = 0.2)
  
  auto_corr_t_med <- apply(auto_corr_t,1,median)
  auto_corr_t_q025 <- apply(auto_corr_t,1,\(x)quantile(x, .025))
  auto_corr_t_q975 <- apply(auto_corr_t,1,\(x)quantile(x, .975))
  
  auto_corr_t_adapt_med <- apply(auto_corr_t_adapt,1,median)
  auto_corr_t_adapt_q025 <- apply(auto_corr_t_adapt,1,\(x)quantile(x, .025))
  auto_corr_t_adapt_q975 <- apply(auto_corr_t_adapt,1,\(x)quantile(x, .975))
  
  ggplot() +
    geom_line(aes(x = seq(0,400,10), y = auto_corr_t_med), col = "red") +
    geom_line(aes(x = seq(0,400,10), y = auto_corr_t_adapt_med), col = "green")+
    geom_ribbon(aes(x = seq(0,400,10), ymin = auto_corr_t_q025, ymax = auto_corr_t_q975), fill = "red", alpha = 0.2)+
    geom_ribbon(aes(x = seq(0,400,10), ymin = auto_corr_t_adapt_q025, ymax = auto_corr_t_adapt_q975), fill = "green", alpha = 0.2)
  
  auto_corr_frank_med <- apply(auto_corr_frank,1,median)
  auto_corr_frank_q025 <- apply(auto_corr_frank,1,\(x)quantile(x, .025))
  auto_corr_frank_q975 <- apply(auto_corr_frank,1,\(x)quantile(x, .975))
  
  auto_corr_frank_adapt_med <- apply(auto_corr_frank_adapt,1,median)
  auto_corr_frank_adapt_q025 <- apply(auto_corr_frank_adapt,1,\(x)quantile(x, .025))
  auto_corr_frank_adapt_q975 <- apply(auto_corr_frank_adapt,1,\(x)quantile(x, .975))
  
  ggplot() +
    geom_line(aes(x = seq(0,400,10), y = auto_corr_frank_med), col = "red") +
    geom_line(aes(x = seq(0,400,10), y = auto_corr_frank_adapt_med), col = "green")+
    geom_ribbon(aes(x = seq(0,400,10), ymin = auto_corr_frank_q025, ymax = auto_corr_frank_q975), fill = "red", alpha = 0.2)+
    geom_ribbon(aes(x = seq(0,400,10), ymin = auto_corr_frank_adapt_q025, ymax = auto_corr_frank_adapt_q975), fill = "green", alpha = 0.2)
  
  
  
}
