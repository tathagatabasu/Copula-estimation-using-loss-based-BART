# packages
library(VineCopula)
library(dplyr)

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

# plot(X_obs[[1]], tau_true_1[[1]], xlab = "Observations", ylab = "tau")
# 
# plot(X_obs[[1]], tau_true_2[[1]], xlab = "Observations", ylab = "tau")

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

n.tree <- 5
test_case <- 2

################################################################################
# source files
################################################################################

source('MCMC_BART_copula.R')
source('import_functions.R')

moves.prob_par <- c(0.1, 0.4, 0.25, 0.25)
lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) 

################################################################################
load(paste0("mcmc_2_tree_5_dat_adapt/gauss_mcmc_",test_case,"_tree_",n.tree,"_dat", ".Rdata"))

load(paste0("mcmc_2_tree_5_dat_adapt/gauss_mcmc_",test_case,"_tree_",n.tree,"_dat_adapt", ".Rdata"))

gauss_theta <- list()

for (i in 1:R) {
  gauss_theta[[i]] <- gauss_dat_list$pred[[i]] %>%
    group_by(obs, theta_true) %>%
    summarise(mean = mean(theta_mean)) 
}

gauss_theta_adapt <- list()

for (i in 1:R) {
  gauss_theta_adapt[[i]] <- gauss_dat_list_adapt$pred[[i]] %>%
    group_by(obs, theta_true) %>%
    summarise(mean = mean(theta_mean)) 
}

rm(gauss_dat_list,gauss_dat_list_adapt)

gauss_pred_cop <- list()

for (i in 1:R) {
  gauss_pred_cop[[i]] <- BiCopSim(n, par = BiCopTau2Par(1, gauss_theta[[i]]$mean), family = 1)
}

gauss_pred_cop_adapt <- list()

for (i in 1:R) {
  gauss_pred_cop_adapt[[i]] <- BiCopSim(n, par = BiCopTau2Par(1, gauss_theta_adapt[[i]]$mean), family = 1)
}

library(cramer)

cram_gauss <- list()

for (j in 1:R) {
  cram_gauss[[j]] <- sapply(1:100, function(i){set.seed(i); cramer.test(get(paste0("copula_uu_gauss_",test_case))[[j]], gauss_pred_cop[[j]], replicates = 100, sim = "permutation")$p.value})
}

mean(unlist(lapply(cram_gauss, mean)))

cram_gauss_adapt <- list()

for (j in 1:R) {
  cram_gauss_adapt[[j]] <- sapply(1:100, function(i){set.seed(i); cramer.test(get(paste0("copula_uu_gauss_",test_case))[[j]], gauss_pred_cop_adapt[[j]], replicates = 100, sim = "permutation")$p.value})
}

mean(unlist(lapply(cram_gauss_adapt, mean)))

library(fasano.franceschini.test)

ff_gauss <- list()

for (j in 1:R) {
  ff_gauss[[j]] <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(get(paste0("copula_uu_gauss_",test_case))[[j]], gauss_pred_cop[[j]], nPermute = 100, verbose = F)$p.value})
}

mean(unlist(lapply(ff_gauss, mean)))

ff_gauss_adapt <- list()

for (j in 1:R) {
  ff_gauss_adapt[[j]] <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(get(paste0("copula_uu_gauss_",test_case))[[j]], gauss_pred_cop_adapt[[j]], nPermute = 100, verbose = F)$p.value})
}

mean(unlist(lapply(ff_gauss_adapt, mean)))

save(cram_gauss, cram_gauss_adapt, ff_gauss, ff_gauss_adapt, file = paste("gauss_gof",test_case,".Rdata"))

################################################################################
load(paste0("mcmc_2_tree_5_dat_adapt/t_mcmc_",test_case,"_tree_",n.tree,"_dat", ".Rdata"))

load(paste0("mcmc_2_tree_5_dat_adapt/t_mcmc_",test_case,"_tree_",n.tree,"_dat_adapt", ".Rdata"))

t_theta <- list()

for (i in 1:R) {
  t_theta[[i]] <- t_dat_list$pred[[i]] %>%
    group_by(obs, theta_true) %>%
    summarise(mean = mean(theta_mean)) 
}

t_theta_adapt <- list()

for (i in 1:R) {
  t_theta_adapt[[i]] <- t_dat_list_adapt$pred[[i]] %>%
    group_by(obs, theta_true) %>%
    summarise(mean = mean(theta_mean)) 
}

rm(t_dat_list,t_dat_list_adapt)

t_pred_cop <- list()

for (i in 1:R) {
  t_pred_cop[[i]] <- BiCopSim(n, par = BiCopTau2Par(2, t_theta[[i]]$mean), family = 2, par2 = 3)
}

t_pred_cop_adapt <- list()

for (i in 1:R) {
  t_pred_cop_adapt[[i]] <- BiCopSim(n, par = BiCopTau2Par(2, t_theta_adapt[[i]]$mean), family = 2, par2 = 3)
}

library(cramer)

cram_t <- list()

for (j in 1:R) {
  cram_t[[j]] <- sapply(1:100, function(i){set.seed(i); cramer.test(get(paste0("copula_uu_t_",test_case))[[j]], t_pred_cop[[j]], replicates = 100, sim = "permutation")$p.value})
}

median(unlist(lapply(cram_t, mean)))

cram_t_adapt <- list()

for (j in 1:R) {
  cram_t_adapt[[j]] <- sapply(1:100, function(i){set.seed(i); cramer.test(get(paste0("copula_uu_t_",test_case))[[j]], t_pred_cop_adapt[[j]], replicates = 100, sim = "permutation")$p.value})
}

median(unlist(lapply(cram_t_adapt, mean)))

library(fasano.franceschini.test)

ff_t <- list()

for (j in 1:R) {
  ff_t[[j]] <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(get(paste0("copula_uu_t_",test_case))[[j]], t_pred_cop[[j]], nPermute = 100, verbose = F)$p.value})
}

mean(unlist(lapply(ff_t, mean)))

ff_t_adapt <- list()

for (j in 1:R) {
  ff_t_adapt[[j]] <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(get(paste0("copula_uu_t_",test_case))[[j]], t_pred_cop_adapt[[j]], nPermute = 100, verbose = F)$p.value})
}

mean(unlist(lapply(ff_t_adapt, mean)))

save(cram_t, cram_t_adapt, ff_t, ff_t_adapt, file = paste("t_gof",test_case,".Rdata"))

################################################################################
load(paste0("mcmc_2_tree_5_dat_adapt/clayton_mcmc_",test_case,"_tree_",n.tree,"_dat", ".Rdata"))

load(paste0("mcmc_2_tree_5_dat_adapt/clayton_mcmc_",test_case,"_tree_",n.tree,"_dat_adapt", ".Rdata"))

clayton_theta <- list()

for (i in 1:R) {
  clayton_theta[[i]] <- clayton_dat_list$pred[[i]] %>%
    group_by(obs, theta_true) %>%
    summarise(mean = mean(theta_mean)) 
}

clayton_theta_adapt <- list()

for (i in 1:R) {
  clayton_theta_adapt[[i]] <- clayton_dat_list_adapt$pred[[i]] %>%
    group_by(obs, theta_true) %>%
    summarise(mean = mean(theta_mean)) 
}

rm(clayton_dat_list,clayton_dat_list_adapt)

clayton_pred_cop <- list()

for (i in 1:R) {
  clayton_pred_cop[[i]] <- BiCopSim(n, par = BiCopTau2Par(3, clayton_theta[[i]]$mean), family = 3)
}

clayton_pred_cop_adapt <- list()

for (i in 1:R) {
  clayton_pred_cop_adapt[[i]] <- BiCopSim(n, par = BiCopTau2Par(3, clayton_theta_adapt[[i]]$mean), family = 3)
}

library(cramer)

cram_clayton <- list()

for (j in 1:R) {
  cram_clayton[[j]] <- sapply(1:100, function(i){set.seed(i); cramer.test(get(paste0("copula_uu_clayton_",test_case))[[j]], clayton_pred_cop[[j]], replicates = 100, sim = "permutation")$p.value})
}

median(unlist(lapply(cram_clayton, mean)))

cram_clayton_adapt <- list()

for (j in 1:R) {
  cram_clayton_adapt[[j]] <- sapply(1:100, function(i){set.seed(i); cramer.test(get(paste0("copula_uu_clayton_",test_case))[[j]], clayton_pred_cop_adapt[[j]], replicates = 100, sim = "permutation")$p.value})
}

median(unlist(lapply(cram_clayton_adapt, mean)))

library(fasano.franceschini.test)

ff_clayton <- list()

for (j in 1:R) {
  ff_clayton[[j]] <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(get(paste0("copula_uu_clayton_",test_case))[[j]], clayton_pred_cop[[j]], nPermute = 100, verbose = F)$p.value})
}

mean(unlist(lapply(ff_clayton, mean)))

ff_clayton_adapt <- list()

for (j in 1:R) {
  ff_clayton_adapt[[j]] <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(get(paste0("copula_uu_clayton_",test_case))[[j]], clayton_pred_cop_adapt[[j]], nPermute = 100, verbose = F)$p.value})
}

mean(unlist(lapply(ff_clayton_adapt, mean)))

save(cram_clayton, cram_clayton_adapt, ff_clayton, ff_clayton_adapt, file = paste("clayton_gof",test_case,".Rdata"))

################################################################################
load(paste0("mcmc_2_tree_5_dat_adapt/gumbel_mcmc_",test_case,"_tree_",n.tree,"_dat", ".Rdata"))

load(paste0("mcmc_2_tree_5_dat_adapt/gumbel_mcmc_",test_case,"_tree_",n.tree,"_dat_adapt", ".Rdata"))

gumbel_theta <- list()

for (i in 1:R) {
  gumbel_theta[[i]] <- gumbel_dat_list$pred[[i]] %>%
    group_by(obs, theta_true) %>%
    summarise(mean = mean(theta_mean)) 
}

gumbel_theta_adapt <- list()

for (i in 1:R) {
  gumbel_theta_adapt[[i]] <- gumbel_dat_list_adapt$pred[[i]] %>%
    group_by(obs, theta_true) %>%
    summarise(mean = mean(theta_mean)) 
}

rm(gumbel_dat_list,gumbel_dat_list_adapt)

gumbel_pred_cop <- list()

for (i in 1:R) {
  gumbel_pred_cop[[i]] <- BiCopSim(n, par = BiCopTau2Par(4, gumbel_theta[[i]]$mean), family = 4)
}

gumbel_pred_cop_adapt <- list()

for (i in 1:R) {
  gumbel_pred_cop_adapt[[i]] <- BiCopSim(n, par = BiCopTau2Par(4, gumbel_theta_adapt[[i]]$mean), family = 4)
}

library(cramer)

cram_gumbel <- list()

for (j in 1:R) {
  cram_gumbel[[j]] <- sapply(1:100, function(i){set.seed(i); cramer.test(get(paste0("copula_uu_gumbel_",test_case))[[j]], gumbel_pred_cop[[j]], replicates = 100, sim = "permutation")$p.value})
}

median(unlist(lapply(cram_gumbel, mean)))

cram_gumbel_adapt <- list()

for (j in 1:R) {
  cram_gumbel_adapt[[j]] <- sapply(1:100, function(i){set.seed(i); cramer.test(get(paste0("copula_uu_gumbel_",test_case))[[j]], gumbel_pred_cop_adapt[[j]], replicates = 100, sim = "permutation")$p.value})
}

median(unlist(lapply(cram_gumbel_adapt, mean)))

library(fasano.franceschini.test)

ff_gumbel <- list()

for (j in 1:R) {
  ff_gumbel[[j]] <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(get(paste0("copula_uu_gumbel_",test_case))[[j]], gumbel_pred_cop[[j]], nPermute = 100, verbose = F)$p.value})
}

mean(unlist(lapply(ff_gumbel, mean)))

ff_gumbel_adapt <- list()

for (j in 1:R) {
  ff_gumbel_adapt[[j]] <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(get(paste0("copula_uu_gumbel_",test_case))[[j]], gumbel_pred_cop_adapt[[j]], nPermute = 100, verbose = F)$p.value})
}

mean(unlist(lapply(ff_gumbel_adapt, mean)))

save(cram_gumbel, cram_gumbel_adapt, ff_gumbel, ff_gumbel_adapt, file = paste("gumbel_gof",test_case,".Rdata"))

################################################################################
# load(paste0("mcmc_2_tree_5_dat_adapt/frank_mcmc_",test_case,"_tree_",n.tree,"_dat", "_new_prop.Rdata"))
# 
# load(paste0("mcmc_2_tree_5_dat_adapt/frank_mcmc_",test_case,"_tree_",n.tree,"_dat_adapt", "_new_prop.Rdata"))

load(paste0("mcmc_2_tree_5_dat_adapt/frank_mcmc_",test_case,"_tree_",n.tree,"_dat", ".Rdata"))

load(paste0("mcmc_2_tree_5_dat_adapt/frank_mcmc_",test_case,"_tree_",n.tree,"_dat_adapt", ".Rdata"))

frank_theta <- list()

for (i in 1:R) {
  frank_theta[[i]] <- frank_dat_list$pred[[i]] %>%
    group_by(obs, theta_true) %>%
    summarise(mean = mean(theta_mean)) 
}

frank_theta_adapt <- list()

for (i in 1:R) {
  frank_theta_adapt[[i]] <- frank_dat_list_adapt$pred[[i]] %>%
    group_by(obs, theta_true) %>%
    summarise(mean = mean(theta_mean)) 
}

rm(frank_dat_list,frank_dat_list_adapt)

frank_pred_cop <- list()

for (i in 1:R) {
  frank_pred_cop[[i]] <- BiCopSim(n, par = BiCopTau2Par(5, frank_theta[[i]]$mean), family = 5)
}

frank_pred_cop_adapt <- list()

for (i in 1:R) {
  frank_pred_cop_adapt[[i]] <- BiCopSim(n, par = BiCopTau2Par(5, frank_theta_adapt[[i]]$mean), family = 5)
}

library(cramer)

cram_frank <- list()

for (j in 1:R) {
  cram_frank[[j]] <- sapply(1:100, function(i){set.seed(i); cramer.test(get(paste0("copula_uu_frank_",test_case))[[j]], frank_pred_cop[[j]], replicates = 100, sim = "permutation")$p.value})
}

median(unlist(lapply(cram_frank, mean)))

cram_frank_adapt <- list()

for (j in 1:R) {
  cram_frank_adapt[[j]] <- sapply(1:100, function(i){set.seed(i); cramer.test(get(paste0("copula_uu_frank_",test_case))[[j]], frank_pred_cop_adapt[[j]], replicates = 100, sim = "permutation")$p.value})
}

median(unlist(lapply(cram_frank_adapt, mean)))

library(fasano.franceschini.test)

ff_frank <- list()

for (j in 1:R) {
  ff_frank[[j]] <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(get(paste0("copula_uu_frank_",test_case))[[j]], frank_pred_cop[[j]], nPermute = 100, verbose = F)$p.value})
}

mean(unlist(lapply(ff_frank, mean)))

ff_frank_adapt <- list()

for (j in 1:R) {
  ff_frank_adapt[[j]] <- sapply(1:100, function(i){set.seed(i); fasano.franceschini.test(get(paste0("copula_uu_frank_",test_case))[[j]], frank_pred_cop_adapt[[j]], nPermute = 100, verbose = F)$p.value})
}

mean(unlist(lapply(ff_frank_adapt, mean)))

save(cram_frank, cram_frank_adapt, ff_frank, ff_frank_adapt, file = paste("frank_gof",test_case,".Rdata"))

