# packages
library(VineCopula)

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
n.born.out.par <- 500
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

################################################################################
# gaussian
################################################################################
if(T){
  n.tree <- 5
  
  for (i in 2) {
    assign(paste0("gauss_mcmc_",i,"_tree_",n.tree,"_adapt"), lapply(1:R, \(x)multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                                                      n.tree = n.tree, n.chain = n.chain_par, n.cores = 1,
                                                                                                      X = X_obs.norm[[x]],
                                                                                                      U1 = get(paste0("copula_uu_gauss_",i))[[x]][,1],
                                                                                                      U2 = get(paste0("copula_uu_gauss_",i))[[x]][,2],
                                                                                                      prior_list = lb.prior.def, 
                                                                                                      moves.prob = moves.prob_par, 
                                                                                                      starting.tree = NULL,
                                                                                                      cont.unif = cont.unif_par,
                                                                                                      include.split = incl.split_par,
                                                                                                      prop_mu = 0, prop_sigma = .2,
                                                                                                      theta_param_1 = 0, theta_param_2 = 1,
                                                                                                      var_param_1 = 1, var_param_2 = 2,
                                                                                                      prior_type = "N",
                                                                                                      cop_type = "gauss",
                                                                                                      adapt = T)))
    
    cat('done case', i, '\n')
    
    save(list = paste0("gauss_mcmc_",i,"_tree_",n.tree,"_adapt"), file = paste0("gauss_mcmc_",i,"_tree_",n.tree,"_adapt", ".Rdata"))
    rm(list = paste0("gauss_mcmc_",i,"_tree_",n.tree,"_adapt"))
    gc()
  }
  gc()
  for (i in 2) {
    assign(paste0("gauss_mcmc_",i,"_tree_",n.tree), lapply(1:R, \(x)multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                                                    n.tree = n.tree, n.chain = n.chain_par, n.cores = 1,
                                                                                                    X = X_obs.norm[[x]],
                                                                                                    U1 = get(paste0("copula_uu_gauss_",i))[[x]][,1],
                                                                                                    U2 = get(paste0("copula_uu_gauss_",i))[[x]][,2],
                                                                                                    prior_list = lb.prior.def, 
                                                                                                    moves.prob = moves.prob_par, 
                                                                                                    starting.tree = NULL,
                                                                                                    cont.unif = cont.unif_par,
                                                                                                    include.split = incl.split_par,
                                                                                                    prop_mu = 0, prop_sigma = .2,
                                                                                                    theta_param_1 = 0, theta_param_2 = 1,
                                                                                                    var_param_1 = 1, var_param_2 = 2,
                                                                                                    prior_type = "N",
                                                                                                    cop_type = "gauss",
                                                                                                    adapt = F)))
    
    cat('done case', i, '\n')
    
    save(list = paste0("gauss_mcmc_",i,"_tree_",n.tree), file = paste0("gauss_mcmc_",i,"_tree_",n.tree, ".Rdata"))
    rm(list = paste0("gauss_mcmc_",i,"_tree_",n.tree))
    gc()
  }
  gc()
}
gc()
################################################################################
# t
################################################################################
if(T){
  n.tree <- 5
  
  for (i in 2) {
    assign(paste0("t_mcmc_",i,"_tree_",n.tree,"_adapt"), lapply(1:R, \(x)multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                                                         n.tree = n.tree, n.chain = n.chain_par, n.cores = 1,
                                                                                                         X = X_obs.norm[[x]],
                                                                                                         U1 = get(paste0("copula_uu_t_",i))[[x]][,1],
                                                                                                         U2 = get(paste0("copula_uu_t_",i))[[x]][,2],
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
                                                                                                         adapt = T)))
    
    cat('done case', i, '\n')
    
    save(list = paste0("t_mcmc_",i,"_tree_",n.tree,"_adapt"), file = paste0("t_mcmc_",i,"_tree_",n.tree,"_adapt", ".Rdata"))
    rm(list = paste0("t_mcmc_",i,"_tree_",n.tree,"_adapt"))
    gc()
  }
  gc()
  for (i in 2) {
    assign(paste0("t_mcmc_",i,"_tree_",n.tree), lapply(1:R, \(x)multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                                                n.tree = n.tree, n.chain = n.chain_par, n.cores = 1,
                                                                                                X = X_obs.norm[[x]],
                                                                                                U1 = get(paste0("copula_uu_t_",i))[[x]][,1],
                                                                                                U2 = get(paste0("copula_uu_t_",i))[[x]][,2],
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
                                                                                                adapt = F)))
    
    cat('done case', i, '\n')
    
    save(list = paste0("t_mcmc_",i,"_tree_",n.tree), file = paste0("t_mcmc_",i,"_tree_",n.tree, ".Rdata"))
    rm(list = paste0("t_mcmc_",i,"_tree_",n.tree))
    gc()
  }
  gc()
}
gc()
################################################################################
# clayton
################################################################################
if(T){
  n.tree <- 5
  
  for (i in 2) {
    assign(paste0("clayton_mcmc_",i,"_tree_",n.tree,"_adapt"), lapply(1:R, \(x)multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                                                               n.tree = n.tree, n.chain = n.chain_par, n.cores = 1,
                                                                                                               X = X_obs.norm[[x]],
                                                                                                               U1 = get(paste0("copula_uu_clayton_",i))[[x]][,1],
                                                                                                               U2 = get(paste0("copula_uu_clayton_",i))[[x]][,2],
                                                                                                               prior_list = lb.prior.def, 
                                                                                                               moves.prob = moves.prob_par, 
                                                                                                               starting.tree = NULL,
                                                                                                               cont.unif = cont.unif_par,
                                                                                                               include.split = incl.split_par,
                                                                                                               prop_mu = 0, prop_sigma = .2,
                                                                                                               theta_param_1 = 0, theta_param_2 = 1,
                                                                                                               var_param_1 = 1, var_param_2 = 2,
                                                                                                               prior_type = "N",
                                                                                                               cop_type = "clayton",
                                                                                                               adapt = T)))
    
    cat('done case', i, '\n')
    
    save(list = paste0("clayton_mcmc_",i,"_tree_",n.tree,"_adapt"), file = paste0("clayton_mcmc_",i,"_tree_",n.tree,"_adapt", ".Rdata"))
    rm(list = paste0("clayton_mcmc_",i,"_tree_",n.tree,"_adapt"))
    gc()
  }
  gc()
  for (i in 2) {
    assign(paste0("clayton_mcmc_",i,"_tree_",n.tree), lapply(1:R, \(x)multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                                                      n.tree = n.tree, n.chain = n.chain_par, n.cores = 1,
                                                                                                      X = X_obs.norm[[x]],
                                                                                                      U1 = get(paste0("copula_uu_clayton_",i))[[x]][,1],
                                                                                                      U2 = get(paste0("copula_uu_clayton_",i))[[x]][,2],
                                                                                                      prior_list = lb.prior.def, 
                                                                                                      moves.prob = moves.prob_par, 
                                                                                                      starting.tree = NULL,
                                                                                                      cont.unif = cont.unif_par,
                                                                                                      include.split = incl.split_par,
                                                                                                      prop_mu = 0, prop_sigma = .2,
                                                                                                      theta_param_1 = 0, theta_param_2 = 1,
                                                                                                      var_param_1 = 1, var_param_2 = 2,
                                                                                                      prior_type = "N",
                                                                                                      cop_type = "clayton",
                                                                                                      adapt = F)))
    
    cat('done case', i, '\n')
    
    save(list = paste0("clayton_mcmc_",i,"_tree_",n.tree), file = paste0("clayton_mcmc_",i,"_tree_",n.tree, ".Rdata"))
    rm(list = paste0("clayton_mcmc_",i,"_tree_",n.tree))
    gc()
  }
  gc()
}
gc()
################################################################################
# Gumbel
################################################################################
if(T){
  n.tree <- 5
  
  for (i in 2) {
    assign(paste0("gumbel_mcmc_",i,"_tree_",n.tree,"_adapt"), lapply(1:R, \(x)multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                                                              n.tree = n.tree, n.chain = n.chain_par, n.cores = 1,
                                                                                                              X = X_obs.norm[[x]],
                                                                                                              U1 = get(paste0("copula_uu_gumbel_",i))[[x]][,1],
                                                                                                              U2 = get(paste0("copula_uu_gumbel_",i))[[x]][,2],
                                                                                                              prior_list = lb.prior.def, 
                                                                                                              moves.prob = moves.prob_par, 
                                                                                                              starting.tree = NULL,
                                                                                                              cont.unif = cont.unif_par,
                                                                                                              include.split = incl.split_par,
                                                                                                              prop_mu = 0, prop_sigma = .2,
                                                                                                              theta_param_1 = 0, theta_param_2 = 1,
                                                                                                              var_param_1 = 1, var_param_2 = 2,
                                                                                                              prior_type = "N",
                                                                                                              cop_type = "gumbel",
                                                                                                              adapt = T)))
    
    cat('done case', i, '\n')
    
    save(list = paste0("gumbel_mcmc_",i,"_tree_",n.tree,"_adapt"), file = paste0("gumbel_mcmc_",i,"_tree_",n.tree,"_adapt", ".Rdata"))
    rm(list = paste0("gumbel_mcmc_",i,"_tree_",n.tree,"_adapt"))
    gc()
  }
  gc()
  for (i in 2) {
    assign(paste0("gumbel_mcmc_",i,"_tree_",n.tree), lapply(1:R, \(x)multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                                                     n.tree = n.tree, n.chain = n.chain_par, n.cores = 1,
                                                                                                     X = X_obs.norm[[x]],
                                                                                                     U1 = get(paste0("copula_uu_gumbel_",i))[[x]][,1],
                                                                                                     U2 = get(paste0("copula_uu_gumbel_",i))[[x]][,2],
                                                                                                     prior_list = lb.prior.def, 
                                                                                                     moves.prob = moves.prob_par, 
                                                                                                     starting.tree = NULL,
                                                                                                     cont.unif = cont.unif_par,
                                                                                                     include.split = incl.split_par,
                                                                                                     prop_mu = 0, prop_sigma = .2,
                                                                                                     theta_param_1 = 0, theta_param_2 = 1,
                                                                                                     var_param_1 = 1, var_param_2 = 2,
                                                                                                     prior_type = "N",
                                                                                                     cop_type = "gumbel",
                                                                                                     adapt = F)))
    
    cat('done case', i, '\n')
    
    save(list = paste0("gumbel_mcmc_",i,"_tree_",n.tree), file = paste0("gumbel_mcmc_",i,"_tree_",n.tree, ".Rdata"))
    rm(list = paste0("gumbel_mcmc_",i,"_tree_",n.tree))
    gc()
  }
  gc()
}
gc()
################################################################################
# frank
################################################################################
if(T){
  n.tree <- 5
  
  for (i in 2) {
    assign(paste0("frank_mcmc_",i,"_tree_",n.tree,"_adapt"), lapply(1:R, \(x)multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                                                             n.tree = n.tree, n.chain = n.chain_par, n.cores = 1,
                                                                                                             X = X_obs.norm[[x]],
                                                                                                             U1 = get(paste0("copula_uu_frank_",i))[[x]][,1],
                                                                                                             U2 = get(paste0("copula_uu_frank_",i))[[x]][,2],
                                                                                                             prior_list = lb.prior.def, 
                                                                                                             moves.prob = moves.prob_par, 
                                                                                                             starting.tree = NULL,
                                                                                                             cont.unif = cont.unif_par,
                                                                                                             include.split = incl.split_par,
                                                                                                             prop_mu = 0, prop_sigma = .2,
                                                                                                             theta_param_1 = 0, theta_param_2 = 1,
                                                                                                             var_param_1 = 1, var_param_2 = 2,
                                                                                                             prior_type = "N",
                                                                                                             cop_type = "frank",
                                                                                                             adapt = T)))
    
    cat('done case', i, '\n')
    
    save(list = paste0("frank_mcmc_",i,"_tree_",n.tree,"_adapt"), file = paste0("frank_mcmc_",i,"_tree_",n.tree,"_adapt", ".Rdata"))
    rm(list = paste0("frank_mcmc_",i,"_tree_",n.tree,"_adapt"))
    gc()
  }
  gc()
  for (i in 2) {
    assign(paste0("frank_mcmc_",i,"_tree_",n.tree), lapply(1:R, \(x)multichain_MCMC_copula(n.iter = n.iter_par, n.burn = n.born.out.par,
                                                                                                    n.tree = n.tree, n.chain = n.chain_par, n.cores = 1,
                                                                                                    X = X_obs.norm[[x]],
                                                                                                    U1 = get(paste0("copula_uu_frank_",i))[[x]][,1],
                                                                                                    U2 = get(paste0("copula_uu_frank_",i))[[x]][,2],
                                                                                                    prior_list = lb.prior.def, 
                                                                                                    moves.prob = moves.prob_par, 
                                                                                                    starting.tree = NULL,
                                                                                                    cont.unif = cont.unif_par,
                                                                                                    include.split = incl.split_par,
                                                                                                    prop_mu = 0, prop_sigma = .2,
                                                                                                    theta_param_1 = 0, theta_param_2 = 1,
                                                                                                    var_param_1 = 1, var_param_2 = 2,
                                                                                                    prior_type = "N",
                                                                                                    cop_type = "frank",
                                                                                                    adapt = F)))
    
    cat('done case', i, '\n')
    
    save(list = paste0("frank_mcmc_",i,"_tree_",n.tree), file = paste0("frank_mcmc_",i,"_tree_",n.tree, ".Rdata"))
    rm(list = paste0("frank_mcmc_",i,"_tree_",n.tree))
    gc()
  }
  gc()
}
gc()