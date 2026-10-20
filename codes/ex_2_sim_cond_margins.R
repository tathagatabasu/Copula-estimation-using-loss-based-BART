# packages
library(VineCopula)

################################################################################
# Simulation with covariate-dependent margins (Reviewer 2, Major comment 3)
#
#   x ~ U(-1, 1),  (U1, U2) | x ~ C(. | tau(x)),  tau(x) as in Section 5
#   (x rescaled to [0,1]),  Y_j = F_{j|x}^{-1}(U_j)  with three marginal designs:
#
#   marg 1 (location, Reviewer 2):  Y1 | x ~ N(x, 1),  Y2 | x ~ N(|x|, 1)
#   marg 2 (nonlinear location-scale):
#            Y1 = 2 sin(pi x) + (0.5 + |x|) Z1,   Y2 = 3 x^2 + exp(x/2) Z2,  Z_j = Phi^{-1}(U_j)
#   marg 3 (non-Gaussian):
#            Y1 = exp(x) * G,  G = Gamma(2,1) quantile of U1   (skewed, scale depends on x)
#            Y2 = 2 x + T,     T = t_3 quantile of U2           (heavy tails)
#
# In every design the global ranks pobs(Y_j) are functions of x, so the
# naive approach mixes the marginal covariate effect into the copula. The true
# conditional pseudo-observations U_j = F_{j|x}(Y_j) are uniform and free of x.
# Design 3 is misspecified for the joint model (Gaussian margins), but not for the
# rank-based two-step estimators.
#
# Methods (each run in parallel over the R replicates):
#   naive  : global ranks pobs(Y)              -> multichain_MCMC_copula
#   bart2s : two-step, BART conditional margins -> multichain_MCMC_copula
#   gam2s  : two-step, GAM  conditional margins -> multichain_MCMC_copula
#   joint  : full model (margins + copula in one MCMC)
################################################################################

################################################################################
# data generation
################################################################################
set.seed(1e3)

n <- 200
R <- 10
X_obs <- lapply(1:R, function(x)matrix(runif(n, -1, 1), ncol = 1))

# normalise predictors
data_norm <- function(X_obs){
  X_obs.norm <- as.data.frame(apply(X_obs, 2, \(x) (x - min(x))/(max(x) - min(x))))
  X_obs.norm <- as.matrix(X_obs.norm)
  rownames(X_obs.norm) <- 1:nrow(X_obs)
  return(X_obs.norm)
}

X_obs.norm <- lapply(1:R, function(i)data_norm(X_obs[[i]]))

# tau with tree structure (on the normalised predictor)
tree_tau_func <- function(X_obs){
  tau_true_1 <- rep(0,nrow(X_obs))
  tau_true_1[X_obs<0.33] <- 0.3
  tau_true_1[(X_obs>=0.33)&(X_obs<0.66)] <- 0.8
  tau_true_1[(X_obs>=0.66)] <- 0.3
  
  tau_true_1 <- matrix(tau_true_1, ncol = 1)
}

tau_true_1 <- lapply(1:R, function(i){tree_tau_func(X_obs.norm[[i]])})

# periodic
tau_true_2 <- lapply(1:R, function(i){0.2*sin(2*pi*X_obs.norm[[i]]) + 0.5})

plot(X_obs.norm[[1]], tau_true_1[[1]], xlab = "Observations", ylab = "tau")

plot(X_obs.norm[[1]], tau_true_2[[1]], xlab = "Observations", ylab = "tau")

# copula data (identical to ex_0_sim_cond_copula.R)

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

# conditional margins: Y_j = F_{j|x}^{-1}(U_j)

marg_transform <- function(uu, x, m){
  x <- x[,1]
  if(m == 1){
    cbind(Y1 = x + qnorm(uu[,1]),
          Y2 = abs(x) + qnorm(uu[,2]))
  } else if(m == 2){
    cbind(Y1 = 2*sin(pi*x) + (0.5 + abs(x))*qnorm(uu[,1]),
          Y2 = 3*x^2 + exp(x/2)*qnorm(uu[,2]))
  } else if(m == 3){
    cbind(Y1 = exp(x)*qgamma(uu[,1], shape = 2, rate = 1),
          Y2 = 2*x + qt(uu[,2], df = 3))
  } else stop("Unknown marginal design")
}

for (cop in c("gauss", "t", "gumbel", "clayton", "frank")) {
  for (i in 1:2) {
    for (m in 1:3) {
      assign(paste0("y_", cop, "_", i, "_m", m), lapply(1:R, function(k){marg_transform(get(paste0("copula_uu_", cop, "_", i))[[k]], X_obs[[k]], m)}))
    }
  }
}

# Y against x, and the naive global ranks against x (they depend on x), for each design
par(mfrow = c(3, 3))
for (m in 1:3) {
  yy <- get(paste0("y_gauss_1_m", m))[[1]]
  plot(X_obs[[1]], yy[,1], xlab = "x", ylab = "Y1", main = paste("design", m))
  plot(X_obs[[1]], yy[,2], xlab = "x", ylab = "Y2", main = paste("design", m))
  plot(X_obs[[1]], pobs(yy[,1]), xlab = "x", ylab = "pobs(Y1)", main = "global ranks vs x")
}
par(mfrow = c(1, 1))

# correlation with x: true conditional U1 (~0) vs global ranks pobs(Y1) (clearly non-zero)
print(sapply(1:3, function(m) c(cor_U1_x = cor(copula_uu_gauss_1[[1]][,1], X_obs[[1]][,1]),
                                cor_pobsY1_x = cor(pobs(get(paste0("y_gauss_1_m", m))[[1]][,1]), X_obs[[1]][,1]))))

# mcmc params

n.chain_par <- 4
n.iter_par <- 3000
n.born.out.par <- 500
n.thin <- 1

incl.split_par <- TRUE
cont.unif_par <- TRUE

n.tree <- 1          # copula trees
n.tree.marg <- 5     # trees per margin (joint model only)
n.cores_par <- 10    # replicates run in parallel

case_list <- 1                                  # 1 = tree tau, 2 = periodic tau
marg_list <- 1:3                                # marginal designs (see top)
adapt_list <- c(TRUE, FALSE)                    # A-C-BART and C-BART
method_list <- c("naive", "bart2s", "gam2s", "joint")

################################################################################
# source files
################################################################################

source('MCMC_BART_copula.R')
source('import_functions.R')
source('MCMC_BART_copula_margins.R')

moves.prob_par <- c(0.1, 0.4, 0.25, 0.25)
lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944))

################################################################################
# one call = one method for one copula / case / adapt, parallel over replicates
################################################################################

run_sim <- function(cop, i, m, method, adapt){
  
  Y_list <- get(paste0("y_", cop, "_", i, "_m", m))
  
  # arguments shared by every method
  common <- list(n.iter = n.iter_par, n.burn = n.born.out.par,
                 n.tree = n.tree, n.chain = n.chain_par,
                 prior_list = lb.prior.def,
                 moves.prob = moves.prob_par,
                 cont.unif = cont.unif_par,
                 include.split = incl.split_par,
                 prop_mu = 0, prop_sigma = .2,
                 theta_param_1 = 0, theta_param_2 = 1,
                 var_param_1 = 1, var_param_2 = 2,
                 prior_type = "N",
                 cop_type = cop,
                 adapt = adapt)
  
  fit_one <- switch(method,
                    
                    naive = \(x){
                      tt <- system.time(
                        fit <- do.call(multichain_MCMC_copula,
                                       c(list(X = X_obs.norm[[x]],
                                              U1 = pobs(Y_list[[x]][,1]),
                                              U2 = pobs(Y_list[[x]][,2]),
                                              starting.tree = NULL, n.cores = 1), common)))
                      fit$time <- tt[["elapsed"]]
                      fit
                    },
                    
                    bart2s = \(x) do.call(two_step_MCMC_copula,
                                          c(list(Y1 = Y_list[[x]][,1], Y2 = Y_list[[x]][,2],
                                                 X = X_obs.norm[[x]], X_marg = X_obs.norm[[x]],
                                                 margin_method = "bart", pobs_type = "rank",
                                                 starting.tree = NULL, n.cores = 1), common)),
                    
                    gam2s = \(x) do.call(two_step_MCMC_copula,
                                         c(list(Y1 = Y_list[[x]][,1], Y2 = Y_list[[x]][,2],
                                                X = X_obs.norm[[x]], X_marg = X_obs.norm[[x]],
                                                margin_method = "gam", pobs_type = "rank",
                                                starting.tree = NULL, n.cores = 1), common)),
                    
                    joint = \(x) do.call(multichain_MCMC_copula_joint,
                                         c(list(Y1 = Y_list[[x]][,1], Y2 = Y_list[[x]][,2],
                                                X = X_obs.norm[[x]], X_marg = X_obs.norm[[x]],
                                                n.tree.marg = n.tree.marg,
                                                prop_sigma_marg = 0.1,
                                                n.cores = 1), common)),
                    
                    stop("Unknown method"))
  
  obj_name <- paste0(cop, "_", method, "_mcmc_", i, "_marg_", m, "_tree_", n.tree, if(adapt) "_adapt" else "")
  
  assign(obj_name, mclapply_manual(1:R, fit_one, mc.cores = n.cores_par, seed.list = rep(1,R)))
  
  cat('done', cop, method, 'case', i, 'marg', m, if(adapt) 'adapt' else '', '\n')
  
  save(list = obj_name, file = paste0(obj_name, ".Rdata"))
  rm(list = obj_name)
  gc()
  
  invisible(obj_name)
}

################################################################################
# gaussian
################################################################################
if(F){
  for (i in case_list) {
    for (m in marg_list) {
      for (adapt in adapt_list) {
        for (method in method_list) {
          run_sim("gauss", i, m, method, adapt)
        }
      }
    }
  }
}

################################################################################
# t
################################################################################
if(F){
  for (i in case_list) {
    for (m in marg_list) {
      for (adapt in adapt_list) {
        for (method in method_list) {
          run_sim("t", i, m, method, adapt)
        }
      }
    }
  }
}

################################################################################
# clayton
################################################################################
if(F){
  for (i in case_list) {
    for (m in marg_list) {
      for (adapt in adapt_list) {
        for (method in method_list) {
          run_sim("clayton", i, m, method, adapt)
        }
      }
    }
  }
}

################################################################################
# Gumbel
################################################################################
if(F){
  for (i in case_list) {
    for (m in marg_list) {
      for (adapt in adapt_list) {
        for (method in method_list) {
          run_sim("gumbel", i, m, method, adapt)
        }
      }
    }
  }
}

################################################################################
# frank
################################################################################
if(F){
  for (i in case_list) {
    for (m in marg_list) {
      for (adapt in adapt_list) {
        for (method in method_list) {
          run_sim("frank", i, m, method, adapt)
        }
      }
    }
  }
}

################################################################################
# quick summary (RMSE, CI length, CI coverage, time) from the saved files
################################################################################
if(T){
  n.discard <- 1500
  
  sim_metrics <- function(cop, i, m, method, adapt){
    obj_name <- paste0(cop, "_", method, "_mcmc_", i, "_marg_", m, "_tree_", n.tree, if(adapt) "_adapt" else "")
    load(paste0(obj_name, ".Rdata"))
    fits <- get(obj_name)
    tau_true <- get(paste0("tau_true_", i))
    
    res <- t(sapply(1:R, function(k){
      ts <- tau_summary(fits[[k]], X_obs.norm[[k]], cop, n.iter_par, n.chain_par, n.discard)
      tt <- tau_true[[k]][,1]
      tm <- if(!is.null(fits[[k]]$time)) unname(fits[[k]]$time[length(fits[[k]]$time)]) else NA
      c(RMSE = sqrt(mean((ts$tau_mean - tt)^2)),
        CI_length = mean(ts$tau_q975 - ts$tau_q025),
        CI_cov = mean(tt >= ts$tau_q025 & tt <= ts$tau_q975),
        time_sec = tm)
    }))
    data.frame(copula = cop, case = i, marg = m, method = method, adapt = adapt, t(colMeans(res, na.rm = TRUE)))
  }
  
  summary_tab <- do.call(rbind, lapply(c("gauss"), function(cop)
    do.call(rbind, lapply(case_list, function(i)
      do.call(rbind, lapply(marg_list, function(m)
        do.call(rbind, lapply(adapt_list, function(adapt)
          do.call(rbind, lapply(method_list, function(method) sim_metrics(cop, i, m, method, adapt)))))))))))
  
  print(summary_tab)
  save(summary_tab, file = "cond_margins_summary.Rdata")
}