################################################################################
# Simulation with covariate-dependent margins (Reviewer 2, Major comment 3)
#
#   x ~ U(-1, 1),   Y1 | x ~ N(x, 1),   Y2 | x ~ N(|x|, 1)
#   (U1, U2) | x ~ C(. | tau(x)),  tau(x) as in Section 5 (on x rescaled to [0,1])
#
# Compares four estimators of tau(x):
#   naive     : global ranks pobs(Y)       (approach of the submitted version)
#   2step-BART: BART conditional margins -> copula BART
#   2step-GAM : GAM  conditional margins -> copula BART
#   joint     : full model, margins and copula estimated in one MCMC
################################################################################

library(VineCopula)

source('MCMC_BART_copula.R')
source('import_functions.R')
source('MCMC_BART_copula_margins.R')

################################################################################
# settings
################################################################################
cop_type   <- "gauss"      # gauss, t, clayton, gumbel, frank
test_case  <- 1            # 1 = tree-structured tau, 2 = periodic tau
n          <- 200
R          <- 4          # replicates (use e.g. R <- 2 for a quick check)

n.chain_par <- 4
n.iter_par  <- 3000
n.burn_par  <- 500         # adaptation starts after n.burn_par
n.discard   <- 1500        # discarded for posterior summaries
n.tree      <- 1           # copula trees
n.tree.marg <- 5           # trees per margin (joint model)
n.cores     <- 1

moves.prob_par <- c(0.1, 0.4, 0.25, 0.25)
lb.prior.def   <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944))

cop_args <- list(n.iter = n.iter_par, n.burn = n.burn_par, n.chain = n.chain_par,
                 n.tree = n.tree, n.cores = n.cores,
                 prior_list = lb.prior.def, moves.prob = moves.prob_par,
                 cont.unif = TRUE, include.split = TRUE,
                 prop_mu = 0, prop_sigma = .5,
                 theta_param_1 = 0, theta_param_2 = 1,
                 var_param_1 = 1, var_param_2 = 2,
                 prior_type = "N", cop_type = cop_type, adapt = TRUE)

################################################################################
# data generation
################################################################################
fam_id <- c(gauss = 1, t = 2, clayton = 3, gumbel = 4, frank = 5)[[cop_type]]
tau2par <- function(tau) if(cop_type == "t") BiCopTau2Par(2, tau) else BiCopTau2Par(fam_id, tau)

tau_fun <- function(x01){
  if(test_case == 1) ifelse(x01 <= 0.33, 0.3, ifelse(x01 <= 0.66, 0.8, 0.3))
  else 0.2 * sin(2 * pi * x01) + 0.5
}

sim_data <- function(seed){
  set.seed(seed)
  x   <- runif(n, -1, 1)
  x01 <- (x - min(x)) / (max(x) - min(x))
  tau <- tau_fun(x01)
  uu  <- t(vapply(seq_len(n), function(i)
    BiCopSim(1, family = fam_id, par = tau2par(tau[i]), par2 = if(cop_type == "t") 3 else 0),
    numeric(2)))
  list(x = x, X = prep_X(matrix(x01, ncol = 1)), tau = tau,
       Y1 = x + qnorm(uu[, 1]),
       Y2 = abs(x) + qnorm(uu[, 2]))
}

metrics <- function(ts, tau_true){
  c(RMSE = sqrt(mean((ts$tau_mean - tau_true)^2)),
    CI_length = mean(ts$tau_q975 - ts$tau_q025),
    CI_cov = mean(tau_true >= ts$tau_q025 & tau_true <= ts$tau_q975))
}

################################################################################
# run
################################################################################
results <- list()

for(r in seq_len(R)){
  d <- sim_data(1000 + r)

  # naive: global ranks (submitted version)
  t_naive <- system.time(
    fit_naive <- do.call(multichain_MCMC_copula,
                         c(list(X = d$X, U1 = pobs(d$Y1), U2 = pobs(d$Y2)), cop_args)))[["elapsed"]]

  # two-step with BART margins
  fit_bart <- do.call(two_step_MCMC_copula,
                      c(list(Y1 = d$Y1, Y2 = d$Y2, X = d$X, X_marg = d$X,
                             margin_method = "bart", pobs_type = "rank"), cop_args))

  # two-step with GAM margins
  fit_gam <- do.call(two_step_MCMC_copula,
                     c(list(Y1 = d$Y1, Y2 = d$Y2, X = d$X, X_marg = d$X,
                            margin_method = "gam", pobs_type = "rank"), cop_args))

  # joint / full model
  joint_args <- cop_args[setdiff(names(cop_args), "n.cores")]
  fit_joint <- do.call(multichain_MCMC_copula_joint,
                       c(list(X = d$X, X_marg = d$X, Y1 = d$Y1, Y2 = d$Y2,
                              n.tree.marg = n.tree.marg, n.cores = n.cores,
                              prop_sigma_marg = 0.1), joint_args))

  fits  <- list(naive = fit_naive, `2step-BART` = fit_bart, `2step-GAM` = fit_gam, joint = fit_joint)
  times <- c(t_naive, fit_bart$time[["total"]], fit_gam$time[["total"]], fit_joint$time)

  results[[r]] <- do.call(rbind, lapply(seq_along(fits), function(k){
    ts <- tau_summary(fits[[k]], d$X, cop_type, n.iter_par, n.chain_par, n.discard)
    data.frame(rep = r, method = names(fits)[k], t(metrics(ts, d$tau)), time_sec = times[k])
  }))
  print(results[[r]])

  save(fits, d, file = paste0("cond_margins_", cop_type, "_case", test_case, "_rep", r, ".Rdata"))
  rm(fits, fit_naive, fit_bart, fit_gam, fit_joint); gc()
}

res_all <- do.call(rbind, results)
summary_tab <- aggregate(cbind(RMSE, CI_length, CI_cov, time_sec) ~ method, data = res_all, FUN = mean)
print(summary_tab)
save(res_all, summary_tab, file = paste0("cond_margins_", cop_type, "_case", test_case, "_summary.Rdata"))
