################################################################################
# Conditional copula BART with conditional margins
#
# This file adds two ways of handling the conditional margins F_{j|X}:
#
#   (A) TWO-STEP   : estimate F_{1|X}, F_{2|X} first (BART or GAM), build the
#                    pseudo-observations U_ij = F_hat_{j|X}(y_ij | x_i) and run
#                    the ORIGINAL sampler multichain_MCMC_copula() on them.
#
#   (B) JOINT/FULL : Y_j | x ~ N(f_j(x), sigma_j^2), with f_j a sum of trees,
#                    and (U1, U2) | x ~ C(. | theta(x)), U_j = Phi((y_j - f_j)/sigma_j).
#                    All three tree ensembles (copula, margin 1, margin 2) and
#                    sigma_1, sigma_2 are updated in one MCMC under the full
#                    joint likelihood
#                       sum_i log c(u_i1, u_i2 | theta(x_i))
#                     + sum_i log phi(y_i1; f_1(x_i), s_1) + log phi(y_i2; f_2(x_i), s_2).
#
# Only MCMC drivers are added here. All tree functions (moves, priors,
# transition probabilities, BART_single_step, tree_step_copula, ...) are
# re-used exactly as they are in MCMC_BART_copula.R and the tree_*.R files.
# The marginal trees use the same RJ-MCMC step: BART_single_step() is simply
# given a log-likelihood closure for the margin instead of the copula one.
#
# Requires: source('MCMC_BART_copula.R'); source('import_functions.R')
# Extra packages: dbarts (two-step BART margins), mgcv (two-step GAM margins)
################################################################################

################################################################################
# Small helpers
################################################################################

# copula log-likelihood on the sum-of-trees (link) scale -- same as in MCMC_copula
get_cop_loglik_fun <- function(cop_type){
  switch(cop_type,
         gauss   = function(rho, u, v) loglik_gauss(link_gauss(rho), u, v),
         t       = function(rho, u, v) loglik_t(link_t(rho), u, v, df = 3),
         clayton = function(rho, u, v) loglik_clayton(link_clayton(rho), u, v),
         gumbel  = function(rho, u, v) loglik_gumbel(link_gumbel(rho), u, v),
         frank   = function(rho, u, v) loglik_frank(link_frank(rho), u, v),
         stop("Unsupported copula type"))
}

# leaf-value prior -- same as in MCMC_copula
get_log_prior_fun <- function(prior_type){
  switch(prior_type,
         B  = function(rho, p1, p2) log_prior_tbeta(rho, bound = 1, p1, p2),
         LN = function(rho, p1, p2) dlnorm(rho, meanlog = p1, sdlog = p2, log = TRUE),
         G  = function(rho, p1, p2) dgamma(rho, shape = p1, rate = p2, log = TRUE),
         N  = function(rho, p1, p2) dnorm(rho, mean = p1, sd = p2, log = TRUE),
         stop("Unsupported prior type"))
}

# map sum-of-trees value to conditional Kendall's tau
eta_to_tau <- function(eta, cop_type){
  switch(cop_type,
         gauss   = BiCopPar2Tau(1, link_gauss(eta)),
         t       = BiCopPar2Tau(2, link_t(eta), par2 = 3),
         clayton = BiCopPar2Tau(3, link_clayton(eta)),
         gumbel  = BiCopPar2Tau(4, link_gumbel(eta)),
         frank   = BiCopPar2Tau(5, link_frank(eta)),
         stop("Unsupported copula type"))
}

# keep pseudo-observations strictly inside (0,1) for BiCopPDF
clamp_u <- function(u, eps = 1e-10) pmin(pmax(u, eps), 1 - eps)

# trees use rownames(X) as observation ids
prep_X <- function(X){
  X <- as.matrix(X)
  rownames(X) <- seq_len(nrow(X))
  X
}

################################################################################
# (A) TWO-STEP: conditional margins
################################################################################

# GAM formula with one smooth per numeric column (linear if few unique values)
gam_rhs <- function(df){
  terms <- vapply(names(df), function(v){
    nu <- length(unique(df[[v]]))
    if(is.numeric(df[[v]]) && nu >= 5) paste0("s(", v, ", k = ", min(10, nu - 1), ")") else v
  }, "")
  paste(terms, collapse = " + ")
}

# Estimate F_{j|X}(y | x) for one response.
#   method   : "bart" (dbarts) or "gam" (mgcv)
#   scale    : TRUE  -> location-scale model, Y = mu(x) + sigma(x) * e
#              FALSE -> location model,       Y = mu(x) + sigma * e
#   pobs_type: "rank"       -> U = rank(e_hat) / (n + 1)   (error distribution left free)
#              "parametric" -> U = Phi(e_hat)               (Gaussian errors)
fit_conditional_margin <- function(y, X, method = c("bart", "gam"),
                                   scale = TRUE,
                                   pobs_type = c("rank", "parametric"),
                                   bart_ntree = 200, bart_ndpost = 1000, bart_nskip = 500,
                                   seed = 1){
  method <- match.arg(method)
  pobs_type <- match.arg(pobs_type)
  X <- as.matrix(X)
  n <- length(y)

  if(method == "bart"){
    if(!requireNamespace("dbarts", quietly = TRUE)) stop("Please install 'dbarts'")
    set.seed(seed)
    fit_mu <- dbarts::bart(x.train = X, y.train = y, ntree = bart_ntree,
                           ndpost = bart_ndpost, nskip = bart_nskip,
                           verbose = FALSE, keeptrees = FALSE)
    mu <- colMeans(fit_mu$yhat.train)
    if(scale){
      # BART on log squared residuals; E[log chi^2_1] = -1.2704
      r2 <- log((y - mu)^2 + 1e-8)
      fit_sd <- dbarts::bart(x.train = X, y.train = r2, ntree = bart_ntree,
                             ndpost = bart_ndpost, nskip = bart_nskip,
                             verbose = FALSE, keeptrees = FALSE)
      sigma <- exp((colMeans(fit_sd$yhat.train) + 1.2704) / 2)
      # re-scale so that the standardised residuals have unit variance
      sigma <- sigma * sd((y - mu) / sigma)
    } else {
      sigma <- rep(mean(fit_mu$sigma), n)
    }
    fit <- list(mu = fit_mu, sd = if(scale) fit_sd else NULL)

  } else {
    if(!requireNamespace("mgcv", quietly = TRUE)) stop("Please install 'mgcv'")
    df <- as.data.frame(X)
    colnames(df) <- paste0("x", seq_len(ncol(df)))
    df$y <- y
    rhs <- gam_rhs(df[, setdiff(names(df), "y"), drop = FALSE])
    if(scale){
      fit <- mgcv::gam(list(as.formula(paste("y ~", rhs)), as.formula(paste("~", rhs))),
                       family = mgcv::gaulss(), data = df, method = "REML")
      pr <- predict(fit, type = "response")
      mu <- pr[, 1]
      sigma <- 1 / pr[, 2]       # gaulss: 2nd column on response scale is 1/sd
    } else {
      fit <- mgcv::gam(as.formula(paste("y ~", rhs)), data = df, method = "REML")
      mu <- as.numeric(fitted(fit))
      sigma <- rep(sqrt(fit$sig2), n)
    }
  }

  e_hat <- (y - mu) / sigma
  U <- if(pobs_type == "rank") rank(e_hat, ties.method = "average") / (n + 1) else pnorm(e_hat)

  list(U = clamp_u(U), mu = mu, sigma = sigma, resid = e_hat,
       method = method, scale = scale, pobs_type = pobs_type, fit = fit)
}

# Two-step estimator: margins first, then the ORIGINAL copula sampler.
#   Y1, Y2  : responses (original scale)
#   X       : covariates for theta(x)
#   X_marg  : covariates for the margins (default: same as X)
#   ...     : every argument of multichain_MCMC_copula() except X, U1, U2
two_step_MCMC_copula <- function(Y1, Y2, X, X_marg = X,
                                 margin_method = c("bart", "gam"),
                                 margin_scale = NULL,
                                 pobs_type = c("rank", "parametric"),
                                 margin_args = list(),
                                 ...){
  margin_method <- match.arg(margin_method)
  pobs_type <- match.arg(pobs_type)
  X <- prep_X(X)
  # default: location-scale GAM (gaulss); location-only BART, since the
  # log-squared-residual BART for sigma(x) is noisy for small n
  if(is.null(margin_scale)) margin_scale <- (margin_method == "gam")

  t_marg <- system.time({
    m1 <- do.call(fit_conditional_margin,
                  c(list(y = Y1, X = X_marg, method = margin_method,
                         scale = margin_scale, pobs_type = pobs_type), margin_args))
    m2 <- do.call(fit_conditional_margin,
                  c(list(y = Y2, X = X_marg, method = margin_method,
                         scale = margin_scale, pobs_type = pobs_type), margin_args))
  })[["elapsed"]]

  t_cop <- system.time({
    fit <- multichain_MCMC_copula(X = X, U1 = m1$U, U2 = m2$U, ...)
  })[["elapsed"]]

  fit$margins <- list(m1 = m1[c("mu", "sigma", "U", "method", "scale", "pobs_type")],
                      m2 = m2[c("mu", "sigma", "U", "method", "scale", "pobs_type")])
  fit$U1 <- m1$U
  fit$U2 <- m2$U
  fit$time <- c(margins = t_marg, copula = t_cop, total = t_marg + t_cop)
  fit
}

################################################################################
# (B) JOINT / FULL MODEL
################################################################################

# One sweep over the trees of one margin.
# Calls the unchanged BART_single_step(); only the log-likelihood closure is new.
#   which = 1 : margin of Y1, closure is log phi(y1) + log c(Phi(.), U2)
#   which = 2 : margin of Y2, closure is log phi(y2) + log c(U1, Phi(.))
update_margin_trees <- function(which, y_s, U_other, eta, s_j,
                                tree_list, pred_mat, sig2_vec,
                                X_marg, cop_ll, log_prior_m,
                                prior_list, moves.prob, cont.unif, include.split,
                                prop_sigma_m, var_param_1_m, var_param_2_m,
                                var_list = NULL, adapt = FALSE){

  n.tree.m <- length(tree_list)
  idx.vec <- seq_len(n.tree.m)

  if(which == 1){
    ll_marg <- function(f, y, u_other){
      u <- clamp_u(pnorm((y - f) / s_j))
      sum(dnorm(y, f, s_j, log = TRUE)) + cop_ll(eta, u, u_other)
    }
  } else {
    ll_marg <- function(f, y, u_other){
      u <- clamp_u(pnorm((y - f) / s_j))
      sum(dnorm(y, f, s_j, log = TRUE)) + cop_ll(eta, u_other, u)
    }
  }

  df_list <- vector("list", n.tree.m)
  for(t in idx.vec){
    res_f <- rowSums(pred_mat[, idx.vec != t, drop = FALSE])
    tv <- if(adapt) var_list[[t]]$cov else NULL

    step <- BART_single_step(X = X_marg,
                             res_theta = res_f,
                             U1 = y_s, U2 = U_other,
                             prop_mu = 0, prop_sigma = prop_sigma_m,
                             prior_list = prior_list,
                             moves.prob = moves.prob,
                             starting.tree = tree_list[[t]],
                             diag = FALSE,
                             cont.unif = cont.unif,
                             include.split = include.split,
                             obs.per.term = 1,
                             empty.count.lim = 10,
                             log_like_fun = ll_marg,
                             log_prior_fun = log_prior_m,
                             theta_param_1 = 0, theta_param_2 = sqrt(sig2_vec[t]),
                             var_param_1 = var_param_1_m, var_param_2 = var_param_2_m,
                             theta_var_list = tv,
                             idx.tree = t)

    tree_list[[t]] <- step$tree
    sig2_vec[t] <- as.numeric(step$df_res$sig2)
    pred_mat[, t] <- get_value_tree(step$tree, X_marg)
    df_list[[t]] <- step$df_res
  }
  list(trees = tree_list, pred = pred_mat, sig2 = sig2_vec, df_res = df_list)
}

# Random-walk MH for log(sigma_j) under the full likelihood,
# prior sigma_j^2 ~ InvGamma(a, b) on the standardised scale.
update_margin_sigma <- function(which, y_s, f, s_j, U_other, eta, cop_ll,
                                a_sig, b_sig, step_sd){
  lp <- function(s){
    u <- clamp_u(pnorm((y_s - f) / s))
    cl <- if(which == 1) cop_ll(eta, u, U_other) else cop_ll(eta, U_other, u)
    sum(dnorm(y_s, f, s, log = TRUE)) + cl +
      (-a_sig * log(s^2) - b_sig / s^2)          # IG prior on s^2 + log-Jacobian
  }
  s_new <- exp(log(s_j) + rnorm(1, 0, step_sd))
  logr <- lp(s_new) - lp(s_j)
  if(is.finite(logr) && log(runif(1)) < logr) list(s = s_new, acc = TRUE) else list(s = s_j, acc = FALSE)
}

MCMC_copula_joint <- function(n.iter, n.burn,
                              n.tree = 10, n.tree.marg = 10,
                              X, X_marg = X,
                              Y1, Y2,
                              prior_list,
                              moves.prob = NULL,
                              cont.unif = TRUE,
                              include.split,
                              prop_mu, prop_sigma,
                              theta_param_1, theta_param_2,
                              var_param_1, var_param_2,
                              prior_type, cop_type, adapt = TRUE,
                              # margin-specific settings
                              prop_sigma_marg = 0.1,
                              var_param_1_marg = 1, var_param_2_marg = 2,
                              a_sig = 1, b_sig = 1,
                              sigma_step = 0.05,
                              init_margins = NULL,
                              store.marg.trees = TRUE,
                              verbose = 0){

  X <- prep_X(X)
  X_marg <- prep_X(X_marg)
  n <- nrow(X)

  cop_ll <- get_cop_loglik_fun(cop_type)
  log_prior_fun <- get_log_prior_fun(prior_type)
  log_prior_m <- get_log_prior_fun("N")

  # standardise responses; everything below is on the standardised scale
  m_y <- c(mean(Y1), mean(Y2)); s_y <- c(sd(Y1), sd(Y2))
  y1 <- (Y1 - m_y[1]) / s_y[1]
  y2 <- (Y2 - m_y[2]) / s_y[2]

  # ---- initial state -------------------------------------------------------
  idx.tree.vec <- seq_len(n.tree)
  idx.m.vec <- seq_len(n.tree.marg)

  last.tree.list <- vector("list", n.tree)
  res_theta_list <- matrix(0, nrow = n, ncol = n.tree)
  sig_2 <- rinvgamma(n.tree, shape = var_param_1, scale = var_param_2)

  trees_m1 <- vector("list", n.tree.marg); trees_m2 <- trees_m1
  pred_m1 <- matrix(0, n, n.tree.marg);    pred_m2 <- pred_m1
  sig2_m1 <- rinvgamma(n.tree.marg, shape = var_param_1_marg, scale = var_param_2_marg)
  sig2_m2 <- rinvgamma(n.tree.marg, shape = var_param_1_marg, scale = var_param_2_marg)
  s1 <- 1; s2 <- 1

  # optional starting values for sigma_1, sigma_2 (original scale),
  # e.g. init_margins = list(sigma1 = 1, sigma2 = 1)
  if(!is.null(init_margins)){
    if(!is.null(init_margins$sigma1)) s1 <- mean(init_margins$sigma1) / s_y[1]
    if(!is.null(init_margins$sigma2)) s2 <- mean(init_margins$sigma2) / s_y[2]
  }

  # ---- storage -------------------------------------------------------------
  trees.for.iter <- vector("list", n.iter)
  df_res.for.iter <- vector("list", n.iter)
  trees_m1.for.iter <- if(store.marg.trees) vector("list", n.iter) else NULL
  trees_m2.for.iter <- if(store.marg.trees) vector("list", n.iter) else NULL
  theta_pred_list <- array(NA_real_, dim = c(n, n.tree, n.burn))
  m1_pred_burn <- array(NA_real_, dim = c(n, n.tree.marg, n.burn))
  m2_pred_burn <- array(NA_real_, dim = c(n, n.tree.marg, n.burn))
  mu1_out <- matrix(NA_real_, n.iter, n); mu2_out <- mu1_out
  sigma_out <- matrix(NA_real_, n.iter, 2, dimnames = list(NULL, c("sigma1", "sigma2")))
  loglik_out <- matrix(NA_real_, n.iter, 3, dimnames = list(NULL, c("copula", "margin1", "margin2")))
  acc_sig <- c(0, 0); step_sd <- c(sigma_step, sigma_step)

  var_list <- NULL; var_m1 <- NULL; var_m2 <- NULL
  st_time <- Sys.time()

  for(idx.iter in seq_len(n.iter)){

    use_adapt <- adapt && (idx.iter > n.burn)

    # ---- current pseudo-observations ---------------------------------------
    f1 <- rowSums(pred_m1); f2 <- rowSums(pred_m2)
    U1 <- clamp_u(pnorm((y1 - f1) / s1))
    U2 <- clamp_u(pnorm((y2 - f2) / s2))

    # ---- 1. copula trees (identical to MCMC_copula) -----------------------
    df_res.list <- vector("list", n.tree)
    for(idx.tree in idx.tree.vec){
      theta_var_list_in <- if(use_adapt) var_list[[idx.tree]]$cov else NULL
      res_theta <- rowSums(res_theta_list[, idx.tree.vec != idx.tree, drop = FALSE])

      new_tree_single <- BART_single_step(X = X,
                                          res_theta = res_theta,
                                          U1 = U1, U2 = U2,
                                          prop_mu = prop_mu, prop_sigma = prop_sigma,
                                          prior_list = prior_list,
                                          moves.prob = moves.prob,
                                          starting.tree = last.tree.list[[idx.tree]],
                                          diag = FALSE,
                                          cont.unif = cont.unif,
                                          include.split = include.split,
                                          obs.per.term = 1,
                                          empty.count.lim = 10,
                                          log_like_fun = cop_ll,
                                          log_prior_fun = log_prior_fun,
                                          theta_param_1 = theta_param_1, theta_param_2 = sqrt(sig_2[idx.tree]),
                                          var_param_1 = var_param_1, var_param_2 = var_param_2,
                                          theta_var_list = theta_var_list_in,
                                          idx.tree = idx.tree)

      last.tree.list[[idx.tree]] <- new_tree_single$tree
      sig_2[idx.tree] <- as.numeric(new_tree_single$df_res$sig2)
      df_res.list[[idx.tree]] <- new_tree_single$df_res
      res_theta_list[, idx.tree] <- get_value_tree(new_tree_single$tree, X)
    }
    eta <- rowSums(res_theta_list)

    # ---- 2. margin 1: trees, then sigma_1 ----------------------------------
    up1 <- update_margin_trees(1, y1, U2, eta, s1, trees_m1, pred_m1, sig2_m1,
                               X_marg, cop_ll, log_prior_m, prior_list, moves.prob,
                               cont.unif, include.split, prop_sigma_marg,
                               var_param_1_marg, var_param_2_marg,
                               var_list = var_m1, adapt = use_adapt)
    trees_m1 <- up1$trees; pred_m1 <- up1$pred; sig2_m1 <- up1$sig2
    f1 <- rowSums(pred_m1)

    ss1 <- update_margin_sigma(1, y1, f1, s1, U2, eta, cop_ll, a_sig, b_sig, step_sd[1])
    s1 <- ss1$s; acc_sig[1] <- acc_sig[1] + ss1$acc
    U1 <- clamp_u(pnorm((y1 - f1) / s1))

    # ---- 3. margin 2: trees, then sigma_2 ----------------------------------
    up2 <- update_margin_trees(2, y2, U1, eta, s2, trees_m2, pred_m2, sig2_m2,
                               X_marg, cop_ll, log_prior_m, prior_list, moves.prob,
                               cont.unif, include.split, prop_sigma_marg,
                               var_param_1_marg, var_param_2_marg,
                               var_list = var_m2, adapt = use_adapt)
    trees_m2 <- up2$trees; pred_m2 <- up2$pred; sig2_m2 <- up2$sig2
    f2 <- rowSums(pred_m2)

    ss2 <- update_margin_sigma(2, y2, f2, s2, U1, eta, cop_ll, a_sig, b_sig, step_sd[2])
    s2 <- ss2$s; acc_sig[2] <- acc_sig[2] + ss2$acc
    U2 <- clamp_u(pnorm((y2 - f2) / s2))

    # tune sigma step size during the burn-in only (target ~0.44)
    if(idx.iter <= n.burn && idx.iter %% 50 == 0){
      rate <- acc_sig / 50
      step_sd <- step_sd * exp(rate - 0.44)
      acc_sig <- c(0, 0)
    }

    # ---- adaptive proposal covariances (same scheme as MCMC_copula) --------
    if(idx.iter <= n.burn){
      theta_pred_list[, , idx.iter] <- res_theta_list
      m1_pred_burn[, , idx.iter] <- pred_m1
      m2_pred_burn[, , idx.iter] <- pred_m2
      if(idx.iter == n.burn){
        var_list <- lapply(idx.tree.vec, function(x) var_adapt(theta_pred_list[, x, 1:n.burn]))
        var_m1 <- lapply(idx.m.vec, function(x) var_adapt(m1_pred_burn[, x, 1:n.burn]))
        var_m2 <- lapply(idx.m.vec, function(x) var_adapt(m2_pred_burn[, x, 1:n.burn]))
        rm(theta_pred_list, m1_pred_burn, m2_pred_burn)
      }
    } else {
      var_list <- lapply(idx.tree.vec, function(x) iter_var(var_list[[x]]$cov, var_list[[x]]$mean, res_theta_list[, x], idx.iter))
      var_m1 <- lapply(idx.m.vec, function(x) iter_var(var_m1[[x]]$cov, var_m1[[x]]$mean, pred_m1[, x], idx.iter))
      var_m2 <- lapply(idx.m.vec, function(x) iter_var(var_m2[[x]]$cov, var_m2[[x]]$mean, pred_m2[, x], idx.iter))
    }

    # ---- store --------------------------------------------------------------
    trees.for.iter[[idx.iter]] <- last.tree.list
    df_res.for.iter[[idx.iter]] <- df_res.list
    if(store.marg.trees){
      trees_m1.for.iter[[idx.iter]] <- trees_m1
      trees_m2.for.iter[[idx.iter]] <- trees_m2
    }
    mu1_out[idx.iter, ] <- m_y[1] + s_y[1] * f1
    mu2_out[idx.iter, ] <- m_y[2] + s_y[2] * f2
    sigma_out[idx.iter, ] <- c(s1 * s_y[1], s2 * s_y[2])
    loglik_out[idx.iter, ] <- c(cop_ll(eta, U1, U2),
                                sum(dnorm(y1, f1, s1, log = TRUE)) - n * log(s_y[1]),
                                sum(dnorm(y2, f2, s2, log = TRUE)) - n * log(s_y[2]))

    if(verbose > 0 && idx.iter %% verbose == 0)
      cat("iter", idx.iter, "| loglik", round(sum(loglik_out[idx.iter, ]), 2),
          "| sigma", round(sigma_out[idx.iter, ], 3), "\n")
  }

  time_diff <- Sys.time() - st_time
  print(time_diff)

  list(trees = trees.for.iter,          # copula trees (same format as MCMC_copula)
       df_res = df_res.for.iter,
       trees_m1 = trees_m1.for.iter,    # marginal trees, standardised scale
       trees_m2 = trees_m2.for.iter,
       mu1 = mu1_out, mu2 = mu2_out,    # n.iter x n, original scale
       sigma = sigma_out,               # n.iter x 2, original scale
       loglik = loglik_out,
       y_center = m_y, y_scale = s_y)
}

multichain_MCMC_copula_joint <- function(n.iter = 6000,
                                         n.burn = 1000,
                                         n.chain = 1,
                                         n.tree = 10,
                                         n.tree.marg = 10,
                                         seed.list = 1:n.chain,
                                         n.cores = 1,
                                         ...){
  st_time <- Sys.time()
  run_chain <- function(x){
    set.seed(seed.list[x])
    MCMC_copula_joint(n.iter = n.iter, n.burn = n.burn,
                      n.tree = n.tree, n.tree.marg = n.tree.marg, ...)
  }
  if(n.cores == 1){
    chain.list <- lapply(1:n.chain, run_chain)
  } else {
    chain.list <- mclapply_manual(1:n.chain, run_chain, mc.cores = n.cores, seed.list = seed.list)
  }

  out <- list(trees = Reduce(c, lapply(chain.list, `[[`, "trees")),
              df_res = Reduce(c, lapply(chain.list, `[[`, "df_res")),
              trees_m1 = Reduce(c, lapply(chain.list, `[[`, "trees_m1")),
              trees_m2 = Reduce(c, lapply(chain.list, `[[`, "trees_m2")),
              mu1 = do.call(rbind, lapply(chain.list, `[[`, "mu1")),
              mu2 = do.call(rbind, lapply(chain.list, `[[`, "mu2")),
              sigma = do.call(rbind, lapply(chain.list, `[[`, "sigma")),
              loglik = do.call(rbind, lapply(chain.list, `[[`, "loglik")),
              chain = rep(seq_len(n.chain), each = n.iter),
              y_center = chain.list[[1]]$y_center,
              y_scale = chain.list[[1]]$y_scale)
  out$time <- as.numeric(difftime(Sys.time(), st_time, units = "secs"))
  print(Sys.time() - st_time)
  out
}

# Predict the margins of a joint fit at new covariate values (original scale)
predict_margins_joint <- function(fit, X_new, iter = seq_along(fit$trees_m1)){
  X_new <- prep_X(X_new)
  if(is.null(fit$trees_m1)) stop("Fit was run with store.marg.trees = FALSE")
  mu1 <- t(vapply(iter, function(k) BART_calculate_pred(fit$trees_m1[[k]], X_new), numeric(nrow(X_new))))
  mu2 <- t(vapply(iter, function(k) BART_calculate_pred(fit$trees_m2[[k]], X_new), numeric(nrow(X_new))))
  list(mu1 = fit$y_center[1] + fit$y_scale[1] * mu1,
       mu2 = fit$y_center[2] + fit$y_scale[2] * mu2,
       sigma = fit$sigma[iter, , drop = FALSE])
}

################################################################################
# Posterior summaries of conditional Kendall's tau (works for all three fits)
################################################################################
tau_summary <- function(fit, X, cop_type, n.iter, n.chain, n.discard){
  X <- prep_X(X)
  keep <- unlist(lapply(seq_len(n.chain), function(c) (c - 1) * n.iter + ((n.discard + 1):n.iter)))
  eta <- t(vapply(keep, function(k) BART_calculate_pred(fit$trees[[k]], X), numeric(nrow(X))))
  tau <- matrix(eta_to_tau(as.vector(eta), cop_type), nrow = nrow(eta))
  data.frame(x = X[, 1],
             tau_mean = colMeans(tau),
             tau_q025 = apply(tau, 2, quantile, 0.025),
             tau_q975 = apply(tau, 2, quantile, 0.975))
}
