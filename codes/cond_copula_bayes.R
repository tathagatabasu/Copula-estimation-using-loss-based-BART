# codes and packages
source('code/import_functions.R')
source('mclapply.R')
# source('MCMC_BART_copula.R')
source('test_MCMC_copula_mult_tree.R')
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

################################################################################
# data generation
################################################################################
set.seed(1e5)

load("analysis_sim_dat_new.RData")

if(T){
  n <- 500
  X_obs <- matrix(runif(n), ncol = 1)
  
  # normalise predictors 
  X_obs.norm <- as.data.frame(apply(X_obs, 2, \(x) (x - min(x))/(max(x) - min(x))))
  X_obs.norm <- as.matrix(X_obs.norm)
  rownames(X_obs.norm) <- 1:nrow(X_obs)
  
  # tau with tree structure
  # initialize list
  tree_ex <- list()
  # one split on the left branch
  tree_ex$left <- list(left = list(left = NULL,
                                   right = NULL),
                       right = list(left = NULL,
                                    right = NULL))
  # no split on the right branch
  tree_ex$right <-  list(left = NULL,
                         right = NULL)
  
  # assign node index
  tree_ex <- assign_node_idx(tree_ex)
  # plot the tree
  get_tree_plot.idx(tree_ex)
  
  # assign first splitting rule
  tree_ex$cond <- list(x.idx = 1, x.val = 0.6)
  # assign second splitting rule
  tree_ex$left$cond <- list(x.idx = 1, x.val = 0.25)
  # set mu_3
  tree_ex$left$left$cond <- NULL
  tree_ex$left$left$value <- .5
  # set mu_4
  tree_ex$left$right$cond <- NULL
  tree_ex$left$right$value <- .7
  # set mu_5
  tree_ex$right$cond <- NULL
  tree_ex$right$value <- .3
  # plot the tree
  get_tree_plot(tree_ex)
  tau_true_1 <- sample_CART(tree_ex, X_obs, sigma_ = 0.001) 
  tau_true_1 <- matrix(tau_true_1, ncol = 1)
  
  # monotone
  tau_true_2 <- 0.3 + 0.2 * sin(3*X_obs) + 0.3*X_obs^2
  # convex
  tau_true_3 <- 0.5 + 0.3 * sin(3*X_obs)
  # non-convex
  tau_true_4 <- 0.6 - 0.3 * sin(2*X_obs) + 0.2 * sin(4*X_obs) + 0.3 * X_obs^2
  
  plot(X_obs, tau_true_1, xlab = "Observations", ylab = "tau")
  plot(X_obs, tau_true_2, xlab = "Observations", ylab = "tau")
  plot(X_obs, tau_true_3, xlab = "Observations", ylab = "tau")
  plot(X_obs, tau_true_4, xlab = "Observations", ylab = "tau")
  
  # gauss
  
  for (i in 1:4) {
    assign(paste0("copula_uu_gauss_",i), BiCopSim(n, family = 1, par = param_gauss(get(paste0("tau_true_",i)))))
  }
  
  # t
  
  for (i in 1:4) {
    assign(paste0("copula_uu_t_",i), BiCopSim(n, family = 2, par = param_gauss(get(paste0("tau_true_",i))), par2 = 3))
  }
  
  # gumbel 
  
  for (i in 1:4) {
    assign(paste0("copula_uu_gumbel_",i), BiCopSim(n, family = 4, par = param_gumbel(get(paste0("tau_true_",i)))))
  }
  
  # clayton
  
  for (i in 1:4) {
    assign(paste0("copula_uu_clayton_",i), BiCopSim(n, family = 3, par = param_clayton(get(paste0("tau_true_",i)))))
  }
  
  # dataset for prediction
  
  X_obs_pred = matrix(runif(n), ncol = 1)
  
  # normalise predictors 
  X_obs_pred.norm <- as.data.frame(apply(X_obs_pred, 2, \(x) (x - min(x))/(max(x) - min(x))))
  X_obs_pred.norm <- as.matrix(X_obs_pred.norm)
  rownames(X_obs_pred.norm) <- 1:nrow(X_obs_pred)
  
  # tau with tree structure
  tau_true_pred_1 <- sample_CART(tree_ex, X_obs_pred, sigma_ = 0.001) 
  tau_true_pred_1 <- matrix(tau_true_pred_1, ncol = 1)
  
  # monotone
  tau_true_pred_2 <- 0.3 + 0.2 * sin(3*X_obs_pred) + 0.3*X_obs_pred^2
  # convex
  tau_true_pred_3 <- 0.5 + 0.3 * sin(3*X_obs_pred)
  # non-convex
  tau_true_pred_4 <- 0.6 - 0.3 * sin(2*X_obs_pred) + 0.2 * sin(4*X_obs_pred) + 0.3 * X_obs_pred^2
  
  # mcmc params
  n.chain_par <- 1
  n.iter_par <- 6000
  n.born.out.par <- 1000
  n.thin <- 1
  incl.split_par <- TRUE
  cont.unif_par <- TRUE
  moves.prob_par <- c(0.4, 0.4, 0.1, 0.1)
  lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) 
  
}

################################################################################
# gaussian
################################################################################
if(T){
  for (i in 1:4) {
    assign(paste0("gauss_mcmc_lb.def_single_",i), MCMC_copula(n.iter = n.iter_par,
                                                                  n.tree = 1,
                                                                  X = X_obs.norm,
                                                                  U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                                  U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                                  prior_list = lb.prior.def, 
                                                                  moves.prob = moves.prob_par, 
                                                                  starting.tree = NULL,
                                                                  cont.unif = cont.unif_par,
                                                                  include.split = incl.split_par,
                                                                  prop_mu = 0, prop_sigma = 1,
                                                                  theta_param_1 = 0, theta_param_2 = 1,
                                                                  prior_type = "N",
                                                                  cop_type = "gauss"))
    
    cat('done case', i, '\n')
  }
  
}

# results

if(F){
  test_case = 3
  
  model.list.def <- list(
    get(paste0("gauss_mcmc_lb.def_single_",test_case))
  )
  
  names(model.list.def) <- c(
    'LB - default - TBeta(1.1,1.1)'
  )
  
  # extract depth, number of terminal nodes and loglik of all the trees
  depth.df <- apply_fun_models(fun_ = get_depth,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  depth.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  nterm.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  like.df <- apply_fun_models(fun_ = \(x) cart_log_lik_copula(tree_top = x, U1 = get(paste0("copula_uu_gauss_",i))[,1],
                                                              U2 = get(paste0("copula_uu_gauss_",i))[,2],
                                                              X = X_obs.norm, 
                                                              log_like_fun = loglik_gaussian), 
                              mcmc.list = model.list.def, 
                              born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  like.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  acc.df <- apply_fun_models_res(fun_ = \(x) (as.factor(x) == TRUE),
                                 mcmc.list = model.list.def,
                                 born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  acc.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  depth_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(depth.df, i, 0, 10)))
  
  nterm_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(nterm.df, i, 0, 10)))
  
  like_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(like.df, i, 0, 10)))
  
  xtable(cbind(depth_conv_all, nterm_conv_all, like_conv_all))
  
  depth.df_thin <- na.omit((depth.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  nterm.df_thin <- na.omit((nterm.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  like.df_thin <- na.omit((like.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  acc.df_thin <- na.omit((acc.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  
  df.sum.def <- data.frame(tree = nterm.df$x,
                           panel.name = nterm.df$panel.name,
                           loglik = like.df$y,
                           nterm = nterm.df$y,
                           depth = depth.df$y,
                           acc_rate = acc.df$y)
  
  df.sum.def_thin <- data.frame(tree = nterm.df_thin$x,
                                panel.name = nterm.df_thin$panel.name,
                                loglik = like.df_thin$y,
                                nterm = nterm.df_thin$y,
                                depth = depth.df_thin$y,
                                acc_rate = acc.df_thin$y)
  
  hist.nl <- ggplot(nterm.df) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab(~n[L]) +
    ylab('PMF') +
    theme_classic()
  
  hist.depth <- ggplot(depth.df) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Depth') +
    ylab('PMF') +
    theme_classic()
  
  trace.nl <- ggplot(nterm.df) +
    geom_line(aes(x/n.thin, y, group = factor(n.chain), colour =factor(n.chain))) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab(~n[L]) +
    theme_classic()  + 
    theme(legend.position = "none")+
    scale_x_continuous(breaks = seq(0,((n.chain_par * (n.iter_par-n.born.out.par))/n.thin),by = (n.chain_par * (n.iter_par-n.born.out.par))/(5*n.thin)))
  
  trace.depth <- ggplot(depth.df) +
    geom_line(aes(x/n.thin, y, group = factor(n.chain), colour =factor(n.chain))) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab('Depth') +
    theme_classic()  + 
    theme(legend.position = "none")+
    scale_x_continuous(breaks = seq(0,((n.chain_par * (n.iter_par-n.born.out.par))/n.thin),by = (n.chain_par * (n.iter_par-n.born.out.par))/(5*n.thin)))
  
  trace.loglik <- ggplot(like.df, aes(x/n.thin, y, group = factor(n.chain), colour =factor(n.chain))) + 
    geom_line() + 
    facet_wrap(facets = ~panel.name) + 
    theme_classic() + 
    theme(legend.position = "none")+
    xlab('Iteration') + 
    ylab('Log likelihood') + 
    scale_x_continuous(breaks = seq(0,((n.chain_par * (n.iter_par-n.born.out.par))/n.thin),by = (n.chain_par * (n.iter_par-n.born.out.par))/(5*n.thin)))
  
  # plots
  hist.nl
  hist.depth
  trace.nl
  trace.depth
  trace.loglik
  
  # prediction
  
  pred_cond = lapply(1:nrow(X_obs_pred.norm), function(i)apply_fun_models(fun_ = function(x)get_value_tree(x, X_obs_pred.norm[i,,drop = FALSE]),
                                                                          mcmc.list = model.list.def,
                                                                          born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par))
  
  
  pred_cond = do.call(rbind,pred_cond)
  
  pred_cond$obs = as.vector(apply(X_obs_pred.norm, 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  pred_cond$theta_true = as.vector(apply(param_gauss(get(paste0("tau_true_pred_",test_case))), 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,(n.thin-1)), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(panel.name, obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025))  
  
  ggplot(pred_cond_mod) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated rho') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod %>%
    group_by(panel.name)%>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = pred_cond_stat %>%
    group_by(panel.name)%>%
    summarise_all(mean)
  
  tree_stat_summary = df.sum.def_thin %>%
    group_by(panel.name)%>%
    summarise_at(c("nterm","depth","acc_rate"),mean)
  
  xtable(cbind(tree_stat_summary, pred_cond_summary[,-1]) , digits = 4)
  
  if(F){
    copula_uu_gauss_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 1, par = pred_cond_mod$theta_mean)
    
    pred_cond_mod$U1 = copula_uu_gauss_pred[,1]
    pred_cond_mod$U2 = copula_uu_gauss_pred[,2]
    
    # Plot 3D histogram
    
    hist_true <- hist2d(get(paste0("copula_uu_gauss_",i))[,1], get(paste0("copula_uu_gauss_",i))[,2], nbins = c(10,10), show = FALSE)
    
    hist_pred <- hist2d(pred_cond_mod$U1, pred_cond_mod$U2, nbins = c(10,10), show = FALSE)
    
    par(mar = c(1,1,1,1), mfrow = c(1,2))
    
    hist3D(
      x = hist_true$x,
      y = hist_true$y,
      z = hist_true$counts,
      colvar = NULL,  # disables color mapping
      col = "lightblue",  # solid monochrome color
      border = "black",
      theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
      lighting = TRUE,
      ltheta = 120, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "",
      main = "Observed copula"
    )
    
    hist3D(
      x = hist_pred$x,
      y = hist_pred$y,
      z = hist_pred$counts,
      colvar = NULL,  # disables color mapping
      col = "lightblue",  # solid monochrome color
      border = "black",
      theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
      lighting = TRUE,
      ltheta = 120, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "",
      main = unique(pred_cond_mod$panel.name)
    )
  }
}

################################################################################
# t
################################################################################
if(F){
  
  for (i in 1:4) {
    assign(paste0("t_mcmc_lb.def_single_",i), MCMC_copula(
                                                              n.iter = n.iter_par,
                                                              n.tree = 1,
                                                              X = X_obs.norm,
                                                              U1 = get(paste0("copula_uu_t_",i))[,1],
                                                              U2 = get(paste0("copula_uu_t_",i))[,2],
                                                              prior_list = lb.prior.def, 
                                                              moves.prob = moves.prob_par, 
                                                              starting.tree = NULL,
                                                              cont.unif = cont.unif_par,
                                                              include.split = incl.split_par,
                                                              prop_mu = 0, prop_sigma = .2,
                                                              theta_param_1 = 1.1, theta_param_2 = 1.1,
                                                              prior_type = "B",
                                                              cop_type = "t"))
  }
  
}

# results

if(F){
  test_case = 1
  
  model.list.def <- list(
    get(paste0("t_mcmc_lb.def_single_",test_case))
  )
  
  names(model.list.def) <- c(
    'LB - default - TBeta(1.1,1.1)'
  )
  
  # extract depth, number of terminal nodes and loglik of all the trees
  depth.df <- apply_fun_models(fun_ = get_depth,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  depth.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  nterm.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  like.df <- apply_fun_models(fun_ = \(x) cart_log_lik_copula(tree_top = x, U1 = get(paste0("copula_uu_t_",i))[,1],
                                                              U2 = get(paste0("copula_uu_t_",i))[,2],
                                                              X = X_obs.norm, 
                                                              log_like_fun = loglik_t), 
                              mcmc.list = model.list.def, 
                              born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  like.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  acc.df <- apply_fun_models_res(fun_ = \(x) (as.factor(x) == TRUE),
                                 mcmc.list = model.list.def,
                                 born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  acc.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  depth_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(depth.df, i, 0, n.thin)))
  
  nterm_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(nterm.df, i, 0, n.thin)))
  
  like_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(like.df, i, 0, n.thin)))
  
  xtable(cbind(depth_conv_all, nterm_conv_all, like_conv_all))
  
  depth.df_thin <- na.omit((depth.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  nterm.df_thin <- na.omit((nterm.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  like.df_thin <- na.omit((like.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  acc.df_thin <- na.omit((acc.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  
  df.sum.def <- data.frame(tree = nterm.df$x,
                           panel.name = nterm.df$panel.name,
                           loglik = like.df$y,
                           nterm = nterm.df$y,
                           depth = depth.df$y,
                           acc_rate = acc.df$y)
  
  df.sum.def_thin <- data.frame(tree = nterm.df_thin$x,
                                panel.name = nterm.df_thin$panel.name,
                                loglik = like.df_thin$y,
                                nterm = nterm.df_thin$y,
                                depth = depth.df_thin$y,
                                acc_rate = acc.df_thin$y)
  
  hist.nl <- ggplot(nterm.df) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab(~n[L]) +
    ylab('PMF') +
    theme_classic()
  
  hist.depth <- ggplot(depth.df) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Depth') +
    ylab('PMF') +
    theme_classic()
  
  trace.nl <- ggplot(nterm.df) +
    geom_line(aes(x/n.thin, y, group = factor(n.chain), colour =factor(n.chain))) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab(~n[L]) +
    theme_classic()  + 
    theme(legend.position = "none")+
    scale_x_continuous(breaks = seq(0,((n.chain_par * (n.iter_par-n.born.out.par))/n.thin),by = (n.chain_par * (n.iter_par-n.born.out.par))/(5*n.thin)))
  
  trace.depth <- ggplot(depth.df) +
    geom_line(aes(x/n.thin, y, group = factor(n.chain), colour =factor(n.chain))) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab('Depth') +
    theme_classic()  + 
    theme(legend.position = "none")+
    scale_x_continuous(breaks = seq(0,((n.chain_par * (n.iter_par-n.born.out.par))/n.thin),by = (n.chain_par * (n.iter_par-n.born.out.par))/(5*n.thin)))
  
  trace.loglik <- ggplot(like.df, aes(x/n.thin, y, group = factor(n.chain), colour =factor(n.chain))) + 
    geom_line() + 
    facet_wrap(facets = ~panel.name) + 
    theme_classic() + 
    theme(legend.position = "none")+
    xlab('Iteration') + 
    ylab('Log likelihood') + 
    scale_x_continuous(breaks = seq(0,((n.chain_par * (n.iter_par-n.born.out.par))/n.thin),by = (n.chain_par * (n.iter_par-n.born.out.par))/(5*n.thin)))
  
  # plots
  hist.nl
  hist.depth
  trace.nl
  trace.depth
  trace.loglik
  
  # prediction
  
  pred_cond = lapply(1:nrow(X_obs_pred.norm), function(i)apply_fun_models(fun_ = function(x)get_value_tree(x, X_obs_pred.norm[i,,drop = FALSE]),
                                                                          mcmc.list = model.list.def,
                                                                          born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par))
  
  
  pred_cond = do.call(rbind,pred_cond)
  
  pred_cond$obs = as.vector(apply(X_obs_pred.norm, 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  pred_cond$theta_true = as.vector(apply(param_t(get(paste0("tau_true_pred_",test_case))), 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,(n.thin-1)), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(panel.name, obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025))  
  
  ggplot(pred_cond_mod) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated rho') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod %>%
    group_by(panel.name)%>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = pred_cond_stat %>%
    group_by(panel.name)%>%
    summarise_all(mean)
  
  tree_stat_summary = df.sum.def_thin %>%
    group_by(panel.name)%>%
    summarise_at(c("nterm","depth","acc_rate"),mean)
  
  xtable(cbind(tree_stat_summary, pred_cond_summary[,-1]) , digits = 4)
  
  if(F){
    copula_uu_t_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 2, par = pred_cond_mod$theta_mean)
    
    pred_cond_mod$U1 = copula_uu_t_pred[,1]
    pred_cond_mod$U2 = copula_uu_t_pred[,2]
    
    # Plot 3D histogram
    
    hist_true <- hist2d(get(paste0("copula_uu_t_",i))[,1], get(paste0("copula_uu_t_",i))[,2], nbins = c(10,10), show = FALSE)
    
    hist_pred <- hist2d(pred_cond_mod$U1, pred_cond_mod$U2, nbins = c(10,10), show = FALSE)
    
    par(mar = c(1,1,1,1), mfrow = c(1,2))
    
    hist3D(
      x = hist_true$x,
      y = hist_true$y,
      z = hist_true$counts,
      colvar = NULL,  # disables color mapping
      col = "lightblue",  # solid monochrome color
      border = "black",
      theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
      lighting = TRUE,
      ltheta = 120, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "",
      main = "Observed copula"
    )
    
    hist3D(
      x = hist_pred$x,
      y = hist_pred$y,
      z = hist_pred$counts,
      colvar = NULL,  # disables color mapping
      col = "lightblue",  # solid monochrome color
      border = "black",
      theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
      lighting = TRUE,
      ltheta = 120, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "",
      main = unique(pred_cond_mod$panel.name)
    )
  }
}

################################################################################
# Gumbel
################################################################################
if(F){
  
  for (i in 1:4) {
    assign(paste0("gumbel_mcmc_lb.def_single_",i), MCMC_copula(
                                                                   n.iter = n.iter_par,
                                                                   n.tree = 1,
                                                                   X = X_obs.norm,
                                                                   U1 = get(paste0("copula_uu_gumbel_",i))[,1],
                                                                   U2 = get(paste0("copula_uu_gumbel_",i))[,2],
                                                                   prior_list = lb.prior.def, 
                                                                   moves.prob = moves.prob_par, 
                                                                   starting.tree = NULL,
                                                                   cont.unif = cont.unif_par,
                                                                   include.split = incl.split_par,
                                                                   prop_mu = 1, prop_sigma = 1,
                                                                   theta_param_1 = 0, theta_param_2 = 2,
                                                                   prior_type = "LN",
                                                                   cop_type = "gumbel"))
  }
  
}

# results

if(F){
  test_case = 4
  
  model.list.def <- list(
    get(paste0("gumbel_mcmc_lb.def_single_",test_case))
  )
  
  names(model.list.def) <- c(
    'LB - default - TBeta(1.1,1.1)'
  )
  
  # extract depth, number of terminal nodes and loglik of all the trees
  depth.df <- apply_fun_models(fun_ = get_depth,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  depth.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  nterm.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  like.df <- apply_fun_models(fun_ = \(x) cart_log_lik_copula(tree_top = x, U1 = get(paste0("copula_uu_gumbel_",i))[,1],
                                                              U2 = get(paste0("copula_uu_gumbel_",i))[,2],
                                                              X = X_obs.norm, 
                                                              log_like_fun = loglik_gumbel), 
                              mcmc.list = model.list.def, 
                              born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  like.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  acc.df <- apply_fun_models_res(fun_ = \(x) (as.factor(x) == TRUE),
                                 mcmc.list = model.list.def,
                                 born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  acc.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  depth_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(depth.df, i, 0, 10)))
  
  nterm_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(nterm.df, i, 0, 10)))
  
  like_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(like.df, i, 0, 10)))
  
  xtable(cbind(depth_conv_all, nterm_conv_all, like_conv_all))
  
  depth.df_thin <- na.omit((depth.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  nterm.df_thin <- na.omit((nterm.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  like.df_thin <- na.omit((like.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  acc.df_thin <- na.omit((acc.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  
  df.sum.def <- data.frame(tree = nterm.df$x,
                           panel.name = nterm.df$panel.name,
                           loglik = like.df$y,
                           nterm = nterm.df$y,
                           depth = depth.df$y,
                           acc_rate = acc.df$y)
  
  df.sum.def_thin <- data.frame(tree = nterm.df_thin$x,
                                panel.name = nterm.df_thin$panel.name,
                                loglik = like.df_thin$y,
                                nterm = nterm.df_thin$y,
                                depth = depth.df_thin$y,
                                acc_rate = acc.df_thin$y)
  
  hist.nl <- ggplot(nterm.df) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab(~n[L]) +
    ylab('PMF') +
    theme_classic()
  
  hist.depth <- ggplot(depth.df) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Depth') +
    ylab('PMF') +
    theme_classic()
  
  trace.nl <- ggplot(nterm.df) +
    geom_line(aes(x/n.thin, y, group = factor(n.chain), colour =factor(n.chain))) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab(~n[L]) +
    theme_classic()  + 
    theme(legend.position = "none")+
    scale_x_continuous(breaks = seq(0,((n.chain_par * (n.iter_par-n.born.out.par))/n.thin),by = (n.chain_par * (n.iter_par-n.born.out.par))/(5*n.thin)))
  
  trace.depth <- ggplot(depth.df) +
    geom_line(aes(x/n.thin, y, group = factor(n.chain), colour =factor(n.chain))) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab('Depth') +
    theme_classic()  + 
    theme(legend.position = "none")+
    scale_x_continuous(breaks = seq(0,((n.chain_par * (n.iter_par-n.born.out.par))/n.thin),by = (n.chain_par * (n.iter_par-n.born.out.par))/(5*n.thin)))
  
  trace.loglik <- ggplot(like.df, aes(x/n.thin, y, group = factor(n.chain), colour =factor(n.chain))) + 
    geom_line() + 
    facet_wrap(facets = ~panel.name) + 
    theme_classic() + 
    theme(legend.position = "none")+
    xlab('Iteration') + 
    ylab('Log likelihood') + 
    scale_x_continuous(breaks = seq(0,((n.chain_par * (n.iter_par-n.born.out.par))/n.thin),by = (n.chain_par * (n.iter_par-n.born.out.par))/(5*n.thin)))
  
  # plots
  hist.nl
  hist.depth
  trace.nl
  trace.depth
  trace.loglik
  
  # prediction
  
  pred_cond = lapply(1:nrow(X_obs_pred.norm), function(i)apply_fun_models(fun_ = function(x)get_value_tree(x, X_obs_pred.norm[i,,drop = FALSE]),
                                                                          mcmc.list = model.list.def,
                                                                          born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par))
  
  
  pred_cond = do.call(rbind,pred_cond)
  
  pred_cond$obs = as.vector(apply(X_obs_pred.norm, 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  pred_cond$theta_true = as.vector(apply(param_gumbel(get(paste0("tau_true_pred_",test_case))), 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,(n.thin-1)), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(panel.name, obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025))  
  
  ggplot(pred_cond_mod) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated rho') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod %>%
    group_by(panel.name)%>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = pred_cond_stat %>%
    group_by(panel.name)%>%
    summarise_all(mean)
  
  tree_stat_summary = df.sum.def_thin %>%
    group_by(panel.name)%>%
    summarise_at(c("nterm","depth","acc_rate"),mean)
  
  xtable(cbind(tree_stat_summary, pred_cond_summary[,-1]) , digits = 4)
  
  if(F){
    copula_uu_gumbel_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 1, par = pred_cond_mod$theta_mean)
    
    pred_cond_mod$U1 = copula_uu_gumbel_pred[,1]
    pred_cond_mod$U2 = copula_uu_gumbel_pred[,2]
    
    # Plot 3D histogram
    
    hist_true <- hist2d(get(paste0("copula_uu_gumbel_",i))[,1], get(paste0("copula_uu_gumbel_",i))[,2], nbins = c(10,10), show = FALSE)
    
    hist_pred <- hist2d(pred_cond_mod$U1, pred_cond_mod$U2, nbins = c(10,10), show = FALSE)
    
    par(mar = c(1,1,1,1), mfrow = c(1,2))
    
    hist3D(
      x = hist_true$x,
      y = hist_true$y,
      z = hist_true$counts,
      colvar = NULL,  # disables color mapping
      col = "lightblue",  # solid monochrome color
      border = "black",
      theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
      lighting = TRUE,
      ltheta = 120, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "",
      main = "Observed copula"
    )
    
    hist3D(
      x = hist_pred$x,
      y = hist_pred$y,
      z = hist_pred$counts,
      colvar = NULL,  # disables color mapping
      col = "lightblue",  # solid monochrome color
      border = "black",
      theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
      lighting = TRUE,
      ltheta = 120, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "",
      main = unique(pred_cond_mod$panel.name)
    )
  }
}

################################################################################
# clayton
################################################################################
if(F){
  
  for (i in 3) {
    assign(paste0("clayton_mcmc_lb.def_single_",i), MCMC_copula(
                                                                    n.iter = n.iter_par,
                                                                    n.tree = 1,
                                                                    X = X_obs.norm,
                                                                    U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                    U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                    prior_list = lb.prior.def, 
                                                                    moves.prob = moves.prob_par, 
                                                                    starting.tree = NULL,
                                                                    cont.unif = cont.unif_par,
                                                                    include.split = incl.split_par,
                                                                    prop_mu = 0, prop_sigma = .2,
                                                                    theta_param_1 = 0, theta_param_2 = 1,
                                                                    prior_type = "LN",
                                                                    cop_type = "clayton"))
  }
  
  # all prop_mu = 0, prop_sigma = .2,
  # 3 prop_mu = 0, prop_sigma = 3,
  
}

# results

if(F){
  test_case <- 3
  
  model.list.def <- list(
    get(paste0("clayton_mcmc_lb.def_single_",test_case))
  )
  
  names(model.list.def) <- c(
    'LB - default - TBeta(1.1,1.1)'
  )
  
  # extract depth, number of terminal nodes and loglik of all the trees
  depth.df <- apply_fun_models(fun_ = get_depth,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  depth.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  nterm.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  like.df <- apply_fun_models(fun_ = \(x) cart_log_lik_copula(tree_top = x, U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                              U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                              X = X_obs.norm, 
                                                              log_like_fun = loglik_clayton), 
                              mcmc.list = model.list.def, 
                              born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  like.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  acc.df <- apply_fun_models_res(fun_ = \(x) (as.factor(x) == TRUE),
                                 mcmc.list = model.list.def,
                                 born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  acc.df$n.chain <- rep(1:n.chain_par, each = (n.iter_par- n.born.out.par))
  
  depth_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(depth.df, i, 0, 10)))
  
  nterm_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(nterm.df, i, 0, 10)))
  
  like_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(like.df, i, 0, 10)))
  
  xtable(cbind(depth_conv_all, nterm_conv_all, like_conv_all))
  
  depth.df_thin <- na.omit((depth.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  nterm.df_thin <- na.omit((nterm.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  like.df_thin <- na.omit((like.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  acc.df_thin <- na.omit((acc.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  
  df.sum.def <- data.frame(tree = nterm.df$x,
                           panel.name = nterm.df$panel.name,
                           loglik = like.df$y,
                           nterm = nterm.df$y,
                           depth = depth.df$y,
                           acc_rate = acc.df$y)
  
  df.sum.def_thin <- data.frame(tree = nterm.df_thin$x,
                                panel.name = nterm.df_thin$panel.name,
                                loglik = like.df_thin$y,
                                nterm = nterm.df_thin$y,
                                depth = depth.df_thin$y,
                                acc_rate = acc.df_thin$y)
  
  hist.nl <- ggplot(nterm.df) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab(~n[L]) +
    ylab('PMF') +
    theme_classic()
  
  hist.depth <- ggplot(depth.df) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Depth') +
    ylab('PMF') +
    theme_classic()
  
  trace.nl <- ggplot(nterm.df) +
    geom_line(aes(x/n.thin, y, group = factor(n.chain), colour =factor(n.chain))) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab(~n[L]) +
    theme_classic()  + 
    theme(legend.position = "none")+
    scale_x_continuous(breaks = seq(0,((n.chain_par * (n.iter_par-n.born.out.par))/n.thin),by = (n.chain_par * (n.iter_par-n.born.out.par))/(5*n.thin)))
  
  trace.depth <- ggplot(depth.df) +
    geom_line(aes(x/n.thin, y, group = factor(n.chain), colour =factor(n.chain))) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab('Depth') +
    theme_classic()  + 
    theme(legend.position = "none")+
    scale_x_continuous(breaks = seq(0,((n.chain_par * (n.iter_par-n.born.out.par))/n.thin),by = (n.chain_par * (n.iter_par-n.born.out.par))/(5*n.thin)))
  
  trace.loglik <- ggplot(like.df, aes(x/n.thin, y, group = factor(n.chain), colour =factor(n.chain))) + 
    geom_line() + 
    facet_wrap(facets = ~panel.name) + 
    theme_classic() + 
    theme(legend.position = "none")+
    xlab('Iteration') + 
    ylab('Log likelihood') + 
    scale_x_continuous(breaks = seq(0,((n.chain_par * (n.iter_par-n.born.out.par))/n.thin),by = (n.chain_par * (n.iter_par-n.born.out.par))/(5*n.thin)))
  
  # plots
  hist.nl
  hist.depth
  trace.nl
  trace.depth
  trace.loglik
  
  # prediction
  
  pred_cond = lapply(1:nrow(X_obs_pred.norm), function(i)apply_fun_models(fun_ = function(x)get_value_tree(x, X_obs_pred.norm[i,,drop = FALSE]),
                                                                          mcmc.list = model.list.def,
                                                                          born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par))
  
  
  pred_cond = do.call(rbind,pred_cond)
  
  pred_cond$obs = as.vector(apply(X_obs_pred.norm, 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  pred_cond$theta_true = as.vector(apply(param_clayton(get(paste0("tau_true_pred_",test_case))), 1, function(x)rep(x, (length(model.list.def) * n.chain_par * (n.iter_par - n.born.out.par)))))
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,(n.thin-1)), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(panel.name, obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025))  
  
  ggplot(pred_cond_mod) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated rho') +
    theme_classic()
  
  pred_cond_stat = pred_cond_mod %>%
    group_by(panel.name)%>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = pred_cond_stat %>%
    group_by(panel.name)%>%
    summarise_all(mean)
  
  tree_stat_summary = df.sum.def_thin %>%
    group_by(panel.name)%>%
    summarise_at(c("nterm","depth","acc_rate"),mean)
  
  xtable(cbind(tree_stat_summary, pred_cond_summary[,-1]) , digits = 4)
  
  if(F){
    copula_uu_clayton_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 1, par = pred_cond_mod$theta_mean)
    
    pred_cond_mod$U1 = copula_uu_clayton_pred[,1]
    pred_cond_mod$U2 = copula_uu_clayton_pred[,2]
    
    # Plot 3D histogram
    
    hist_true <- hist2d(get(paste0("copula_uu_clayton_",i))[,1], get(paste0("copula_uu_clayton_",i))[,2], nbins = c(10,10), show = FALSE)
    
    hist_pred <- hist2d(pred_cond_mod$U1, pred_cond_mod$U2, nbins = c(10,10), show = FALSE)
    
    par(mar = c(1,1,1,1), mfrow = c(1,2))
    
    hist3D(
      x = hist_true$x,
      y = hist_true$y,
      z = hist_true$counts,
      colvar = NULL,  # disables color mapping
      col = "lightblue",  # solid monochrome color
      border = "black",
      theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
      lighting = TRUE,
      ltheta = 120, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "",
      main = "Observed copula"
    )
    
    hist3D(
      x = hist_pred$x,
      y = hist_pred$y,
      z = hist_pred$counts,
      colvar = NULL,  # disables color mapping
      col = "lightblue",  # solid monochrome color
      border = "black",
      theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
      lighting = TRUE,
      ltheta = 120, ticktype = "detailed",
      xlab = "", ylab = "", zlab = "",
      main = unique(pred_cond_mod$panel.name)
    )
  }
}
