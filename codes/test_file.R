if(T){
  lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944))
  n.iter_par <- 500
  n.chain_par <- 1
  test_case = 3
  for (i in test_case) {
    assign(paste0("t_mcmc_lb.def_unif_",i), MCMC_copula(#n.chain = n.chain_par,
                                                                   n.iter = n.iter_par,
                                                                   n.tree = 10,
                                                                   X = X_obs.norm,
                                                                   U1 = get(paste0("copula_uu_t_",i))[,1],
                                                                   U2 = get(paste0("copula_uu_t_",i))[,2],
                                                                   prior_list = lb.prior.def, 
                                                                   moves.prob = moves.prob_par, 
                                                                   starting.tree = NULL,
                                                                   cont.unif = cont.unif_par,
                                                                   include.split = incl.split_par,
                                                                   prop_mu = 0, prop_sigma = 1,
                                                                   theta_param_1 = 0, theta_param_2 = 1,
                                                                   prior_type = "N",
                                                                   cop_type = "t"))
  }
}

if(T){
  model.list.def <- list(get(paste0("t_mcmc_lb.def_unif_",test_case))#,
    # get(paste0("t_mcmc_lb.def_half_",test_case)),
    # get(paste0("t_mcmc_lb.def_jeff_",test_case)),
    # get(paste0("t_mcmc_lb.def_two_",test_case))
  )
  
  names(model.list.def) <- c(
    'LB - default - unif'#,
    # 'LB - default - half',
    # 'LB - default - jeff',
    # 'LB - default - two'
  )
  
  
  list_pred_lb_10 <- lapply(1:length(model.list.def$`LB - default - unif`$trees), \(idx) BART_calculate_pred(model.list.def$`LB - default - unif`$trees[[idx]], X_obs_pred.norm))
  
  pred_val = do.call(rbind,list_pred_lb_10)
  
  n.born.out.par <- 250
  n.thin <- 1
  
  pred_val_vec = as.vector(pred_val[(n.born.out.par+1):nrow(pred_val),])
  
  pred_obs = rep(X_obs_pred.norm, each = (n.chain_par * (n.iter_par - n.born.out.par)))
  
  theta_true = rep(param_t(get(paste0("tau_true_pred_",test_case))), each = (n.chain_par * (n.iter_par - n.born.out.par)))
  
  pred_cond <- data.frame("obs" = pred_obs)
  pred_cond$obs = pred_obs
  pred_cond$theta_true = theta_true
  pred_cond$y = link_t(pred_val_vec)
  
  pred_cond_thin = na.omit(pred_cond[c(rep(NA,(n.thin-1)), TRUE),])
  
  pred_cond_mod = pred_cond_thin %>%
    group_by(obs, theta_true) %>%
    summarise(theta_mean = mean(y), theta_q975 = quantile(y, .975), theta_q025 = quantile(y, .025)) 
  
  ggplot(pred_cond_mod) +
    geom_line(aes(obs, theta_mean)) +
    geom_line(aes(obs, theta_true), col = 2) +
    geom_line(aes(obs, theta_q975), col = 3) +
    geom_line(aes(obs, theta_q025), col = 3) +
    # facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('X') +
    ylab('estimated rho') +
    theme_classic()
  
  
  pred_cond_stat = pred_cond_mod %>%
    mutate(RMSE = mean((theta_true - theta_mean)^2)) %>%
    mutate(CI.length = mean(theta_q975 - theta_q025)) %>%
    mutate(CI.cov = mean((theta_true < theta_q975) & (theta_true > theta_q025))) %>%
    dplyr::select(c(RMSE, CI.length, CI.cov))
  
  pred_cond_summary = colMeans(pred_cond_stat[,-1])
  
  # nterm
  
  nt_lb_10.df <- nterm_BART(model.list.def$`LB - default - unif`)
  nt_lb_10.df_burn <- nt_lb_10.df[nt_lb_10.df$idx >= n.born.out.par,]
  
  pl_nl_10 <- ggplot(nt_lb_10.df_burn, aes(trees, nn)) + 
    geom_boxplot() + 
    ylab('') + 
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  pl_nl_10
  
  # depth
  
  depth_lb_10.df <- depth_BART(model.list.def$`LB - default - unif`)
  depth_lb_10.df_burn <- depth_lb_10.df[depth_lb_10.df$idx >= n.born.out.par,]
  
  ggplot(depth_lb_10.df_burn, aes(trees, nn)) + 
    geom_boxplot() + 
    ylab('Depth') +
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  
  n.born.out.par <- 0
  # extract depth, number of terminal nodes and loglik of all the trees
  depth.df <- apply_fun_models(fun_ = get_depth,
                               mcmc.list = model.list.def,
                               born.out.pc = 0, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes,
                               mcmc.list = model.list.def,
                               born.out.pc = 0, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  like.df <- apply_fun_models(fun_ = \(x) cart_log_lik_copula(tree_top = x, U1 = get(paste0("copula_uu_t_",i))[,1],
                                                              U2 = get(paste0("copula_uu_t_",i))[,2],
                                                              X = X_obs.norm, 
                                                              log_like_fun = function(rho, u, v) loglik_t(exp(rho), u, v)), #((2*sigmoid(rho)-1), u, v) 
                              mcmc.list = model.list.def, 
                              born.out.pc = 0, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  acc.df <- apply_fun_models_res(fun_ = \(x) (as.factor(x) == TRUE),
                                 mcmc.list = model.list.def,
                                 born.out.pc = 0, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  
  depth_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(depth.df, i, 0, 1)))
  depth_conv_red <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(depth.df, i, n.born.out.par, n.thin)))
  
  nterm_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(nterm.df, i, 0, 1)))
  nterm_conv_red <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(nterm.df, i, n.born.out.par, n.thin)))
  
  like_conv_all <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(like.df, i, 0, 1)))
  like_conv_red <- do.call(rbind,lapply(names(model.list.def), function(i)conv_diag(like.df, i, n.born.out.par, n.thin)))
  
  xtable(cbind(depth_conv_all, nterm_conv_all, like_conv_all))
  xtable(cbind(depth_conv_red, nterm_conv_red, like_conv_red))
  
  depth.df_thin <- na.omit((depth.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  nterm.df_thin <- na.omit((nterm.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  like.df_thin <- na.omit((like.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  acc.df_thin <- na.omit((acc.df[(n.born.out.par+1:n.iter_par),])[c(rep(NA,(n.thin-1)), TRUE),])
  
  df.sum.def <- data.frame(tree = nterm.df$x,
                           panel.name = nterm.df$panel.name,
                           loglik = like.df$y,
                           nterm = nterm.df$y)
  
  df.sum.def_thin <- data.frame(tree = nterm.df_thin$x,
                                panel.name = nterm.df_thin$panel.name,
                                loglik = like.df_thin$y,
                                nterm = nterm.df_thin$y,
                                depth = depth.df_thin$y,
                                acc_rate = acc.df_thin$y)
  # df.sum.def_thin <- df.sum.def
  hist.nl <- ggplot(nterm.df_thin) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab(~n[L]) +
    ylab('PMF') +
    theme_classic()
  
  hist.depth <- ggplot(depth.df_thin) +
    geom_histogram(aes(x = y, y = after_stat(density)),
                   binwidth = 1, color = 'black', fill = 'white') +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Depth') +
    ylab('PMF') +
    theme_classic()
  
  trace.nl <- ggplot(nterm.df_thin) +
    geom_line(aes(x/n.thin, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab(~n[L]) +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,((n.iter_par-n.born.out.par)/n.thin),by = (n.iter_par-n.born.out.par)/(5*n.thin)))
  
  trace.depth <- ggplot(depth.df_thin) +
    geom_line(aes(x/n.thin, y)) +
    facet_wrap(facets = ~panel.name, ncol = 2) +
    xlab('Iteration') +
    ylab('Depth') +
    theme_classic() +
    scale_x_continuous(breaks = seq(0,((n.iter_par-n.born.out.par)/n.thin),by = (n.iter_par-n.born.out.par)/(5*n.thin)))
  
  trace.loglik <- ggplot(df.sum.def_thin, aes(tree/n.thin, loglik)) + 
    geom_line() + 
    facet_wrap(facets = ~panel.name) + 
    theme_classic() + 
    xlab('Iteration') + 
    ylab('Log likelihood')
  
  # plots
  
  hist.nl
  hist.depth
  trace.nl
  trace.depth
  trace.loglik
  
  
  tree_stat_summary = df.sum.def_thin %>%
    group_by(panel.name)%>%
    summarise_at(c("nterm","depth","acc_rate"),mean)
    
  xtable(cbind(tree_stat_summary, pred_cond_summary[,-1]) , digits = 4)
}
