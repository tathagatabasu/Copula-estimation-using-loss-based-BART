if(T){
  lb.prior.def <- list(fun = joint.prior.new.tree, param = c(1.5618883, 0.6293944)) # c(1.5618883, 0.6293944)
  n.iter_par <- 1000
  for (i in 1) {
    assign(paste0("clayton_mcmc_lb.def_unif_",i), MCMC_copula(#n.chain = n.chain_par,
                                                                   n.iter = n.iter_par,
                                                                   X = X_obs.norm,
                                                                   U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                                   U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                                   prior_list = lb.prior.def, 
                                                                   moves.prob = moves.prob_par, 
                                                                   starting.tree = NULL,
                                                                   cont.unif = cont.unif_par,
                                                                   include.split = incl.split_par,
                                                                   prop_mu = 1, prop_sigma = .25,
                                                                   theta_param_1 = 1, theta_param_2 = 1,
                                                                   prior_type = "LN",
                                                                   cop_type = "clayton"))
  }
  
  # for (i in 1) {
  #   assign(paste0("t_mcmc_lb.def_jeff_",i), multichain_MCMC_copula(n.chain = n.chain_par,
  #                                                                  n.iter = n.iter_par,
  #                                                                  X = X_obs.norm,
  #                                                                  U1 = get(paste0("copula_uu_t_",i))[,1],
  #                                                                  U2 = get(paste0("copula_uu_t_",i))[,2],
  #                                                                  mu = 0, cop_type = "t", 
  #                                                                  sigma = .005, alpha_val = 0, beta_val = 0, 
  #                                                                  log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
  #                                                                  prior_list = lb.prior.def, 
  #                                                                  moves.prob = moves.prob_par, 
  #                                                                  starting.tree = NULL,
  #                                                                  cont.unif = cont.unif_par,
  #                                                                  include.split = incl.split_par))
  # }
  # 
  # 
  # for (i in 1) {
  #   assign(paste0("t_mcmc_lb.def_half_",i), multichain_MCMC_copula(n.chain = n.chain_par,
  #                                                                  n.iter = n.iter_par,
  #                                                                  X = X_obs.norm,
  #                                                                  U1 = get(paste0("copula_uu_t_",i))[,1],
  #                                                                  U2 = get(paste0("copula_uu_t_",i))[,2],
  #                                                                  mu = 0, cop_type = "t", 
  #                                                                  sigma = .005, alpha_val = .5, beta_val = .5, 
  #                                                                  log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
  #                                                                  prior_list = lb.prior.def, 
  #                                                                  moves.prob = moves.prob_par, 
  #                                                                  starting.tree = NULL,
  #                                                                  cont.unif = cont.unif_par,
  #                                                                  include.split = incl.split_par))
  # }
  # 
  # for (i in 1) {
  #   assign(paste0("t_mcmc_lb.def_two_",i), multichain_MCMC_copula(n.chain = n.chain_par,
  #                                                                 n.iter = n.iter_par,
  #                                                                 X = X_obs.norm,
  #                                                                 U1 = get(paste0("copula_uu_t_",i))[,1],
  #                                                                 U2 = get(paste0("copula_uu_t_",i))[,2],
  #                                                                 mu = 0, cop_type = "t", 
  #                                                                 sigma = .005, alpha_val = 2, beta_val = 2, 
  #                                                                 log_nor_mu = 0, log_nor_sigma = 1, prior_type = "B",
  #                                                                 prior_list = lb.prior.def, 
  #                                                                 moves.prob = moves.prob_par, 
  #                                                                 starting.tree = NULL,
  #                                                                 cont.unif = cont.unif_par,
  #                                                                 include.split = incl.split_par))
  # }
  
}

# results

if(T){
  test_case = 1
  
  n.born.out.par <- 0
  
  model.list.def <- list(get(paste0("clayton_mcmc_lb.def_unif_",test_case))#,
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
  
  # extract depth, number of terminal nodes and loglik of all the trees
  depth.df <- apply_fun_models(fun_ = get_depth,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  nterm.df <- apply_fun_models(fun_ = get_num_terminal_nodes,
                               mcmc.list = model.list.def,
                               born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  like.df <- apply_fun_models(fun_ = \(x) cart_log_lik_copula(tree_top = x, U1 = get(paste0("copula_uu_clayton_",i))[,1],
                                                              U2 = get(paste0("copula_uu_clayton_",i))[,2],
                                                              X = X_obs.norm, 
                                                              log_like_fun = loglik_clayton), 
                              mcmc.list = model.list.def, 
                              born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  acc.df <- apply_fun_models_res(fun_ = \(x) (as.factor(x) == TRUE),
                                 mcmc.list = model.list.def,
                                 born.out.pc = n.born.out.par, n.chain = n.chain_par, sample.pc = n.iter_par)
  
  n.born.out.par <- 500
  n.thin <- 1
  
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
  
  copula_uu_t_pred <- BiCopSim(N = nrow(pred_cond_mod), family = 2, par = pred_cond_mod$theta_mean, par2 = 3)
  
  pred_cond_mod$U1 = copula_uu_t_pred[,1]
  pred_cond_mod$U2 = copula_uu_t_pred[,2]
  
  pred_cond_mod_p1 = pred_cond_mod %>% filter(panel.name == "LB - default - unif")
  pred_cond_mod_p2 = pred_cond_mod %>% filter(panel.name == "LB - default - half")
  pred_cond_mod_p3 = pred_cond_mod %>% filter(panel.name == "LB - default - jeff")
  pred_cond_mod_p4 = pred_cond_mod %>% filter(panel.name == "LB - default - two")
  
  # 3d histogram
  
  hist_true <- hist2d(get(paste0("copula_uu_t_",i))[,1], get(paste0("copula_uu_t_",i))[,2], nbins = c(10,10), show = FALSE)
  
  hist_pred1 <- hist2d(pred_cond_mod_p1$U1, pred_cond_mod_p1$U2, nbins = c(10,10), show = FALSE)
  hist_pred2 <- hist2d(pred_cond_mod_p2$U1, pred_cond_mod_p2$U2, nbins = c(10,10), show = FALSE)
  hist_pred3 <- hist2d(pred_cond_mod_p3$U1, pred_cond_mod_p3$U2, nbins = c(10,10), show = FALSE)
  hist_pred4 <- hist2d(pred_cond_mod_p4$U1, pred_cond_mod_p4$U2, nbins = c(10,10), show = FALSE)
  
  par(mar = c(1,1,1,1), mfrow = c(2,3))
  
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
    x = hist_pred1$x,
    y = hist_pred1$y,
    z = hist_pred1$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = unique(pred_cond_mod_p1$panel.name)
  )
  
  hist3D(
    x = hist_pred2$x,
    y = hist_pred2$y,
    z = hist_pred2$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = unique(pred_cond_mod_p2$panel.name)
  )
  
  hist3D(
    x = hist_pred3$x,
    y = hist_pred3$y,
    z = hist_pred3$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = unique(pred_cond_mod_p3$panel.name)
  )
  
  hist3D(
    x = hist_pred4$x,
    y = hist_pred4$y,
    z = hist_pred4$counts,
    colvar = NULL,  # disables color mapping
    col = "lightblue",  # solid monochrome color
    border = "black",
    theta = -30, scale = FALSE, expand = 0.02, bty = "g", phi = 45,    # shading gives 3D effect
    lighting = TRUE,
    ltheta = 120, ticktype = "detailed",
    xlab = "", ylab = "", zlab = "",
    main = unique(pred_cond_mod_p4$panel.name)
  )
  
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
}
