library(parallel)
library(mc2d)
library(pracma)
library(invgamma)
multichain_MCMC_copula <- function(n.iter, n.chain,n.tree = 10,
                                   X, U1, U2,
                                   prior_list,
                                   moves.prob = NULL, starting.tree = NULL,
                                   n.cores = min(n.chain, 5),
                                   cont.unif = TRUE,
                                   include.split,
                                   prop_mu, prop_sigma,
                                   theta_param_1, theta_param_2,
                                   prior_type, cop_type){
  chain.list <- mclapply(1:n.chain,
                         \(x) MCMC_copula(n.iter = n.iter, n.tree = n.tree,
                                          X = X,
                                          U1 = U1,
                                          U2 = U2,
                                          prior_list = prior_list,
                                          moves.prob = moves.prob,
                                          starting.tree = starting.tree,
                                          cont.unif = cont.unif,
                                          include.split = include.split,
                                          prop_mu = prop_mu, prop_sigma = prop_sigma,
                                          theta_param_1 = theta_param_1, theta_param_2 = theta_param_2,
                                          prior_type = prior_type, cop_type = cop_type),
                         mc.cores = n.cores)
  tree_list <- lapply(chain.list, function(x) x$trees)
  tree_list_comb <- Reduce(c, tree_list)

  df_list <- lapply(chain.list,
                    function(x) x$df_res)
  df.res <- Reduce(c, df_list)
  return(list(trees = tree_list_comb, df.res = df.res))
}

MCMC_copula <- function(n.iter, n.tree = 10, X, U1, U2, prior_list,
                        moves.prob = NULL, starting.tree = NULL, 
                        diag = FALSE, cont.unif = TRUE, include.split,
                        prop_mu, prop_sigma, theta_param_1, theta_param_2, prior_type, cop_type){
  
  ##############################################################################
  
  if(cop_type == "gauss"){
    log_like_fun <- function(rho, u, v) loglik_gauss(link_gauss(rho), u, v) 
  }else if(cop_type == "t"){
    log_like_fun <- function(rho, u, v) loglik_t(link_t(rho), u, v, df = 3)
  }else if(cop_type == "clayton"){
    log_like_fun <- function(rho, u, v) loglik_clayton(link_clayton(rho), u, v)
  }else if(cop_type == "gumbel"){
    log_like_fun <- function(rho, u, v) loglik_gumbel(link_gumbel(rho), u, v)
  }
  
  if(cop_type == "gauss"){
    fisher_fun_log <- FI_gauss 
  }else if(cop_type == "t"){
    fisher_fun_log <- function(rho, u, v) FI_t(rho, u, v, df = 3)
  }else if(cop_type == "clayton"){
    fisher_fun_log <- FI_clayton
  }else if(cop_type == "gumbel"){
    fisher_fun_log <- FI_gumbel
  }
  
  if(cop_type == "gauss"){
    fd_fun_log <- FD_log_gauss 
  }else if(cop_type == "t"){
    fd_fun_log <- function(rho, u, v) FD_log_t(rho, u, v, df = 3)
  }else if(cop_type == "clayton"){
    fd_fun_log <- FD_log_clayton
  }else if(cop_type == "gumbel"){
    fd_fun_log <- FD_log_gumbel
  }
  
  if(prior_type == "B"){
    log_prior_fun <- function(rho) log_prior_tbeta(rho, bound = 1, theta_param_1, theta_param_2)
  }else if(prior_type == "LN"){
    log_prior_fun <- function(rho) dlnorm(rho, meanlog = theta_param_1, sdlog = theta_param_2, log = TRUE)
  }else if(prior_type == "G"){
    log_prior_fun <- function(rho) dgamma(rho, shape = theta_param_1, rate = theta_param_2, log = TRUE)
  }else if(prior_type == "N"){
    log_prior_fun <- function(rho) dnorm(rho, mean = theta_param_1, sd = theta_param_2, log = TRUE)
  }
  ##
  
  # initialise output
  starting.tree.list <- lapply(1:n.tree, \(x) NULL)
  last.tree.list <- starting.tree.list
  df_res.list <- list()
  idx.tree.vec <- 1:n.tree
  res_theta_list <- lapply(1:n.tree, \(x) rep(0, length(X)))
  
  sig_2 <- rinvgamma(n.tree, shape = 0.01, rate = 0.01)
  
  # params for computation
  
  trees.for.iter <- list()
  df_res.for.iter <- list()
  
  st_time = Sys.time()
  
  if(n.tree>1){
    for(idx.iter in 1:n.iter){
      cat('Iteration: ', idx.iter, '\n')
      pred_at_obs <- Reduce('+', res_theta_list)
      
      for(idx.tree in idx.tree.vec){
        res_theta <- Reduce('+', res_theta_list[idx.tree.vec != idx.tree])
        
        cat('Tree: ', idx.tree, '\n')
        old_tree <- last.tree.list[[idx.tree]]
        
        new_tree_single <- BART_single_step(X = X, res_theta = res_theta, U1 = U1, U2 = U2, 
                                            prop_mu = prop_mu, prop_sigma = prop_sigma, prior_list = prior_list,
                                            moves.prob = moves.prob, starting.tree = old_tree, 
                                            diag = diag, cont.unif = cont.unif, include.split = include.split,
                                            obs.per.term = 1, empty.count.lim = 10, cop_type =  cop_type,
                                            log_like_fun = log_like_fun, 
                                            log_prior_fun = log_prior_fun,
                                            fisher_fun_log = fisher_fun_log,
                                            fd_fun_log = fd_fun_log,
                                            theta_param_1 = theta_param_1, theta_param_2 =  sqrt(sig_2[idx.tree]))
        
        new_tree <- new_tree_single$tree
        
        sig_2[idx.tree] <- as.numeric(new_tree_single$df_res$sig2)
        
        last.tree.list[[idx.tree]] <- new_tree
        df_res.list[[idx.tree]] <- new_tree_single$df_res
        theta_pred <- vapply(1:nrow(X), \(x) g.T(new_tree, X[x,]), 0)
        res_theta_list[[idx.tree]] <- theta_pred
        
      }
      trees.for.iter[[idx.iter]] <- last.tree.list
      df_res.for.iter[[idx.iter]] <- df_res.list
      
    }
  }else{
    for(idx.iter in 1:n.iter){
      cat('Iteration: ', idx.iter, '\n')
      pred_at_obs <- Reduce('+', res_theta_list)
      
      for(idx.tree in idx.tree.vec){
        res_theta <- rep(0,length(X))
        
        cat('Tree: ', idx.tree, '\n')
        old_tree <- last.tree.list[[idx.tree]]
        
        new_tree_single <- BART_single_step(X = X, res_theta = res_theta, U1 = U1, U2 = U2, 
                                            prop_mu = prop_mu, prop_sigma = prop_sigma, prior_list = prior_list,
                                            moves.prob = moves.prob, starting.tree = old_tree, 
                                            diag = diag, cont.unif = cont.unif, include.split = include.split,
                                            obs.per.term = 1, empty.count.lim = 10, cop_type =  cop_type,
                                            log_like_fun = log_like_fun, 
                                            log_prior_fun = log_prior_fun,
                                            fisher_fun_log = fisher_fun_log,
                                            fd_fun_log = fd_fun_log,
                                            theta_param_1 = theta_param_1, theta_param_2 =  sqrt(sig_2[idx.tree]))
        
        new_tree <- new_tree_single$tree
        
        sig_2[idx.tree] <- as.numeric(new_tree_single$df_res$sig2)
        
        last.tree.list[[idx.tree]] <- new_tree
        df_res.list[[idx.tree]] <- new_tree_single$df_res
        theta_pred <- vapply(1:nrow(X), \(x) g.T(new_tree, X[x,]), 0)
        res_theta_list[[idx.tree]] <- theta_pred
        
      }
      trees.for.iter[[idx.iter]] <- last.tree.list
      df_res.for.iter[[idx.iter]] <- df_res.list
      
    }
  }
  
  time_diff <- Sys.time() - st_time
  print(time_diff)
  return(list(trees = trees.for.iter, 
              df_res = df_res.for.iter))
}

BART_single_step <- function(X, res_theta, U1, U2, 
                             prop_mu, prop_sigma, prior_list,
                             moves.prob = NULL, starting.tree = NULL, 
                             diag = FALSE, cont.unif = TRUE, include.split,
                             obs.per.term = 1, empty.count.lim = 10,
                             log_like_fun, cop_type,
                             log_prior_fun,
                             fisher_fun_log,
                             fd_fun_log,
                             theta_param_1, theta_param_2){
  if(is.null(starting.tree)){
    rt_old <- generate_random_binary_tree_depth_free(1)
    rt_old <- assign_node_idx(rt_old)
    rt_old <- assign_split_rules(rt_old, X)
    rt_old <- assign_term_node_values_copula(rt_old, prop_mu, prop_sigma, cop_type) # check input
  } else {
    rt_old <- starting.tree
  }
  
  old_depth <- get_depth(rt_old)
  if(old_depth == 0){
    move.type <- 'grow'
  } else if(old_depth == 1){
    move.type <- sample(c('change', 'grow'), size = 1, prob = c(0.8, 0.2))
  } else{
    if(is.null(moves.prob)){
      move.type <- sample(c('swap', 'change', 'grow','prune'), size = 1, 
                          prob = c(0.4, 0.4, 0.1, 0.1))  
    } else{
      move.type <- sample(c('swap', 'change', 'grow','prune'), size = 1, 
                          prob = moves.prob) 
    }
  }
  
  move.flag = TRUE
  while(move.flag){
    #calculate obs at nodes with last tree
    if(diag == TRUE){
      nobs_diag <- vapply(get_terminal_nodes_idx(rt_old),
                          \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                    tree_top = rt_old, X.orig = X)), 0)
      #cat('Before move nobs at nodes:', nobs_diag,'\n')
    }
    
    new.tree.list <- tree_step_copula(move.type = move.type, # check input
                                      old_tree = rt_old, res_theta = res_theta, 
                                      X = X, U1 = U1, U2 = U2, 
                                      prior_list = prior_list, 
                                      cont.unif = cont.unif, 
                                      include.split = include.split,
                                      prop_mu = prop_mu, prop_sigma = prop_sigma, 
                                      log_like_fun = log_like_fun, 
                                      log_prior_fun = log_prior_fun,
                                      fisher_fun_log = fisher_fun_log,
                                      fd_fun_log = fd_fun_log,
                                      theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
    #cat('is tree char?' , is.character(new.tree.list), '\n')
    if(is.list(new.tree.list)){
      move.flag = FALSE} else{
        if(old_depth == 0){
          move.type <- 'grow'
        } else if(old_depth == 1){
          move.type <- sample(c('change', 'grow'), size = 1, prob = c(0.8, 0.2))
        } else{
          if(is.null(moves.prob)){
            move.type <- sample(c('swap', 'change', 'grow','prune'), size = 1, 
                                prob = c(0.4, 0.4, 0.1, 0.1))  
          } else{
            move.type <- sample(c('swap', 'change', 'grow','prune'), size = 1, 
                                prob = moves.prob) 
          }
        } 
      }
  }
  
  new.tree_mu.old <- new.tree.list$tree
  term_node_old <- get_terminal_nodes_idx(new.tree_mu.old)
  term_val_old <- sapply(term_node_old, function(i)unique(get_value_tree(new.tree_mu.old, get_obs_at_node(i, X, new.tree_mu.old,1,X))))


  new.tree_mu.new <- assign_term_node_values_cond_copula(tree_top = new.tree_mu.old,
                                                         X = X, U1 = U1, U2 = U2, res_theta = res_theta,
                                                         prop_mu = prop_mu,
                                                         prop_sigma = prop_sigma,
                                                         log_like_fun = log_like_fun,
                                                         log_prior_fun = log_prior_fun,
                                                         fisher_fun_log = fisher_fun_log,
                                                         fd_fun_log = fd_fun_log,
                                                         theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)

  term_node_new <- get_terminal_nodes_idx(new.tree_mu.new)
  term_val_new <- sapply(term_node_new, function(i)unique(get_value_tree(new.tree_mu.new, get_obs_at_node(i, X, new.tree_mu.new,1,X))))

  term_val_prior_old <- (sum(sapply(term_val_old, log_prior_fun)))
  term_val_prior_new <- (sum(sapply(term_val_new, log_prior_fun)))

  HR <- exp(cart_log_lik_copula(new.tree_mu.new, U1, U2, res_theta, X, log_like_fun = log_like_fun) - cart_log_lik_copula(new.tree_mu.old, U1, U2, res_theta, X, log_like_fun = log_like_fun) + term_val_prior_new - term_val_prior_old)

  if(runif(1)<HR){
    new.tree <- new.tree_mu.new
  }else{
    new.tree <- new.tree_mu.old
  }
  
  term_nodes <- get_terminal_nodes_idx(new.tree)
  
  term_vals <- sapply(term_nodes, function(i)unique(get_value_tree(new.tree, get_obs_at_node(i, X, new.tree,1,X))))
  
  sig_2 <- rinvgamma(1, shape = (0.01 + length(term_vals)/2), rate = (0.01 + sum(term_vals^2)/2))
  
  matrix.res <- matrix(c(move.type, get_depth(rt_old), get_num_terminal_nodes(rt_old),
                         new.tree.list$prior.ratio, new.tree.list$lik.ratio, 
                         new.tree.list$trans.ratio, new.tree.list$acc.prob, 
                         new.tree.list$accepted, sig_2),nrow = 1)
  
  colnames(matrix.res) <- c('move', 'old.depth', 'old.nterm', 'prior.ratio', 
                            'lik.ratio', 'trans.ratio', 'acc.prob', 'acceptance','sig2')
  
  return(list("tree" = new.tree, "df_res" = data.frame(matrix.res)))
  
}


tree_step_copula <- function(move.type, old_tree, X, U1, U2, res_theta, 
                             prior_list, cont.unif = TRUE, include.split,
                             obs.per.term = 1, empty.count.lim = 10,
                             prop_mu, prop_sigma,  
                             log_like_fun, 
                             log_prior_fun,
                             fisher_fun_log,
                             fd_fun_log,
                             theta_param_1, theta_param_2){
  empty.flag = TRUE
  empty.count = 0
  if(move.type == 'swap'){
    while(empty.flag & (empty.count <= empty.count.lim)){
      
      move.list <- swap_move_copula(old_tree, X, U1, U2, res_theta, 
                                    prop_mu = mean(get_value_tree(tree_top, obs.at.node_1)),
                                    prop_sigma = prop_sigma, 
                                    log_like_fun = log_like_fun, 
                                    log_prior_fun = log_prior_fun,
                                    fisher_fun_log = fisher_fun_log,
                                    fd_fun_log= fd_fun_log,
                                    theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
      
      cat('move=',move.list$move, ',idx = ',move.list$node.idx[1,1],',',move.list$node.idx[1,2],'\n')
      
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        cat('empty count: ', empty.count, '\n') }
      else{empty.flag = FALSE}
    }
    
  }
  if(move.type == 'change'){
    while(empty.flag& (empty.count <= empty.count.lim)){
      
      move.list <- change_move_copula(tree_top = old_tree, X = X, obs.per.term = obs.per.term, 
                                      cont.unif = cont.unif, U1, U2, 
                                      prop_mu = prop_mu, prop_sigma = prop_sigma, 
                                      log_like_fun = log_like_fun, 
                                      log_prior_fun = log_prior_fun,
                                      fisher_fun_log = fisher_fun_log,
                                      fd_fun_log = fd_fun_log,
                                      theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
      
      cat('move=',move.list$move, ',idx = ',move.list$node.idx,'\n')
      
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        cat('empty count: ', empty.count, '\n') }
      else{empty.flag = FALSE}  
    }
  }
  
  if(move.type == 'grow'){
    while(empty.flag& (empty.count <= empty.count.lim)){
      
      move.list <- grow_move_copula(tree_top = old_tree, X = X, U1 = U1, U2 = U2, res_theta = res_theta, # check input
                                    obs.per.term = obs.per.term,
                                    cont.unif = cont.unif,
                                    prop_mu = prop_mu, prop_sigma = prop_sigma, 
                                    log_like_fun = log_like_fun, 
                                    log_prior_fun = log_prior_fun,
                                    fisher_fun_log = fisher_fun_log,
                                    fd_fun_log = fd_fun_log,
                                    theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
      
      cat('move=',move.list$move, ',idx = ',move.list$node.idx,'\n')
      
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        cat('empty count: ', empty.count, '\n') }
      else{empty.flag = FALSE}  
    }
  }
  if(move.type == 'prune'& (empty.count <= empty.count.lim)){
    while(empty.flag){
      
      move.list <- prune_move_copula(tree_top = old_tree, 
                                     X = X, U1 = U1, U2 = U2, res_theta = res_theta,
                                     prop_mu = prop_mu, prop_sigma = prop_sigma, 
                                     log_like_fun = log_like_fun, 
                                     log_prior_fun = log_prior_fun,
                                     fisher_fun_log = fisher_fun_log,
                                     fd_fun_log = fd_fun_log,
                                     theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
      
      cat('move=',move.list$move, ',idx = ',move.list$node.idx,'\n')
      
      nobs.per.term <- vapply(get_terminal_nodes_idx(move.list$tree), 
                              \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                        tree_top = move.list$tree, X.orig = X)), 0)
      if(sum(nobs.per.term == 0) > 0){
        empty.count = empty.count + 1
        cat('empty count: ', empty.count, '\n') }
      else{empty.flag = FALSE}
    }
    
  }
  
  if(empty.count > empty.count.lim){
    print('empty count exceeded')
    return('empty count exceeded')
  }
  
  acc.prob <- acceptance.prob.list_copula(move_list = move.list, 
                                          old.tree = old_tree,
                                          X = X,
                                          U1 = U1,
                                          U2 = U2, res_theta,
                                          prior_input_list = prior_list,
                                          include.split = include.split,
                                          prop_mu = prop_mu, prop_sigma = prop_sigma,  
                                          log_like_fun = log_like_fun, 
                                          log_prior_fun = log_prior_fun,
                                          fisher_fun_log = fisher_fun_log,
                                          fd_fun_log = fd_fun_log,
                                          theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
  
  acceptance <- runif(1) <= acc.prob$alpha
  
  if(acceptance){
    
    return(list(tree = move.list$tree, 
                prior.ratio = acc.prob$prior.ratio,
                lik.ratio = acc.prob$lik.ratio,
                trans.ratio = acc.prob$trans.ratio,
                acc.prob = acc.prob$alpha,
                accepted = acceptance) ) 
    
  } else{
    return(list(tree = old_tree, 
                prior.ratio = acc.prob$prior.ratio,
                lik.ratio = acc.prob$lik.ratio,
                trans.ratio = acc.prob$trans.ratio,
                acc.prob = acc.prob$alpha,
                accepted = FALSE))
  }
}

assign_term_node_values_copula <- function(tree_top, prop_mu, prop_sigma, cop_type){ # check input
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  for(node.idx in term.node.idx){
    tree_top <- set_term_node_value_copula(node.idx, tree_top, prop_mu, prop_sigma, cop_type) # check input
  }
  return(tree_top)
}

assign_term_node_values_cond_copula <- function(tree_top, # check input
                                                X, U1, U2, res_theta, 
                                                prop_mu, prop_sigma, 
                                                log_like_fun, 
                                                log_prior_fun,
                                                fisher_fun_log,
                                                fd_fun_log,
                                                theta_param_1, theta_param_2){
  
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  for(node.idx in term.node.idx){
    obs.at.node <- get_obs_at_node(node.idx = node.idx, X = X, tree_top = tree_top, X.orig = X)
    U1.at.node <- U1[as.numeric(rownames(obs.at.node))]
    U2.at.node <- U2[as.numeric(rownames(obs.at.node))]
    res_theta.at.node <- res_theta[as.numeric(rownames(obs.at.node))]
    tree_top <- set_term_node_value_cond_copula(node.idx = node.idx, tree_top = tree_top, # check input
                                                X = X, U1 =U1, U2=U2, res_theta, U1.at.node = U1.at.node, 
                                                U2.at.node = U2.at.node, res_theta.at.node = res_theta.at.node,
                                                prop_mu = mean(get_value_tree(tree_top, obs.at.node)),
                                                prop_sigma = prop_sigma,  
                                                log_like_fun = log_like_fun, 
                                                log_prior_fun = log_prior_fun,
                                                fisher_fun_log = fisher_fun_log,
                                                fd_fun_log = fd_fun_log,
                                                theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
  }
  return(tree_top)
}

set_term_node_value_cond_copula <- function(node.idx, tree_top, # check input
                                            X, U1, U2, res_theta, res_theta.at.node = NULL,
                                            U1.at.node=NULL, U2.at.node=NULL, binary = FALSE, 
                                            prop_mu, prop_sigma, log_like_fun, 
                                            log_prior_fun,
                                            fisher_fun_log,
                                            fd_fun_log,
                                            theta_param_1, theta_param_2){
  
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      tree_top$value = sample.cond.mu.copula(tree_top = tree_top, # check input
                                             node.idx = node.idx,
                                             X = X, 
                                             U1 = U1, U2 = U2, res_theta, U1.at.node=U1.at.node, 
                                             U2.at.node=U2.at.node, res_theta.at.node = res_theta.at.node,
                                             prop_mu = prop_mu, prop_sigma = prop_sigma, 
                                             log_like_fun = log_like_fun, 
                                             log_prior_fun = log_prior_fun,
                                             fisher_fun_log = fisher_fun_log,
                                             fd_fun_log = fd_fun_log,
                                             theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)$val
      return(tree_top)
    } else {
      return(tree_top)  
    }
  }  else {
    tree.left <- set_term_node_value_cond_copula(node.idx, tree_top$left, 
                                                 X, U1, U2, res_theta, 
                                                 U1.at.node, U2.at.node, res_theta.at.node = res_theta.at.node,
                                                 prop_mu = prop_mu, prop_sigma = prop_sigma, 
                                                 log_like_fun = log_like_fun, 
                                                 log_prior_fun = log_prior_fun,
                                                 fisher_fun_log = fisher_fun_log,
                                                 fd_fun_log = fd_fun_log,
                                                 theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
    
    tree.right <- set_term_node_value_cond_copula(node.idx, tree_top$right, 
                                                  X, U1, U2, res_theta, 
                                                  U1.at.node, U2.at.node, res_theta.at.node = res_theta.at.node,
                                                  prop_mu = prop_mu, prop_sigma = prop_sigma, 
                                                  log_like_fun = log_like_fun, 
                                                  log_prior_fun = log_prior_fun,
                                                  fisher_fun_log = fisher_fun_log,
                                                  fd_fun_log = fd_fun_log,
                                                  theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
    
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}

# sample conditional mu assuming binary observations and beta prior
sample.cond.mu.copula <- function(tree_top = NULL, 
                                  node.idx = NULL, 
                                  X = NULL, 
                                  U1 = NULL, 
                                  U2 = NULL, res_theta, 
                                  U1.at.node = NULL,
                                  U2.at.node = NULL,
                                  res_theta.at.node = NULL,
                                  prop_mu, prop_sigma, 
                                  log_like_fun, 
                                  log_prior_fun,
                                  fisher_fun_log,
                                  fd_fun_log,
                                  theta_param_1, theta_param_2){
  
  if(is.null(U1.at.node)){
    obs.at.node <- get_obs_at_node(node.idx = node.idx, X = X, tree_top = tree_top, X.orig = X)
    U1.at.node <- U1[as.numeric(rownames(obs.at.node))]
    U2.at.node <- U2[as.numeric(rownames(obs.at.node))]
    res_theta.at.node <- res_theta[as.numeric(rownames(obs.at.node))]
    prop_mu = mean(get_value_tree(tree_top, obs.at.node))
  } 
  
  ##############################################################################
  
  log_post <- function(x)logposterior(x, U1.at.node, U2.at.node, 
                                      log_like_fun = function(rho, u, v)log_like_fun((rho+res_theta.at.node),u,v), 
                                      log_prior_fun = log_prior_fun)
  
  prop_params <- laplace_approx(prop_mu, theta_param_2, U1.at.node, U2.at.node, function(rho,u,v)fisher_fun_log((rho+res_theta.at.node),u,v), function(rho,u,v)fd_fun_log((rho+res_theta.at.node),u,v), log_post)
  
  prop_mu_new <- prop_params$prop_mu
  
  prop_sd <- prop_params$prop_sd
  
  proposal = rnorm(1, prop_mu_new, prop_sd)
  
  return(list("val"= proposal, "mu" = prop_mu_new, "sigma" = prop_sd))
}

grow_move_copula <- function(tree_top, X, U1, U2, res_theta, cont.unif = TRUE, obs.per.term = 1, 
                             prop_mu, prop_sigma, 
                             log_like_fun, 
                             log_prior_fun,
                             fisher_fun_log,
                             fd_fun_log,
                             theta_param_1, theta_param_2){ # check input
  
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  
  nobs.at.nodes <- vapply(term.node.idx, \(x) nrow(get_obs_at_node(node.idx = x, X = X, tree_top = tree_top, X.orig = X)),0)
  term.node.idx <- term.node.idx[nobs.at.nodes > 1]
  check.valid.idx <- vapply(term.node.idx, \(x) check_grow_valid(x, tree_top, X), TRUE)
  term.node.idx <- term.node.idx[check.valid.idx]
  
  if(length(term.node.idx) == 0){
    return(list(tree = tree_top, node.idx = NULL, move = 'grow', 
                valid.pred = NULL, 
                valid.split = NULL))
  }  else if(length(term.node.idx) > 1){ 
    grow.idx <- sample(term.node.idx, 1, replace = FALSE)
  } else if(length(term.node.idx) == 1){
    grow.idx <- term.node.idx
  }
  
  new.cond.list <- gen_node_condition(grow.idx, tree_top, X, obs.per.term, cont.unif = cont.unif, for.grow = TRUE)
  obs.at.node <- new.cond.list$obs.at.node
  U1.at.parent <- U1[as.numeric(rownames(obs.at.node))]
  U2.at.parent <- U2[as.numeric(rownames(obs.at.node))]
  res_theta.at.parent <- res_theta[as.numeric(rownames(obs.at.node))]
  if(is.numeric(new.cond.list$cond$x.val)){
    idx.left <- obs.at.node[,new.cond.list$cond$x.idx] <= new.cond.list$cond$x.val
  } else {
    idx.left <- obs.at.node[, new.cond.list$cond$x.idx] %in% new.cond.list$cond$x.val
  }
  idx.right <- !idx.left
  
  U1.at.left <- U1.at.parent[idx.left]
  U1.at.right <- U1.at.parent[idx.right]
  U2.at.left <- U2.at.parent[idx.left]
  U2.at.right <- U2.at.parent[idx.right]
  res_theta.at.right <- res_theta.at.parent[idx.right]
  res_theta.at.left <- res_theta.at.parent[idx.left]
  
  term.node.value.left <- sample.cond.mu.copula(U1.at.node = U1.at.left, U2.at.node = U2.at.left, res_theta = res_theta,
                                                res_theta.at.node = res_theta.at.left,
                                                prop_mu = mean(get_value_tree(tree_top, obs.at.node[idx.left, new.cond.list$cond$x.idx, drop = FALSE])),
                                                prop_sigma = prop_sigma, 
                                                log_like_fun = log_like_fun, 
                                                log_prior_fun = log_prior_fun,
                                                fisher_fun_log = fisher_fun_log,
                                                fd_fun_log = fd_fun_log,
                                                theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
  
  term.node.value.right <- sample.cond.mu.copula(U1.at.node = U1.at.right, U2.at.node = U2.at.right, res_theta = res_theta,
                                                 res_theta.at.node = res_theta.at.right,
                                                 prop_mu = mean(get_value_tree(tree_top, obs.at.node[idx.right, new.cond.list$cond$x.idx, drop = FALSE])),
                                                 prop_sigma = prop_sigma, 
                                                 log_like_fun = log_like_fun, 
                                                 log_prior_fun = log_prior_fun,
                                                 fisher_fun_log = fisher_fun_log,
                                                 fd_fun_log = fd_fun_log,
                                                 theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
  
  tree_top_grow <- grow_terminal(node.idx = grow.idx, 
                                 tree_top = tree_top, 
                                 new.cond = new.cond.list$cond, 
                                 value.l = term.node.value.left$val,
                                 value.r = term.node.value.right$val)
  tree_top_grow <- assign_node_idx(tree_top_grow)
  
  return(list(tree = tree_top_grow, node.idx = grow.idx, move = 'grow', 
              valid.pred = new.cond.list$valid.pred, 
              valid.split = new.cond.list$valid.split,
              prop_right = term.node.value.right,
              prop_left = term.node.value.left))
}

prune_move_copula <- function(tree_top, X, U1, U2, res_theta, 
                              prop_mu, prop_sigma, 
                              log_like_fun, 
                              log_prior_fun,
                              fisher_fun_log,
                              fd_fun_log,
                              theta_param_1, theta_param_2){ 
  
  prune.node.idx <- get_prune_idx(tree_top)
  
  if(length(prune.node.idx) == 1){
    prune.idx <- prune.node.idx
  } else{
    prune.idx <- sample(prune.node.idx, 1, replace = FALSE)
  }
  
  prune.value <- sample.cond.mu.copula(tree_top = tree_top, # check input
                                       node.idx = prune.idx,
                                       res_theta = res_theta,
                                       X = X, 
                                       U1 = U1,
                                       U2 = U2,
                                       prop_mu = prop_mu, prop_sigma = prop_sigma, 
                                       log_like_fun = log_like_fun, 
                                       log_prior_fun = log_prior_fun,
                                       fisher_fun_log = fisher_fun_log,
                                       fd_fun_log = fd_fun_log,
                                       theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
  
  tree_top_prune <- prune_terminal(node.idx = prune.idx, 
                                   tree_top = tree_top, 
                                   value.node = prune.value$val)
  tree_top_prune <- assign_node_idx(tree_top_prune)
  return(list(tree = tree_top_prune, move = 'prune', node.idx = prune.idx,
              prune.node.idx = prune.node.idx,
              list_prop = prune.value))
}

# perform swap move
swap_move_copula <- function(tree_top, X, U1, U2, res_theta, 
                             prop_mu, prop_sigma, 
                             log_like_fun, 
                             log_prior_fun,
                             fisher_fun_log,
                             fd_fun_log,
                             theta_param_1, theta_param_2){
  
  swap.node.idx <- get_valid_swap_idx(tree_top)
  
  if(nrow(swap.node.idx) == 0){
    return(list(tree = tree_top, move = 'swap', valid.pred = NULL, valid.split = NULL))
  } else {
    swap.idx.row <- sample(1:nrow(swap.node.idx), 1)
    swap.idx <- swap.node.idx[swap.idx.row, ]
    swap.tree_top <- swap_node_condition(swap.idx[1,1], swap.idx[1,2], tree_top)
    
    return(list(tree = swap.tree_top, move = 'swap', node.idx = swap.idx, valid.pred = NULL, valid.split = NULL))
  }
}

# perform change move
change_move_copula <- function(tree_top, X, obs.per.term = 2, cont.unif = TRUE, U1,U2,
                               prop_mu, prop_sigma, 
                               log_like_fun, 
                               log_prior_fun,
                               fisher_fun_log,
                               fd_fun_log,
                               theta_param_1, theta_param_2){
  
  internal.node.idx <- get_internal_nodes_idx(tree_top)
  change.idx <- sample(internal.node.idx, 1)
  new.cond <- gen_node_condition(change.idx, tree_top, X, obs.per.term, cont.unif)
  change.tree_top <- set_node_condition(change.idx, tree_top, new.cond$cond)
  
  return(list(tree = change.tree_top, move = 'change', 
              node.idx = change.idx, 
              valid.pred = new.cond$valid.pred, 
              valid.split = new.cond$valid.split,
              n.left = new.cond$n.left,
              n.right = new.cond$n.right))
}


acceptance.prob.list_copula <- function(move_list, old.tree, X, U1, U2, res_theta, prior_input_list,
                                        include.split, cont.unif = TRUE, 
                                        prop_mu, prop_sigma, 
                                        log_like_fun, 
                                        log_prior_fun,
                                        fisher_fun_log,
                                        fd_fun_log,
                                        theta_param_1, theta_param_2){ 
  
  if(identical(move_list$tree, old.tree)){
    return(list(prior.ratio = 1, lik.ratio = 1, 
                trans.ratio = 1, alpha = 1))
  }
  
  prior.tree.old <- prior_input_list$fun(old.tree, prior_input_list$param[1], prior_input_list$param[2])
  prior.tree.new <- prior_input_list$fun(move_list$tree, prior_input_list$param[1], prior_input_list$param[2])
  
  prior.old <- prior.tree.old*prior.split.rule(old.tree, X, cont.unif = cont.unif)
  prior.new <- prior.tree.new*prior.split.rule(move_list$tree, X, cont.unif = cont.unif)
  
  # likelihood 
  log.lik.old <- cart_log_lik_copula(tree_top = old.tree, U1, U2, res_theta, X, log_like_fun = log_like_fun)
  log.lik.new <- cart_log_lik_copula(tree_top = move_list$tree, U1, U2, res_theta, X, log_like_fun = log_like_fun)
  
  
  #transition probabilities
  if(move_list$move == 'swap'){
    prob.old.to.new <- 1#transition.p.swap(old.tree)
    prob.new.to.old <- 1#transition.p.swap(move_list$tree)
    
    term_val_prior_new <- 1
    term_val_prior_old <- 1
    prop_dist_old_val <- 1
    prop_dist_new_val <- 1
  } else if(move_list$move == 'change'){
    prob.old.to.new <- 1#transition.p.change(move_list)
    prob.new.to.old <- 1#transition.p.change.reverse(move_list, X, old.tree)
    
    term_val_prior_new <- 1
    term_val_prior_old <- 1
    prop_dist_old_val <- 1
    prop_dist_new_val <- 1
  } else if(move_list$move == 'grow'){
    prob.old.to.new <- transition.p.grow(move_list, include.split = include.split)
    prob.new.to.old <- transition.p.grow.reverse(move_list)
    
    term_val_old <- sample.cond.mu.copula(tree_top = old.tree, # check input
                                                         node.idx = move_list$node.idx,
                                                         res_theta = res_theta,
                                                         X = X, 
                                                         U1 = U1,
                                                         U2 = U2,
                                                         prop_mu = prop_mu, prop_sigma = prop_sigma, 
                                                         log_like_fun = log_like_fun, 
                                                         log_prior_fun = log_prior_fun,
                                                         fisher_fun_log = fisher_fun_log,
                                                         fd_fun_log = fd_fun_log,
                                                         theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
    
    term_val_new_right <- move_list$prop_right
    term_val_new_left <- move_list$prop_left
    
    term_val_prior_old <- (sum(sapply(term_val_old$val, log_prior_fun)))
    term_val_prior_new <- (sum(sapply(c(term_val_new_left$val,term_val_new_right$val), log_prior_fun)))
    
    prop_dist_old_val <- sum(dnorm(term_val_old$val, mean = term_val_old$mu, sd = term_val_old$sigma, log = T))
    prop_dist_new_val <- sum(dnorm(c(term_val_new_left$val,term_val_new_right$val), mean = c(term_val_new_left$mu,term_val_new_right$mu), sd = c(term_val_new_left$sigma,term_val_new_right$sigma), log = T))
    
  } else if(move_list$move == 'prune'){
    prob.old.to.new <- transition.p.prune(move_list)
    prob.new.to.old <- transition.p.prune.reverse(move_list, old.tree, X, obs.per.term = 1,
                                                  include.split = include.split)
    
    ######
    new.cond.list <- get_node_condition(move_list$node.idx, old.tree)
    obs.at.node <- get_obs_at_node(move_list$node.idx, X, old.tree,1,X)
    U1.at.parent <- U1[as.numeric(rownames(obs.at.node))]
    U2.at.parent <- U2[as.numeric(rownames(obs.at.node))]
    res_theta.at.parent <- res_theta[as.numeric(rownames(obs.at.node))]
    if(is.numeric(new.cond.list$x.val)){
      idx.left <- obs.at.node[,new.cond.list$x.idx] <= new.cond.list$x.val
    } else {
      idx.left <- obs.at.node[, new.cond.list$x.idx] %in% new.cond.list$x.val
    }
    idx.right <- !idx.left
    
    U1.at.left <- U1.at.parent[idx.left]
    U1.at.right <- U1.at.parent[idx.right]
    U2.at.left <- U2.at.parent[idx.left]
    U2.at.right <- U2.at.parent[idx.right]
    res_theta.at.right <- res_theta.at.parent[idx.right]
    res_theta.at.left <- res_theta.at.parent[idx.left]
    
    term_val_old_left <- sample.cond.mu.copula(U1.at.node = U1.at.left, U2.at.node = U2.at.left, res_theta = res_theta,
                                                  res_theta.at.node = res_theta.at.left,
                                                  prop_mu = mean(get_value_tree(old.tree, obs.at.node[idx.left, new.cond.list$x.idx, drop = FALSE])),
                                                  prop_sigma = prop_sigma, 
                                                  log_like_fun = log_like_fun, 
                                                  log_prior_fun = log_prior_fun,
                                                  fisher_fun_log = fisher_fun_log,
                                                  fd_fun_log = fd_fun_log,
                                                  theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
    
    term_val_old_right <- sample.cond.mu.copula(U1.at.node = U1.at.right, U2.at.node = U2.at.right, res_theta = res_theta,
                                                   res_theta.at.node = res_theta.at.right,
                                                   prop_mu = mean(get_value_tree(old.tree, obs.at.node[idx.right, new.cond.list$x.idx, drop = FALSE])),
                                                   prop_sigma = prop_sigma, 
                                                   log_like_fun = log_like_fun, 
                                                   log_prior_fun = log_prior_fun,
                                                   fisher_fun_log = fisher_fun_log,
                                                   fd_fun_log = fd_fun_log,
                                                   theta_param_1 = theta_param_1, theta_param_2 = theta_param_2)
    
    term_val_new <- move_list$list_prop
    
    term_val_prior_old <- (sum(sapply(c(term_val_old_right$val,term_val_old_left$val), log_prior_fun)))
    term_val_prior_new <- (sum(sapply(term_val_new$val, log_prior_fun)))
    
    prop_dist_old_val <- sum(dnorm(c(term_val_old_right$val,term_val_old_left$val), mean = c(term_val_old_right$mu,term_val_old_left$mu), sd = c(term_val_old_right$sigma,term_val_old_left$sigma), log = T))
    prop_dist_new_val <- sum(dnorm(term_val_new$val, mean = term_val_new$mu, sd = term_val_new$sigma, log = T))
  } else {
    stop('Unknown move')
  }
  prior.ratio <- (prior.new/prior.old) * exp(term_val_prior_new - term_val_prior_old)
  lik.ratio <- exp(log.lik.new - log.lik.old)
  trans.ratio <- (prob.new.to.old/prob.old.to.new) * exp(prop_dist_old_val-prop_dist_new_val)
  acc.prob <- (prior.new/prior.old) * (prob.new.to.old/prob.old.to.new) *
    exp(log.lik.new - log.lik.old + term_val_prior_new - term_val_prior_old + prop_dist_old_val-prop_dist_new_val)
  
  return(list(prior.ratio = prior.ratio, lik.ratio = lik.ratio, 
              trans.ratio = trans.ratio, alpha = min(1, acc.prob)))
}

cart_log_lik_copula <- function(tree_top, U1, U2, res_theta, X, log_like_fun = log_like_fun){ # check input
  tree.at.obs <- get_value_tree(tree_top, X)
  
  log.prob.obs <- log_like_fun(tree.at.obs + res_theta, U1, U2)
  return(log.prob.obs)
}

set_term_node_value_copula <- function(node.idx, tree_top, 
                                       prop_mu, prop_sigma, cop_type){
  
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      if(cop_type == "gauss" || cop_type == "t"){
        tree_top$value = link_gauss(rnorm(1, mean = prop_mu, sd = prop_sigma))
      } else if(cop_type == "gumbel") {
        tree_top$value = link_gumbel(rnorm(1, mean = prop_mu, sd = prop_sigma))
      } else if(cop_type == "clayton") {
        tree_top$value = link_clayton(rnorm(1, mean = prop_mu, sd = prop_sigma)) 
      }
      return(tree_top)
    } else {
      return(tree_top)  
    }
  }  else {
    
    tree.left <- set_term_node_value_copula(node.idx, tree_top$left, prop_mu, prop_sigma, cop_type)
    
    tree.right <- set_term_node_value_copula(node.idx, tree_top$right, prop_mu, prop_sigma, cop_type)
    
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}

loglik_gauss <- function(rho, u, v) {
  
  if(any(abs(rho)>=1)) return(-Inf)
  
  log_lik <- sum(log(BiCopPDF(u,v,1,rho)))
  return(log_lik)
}

FI_gauss <- function(theta, u, v) {
  
  if(any(abs(link_gauss(theta))>=1)) return(0)
  
  f_deriv <- BiCopDeriv(u,v,1,link_gauss(theta))
  f_deriv2 <- BiCopDeriv2(u,v,1,link_gauss(theta))
  
  f<- BiCopPDF(u,v,1,link_gauss(theta))
  
  FI <- -(f_deriv2 / f - f_deriv^2/f^2) * (2*sigmoid(theta) * (1-sigmoid(theta)))^2 - (f_deriv/f) * (2*sigmoid(theta) * (1-sigmoid(theta)) * (1-2*sigmoid(theta)))
  
  return(sum(FI))
}

FD_log_gauss <- function(theta, u, v) {
  
  if(any(abs(link_gauss(theta))>=1)) return(0)
  
  f_deriv <- BiCopDeriv(u,v,1,link_gauss(theta))
  
  f<- BiCopPDF(u,v,1,link_gauss(theta))
  
  FD <- (f_deriv / f) * (2*sigmoid(theta) * (1-sigmoid(theta)))
  
  return(sum(FD))
}

loglik_t <- function(rho, u, v, df = 3) {
  
  if(any(abs(rho)>=1)) return(-Inf)
  
  log_lik <- sum(log(BiCopPDF(u,v,2,rho, par2 = df)))
  
  return(log_lik)
}

FI_t <- function(theta, u, v, df = 3) {
  
  if(any(abs(link_t(theta))>=1)) return(0)
  
  f_deriv <- BiCopDeriv(u,v,2,link_t(theta), par2 = df)
  f_deriv2 <- BiCopDeriv2(u,v,2,link_t(theta), par2 = df)
  
  f<- BiCopPDF(u,v,2,link_t(theta), par2 = df)
  
  FI <- -(f_deriv2 / f - f_deriv^2/f^2) * (2*sigmoid(theta) * (1-sigmoid(theta)))^2 - (f_deriv/f) * (2*sigmoid(theta) * (1-sigmoid(theta)) * (1-2*sigmoid(theta)))
  
  return(sum(FI))
}

FD_log_t <- function(theta, u, v, df = 3) {
  
  if(any(abs(link_t(theta))>=1)) return(0)
  
  f_deriv <- BiCopDeriv(u,v,2,link_t(theta), par2 = df)
  
  f<- BiCopPDF(u,v,2,link_t(theta), par2 = df)
  
  FD <- (f_deriv / f) * (2*sigmoid(theta) * (1-sigmoid(theta)))
  
  return(sum(FD))
}

loglik_clayton <- function(theta, u, v) {
  if(any(theta<=0)||any(theta>=28)) return(-Inf)
  
  log_density <- log(BiCopPDF(u,v,3,theta))
  
  return(sum(log_density))
}

FI_clayton <- function(theta, u, v) {
  if(any(link_clayton(theta)<=1e-10)||any(link_clayton(theta)>=28)) return(0)
  
  f_deriv <- BiCopDeriv(u,v,3,link_clayton(theta))
  f_deriv2 <- BiCopDeriv2(u,v,3,link_clayton(theta))
  
  f<- BiCopPDF(u,v,3,link_clayton(theta))
  
  FI <- -(f_deriv2 / f - f_deriv^2/f^2) * (28*sigmoid(theta) * (1-sigmoid(theta)))^2 - (f_deriv/f) * (28*sigmoid(theta) * (1-sigmoid(theta)) * (1-2*sigmoid(theta)))
  
  return(sum(FI))
}

FD_log_clayton <- function(theta, u, v) {
  
  if(any(link_clayton(theta)<=1e-10)||any(link_clayton(theta)>=28)) return(0)
  
  f_deriv <- BiCopDeriv(u,v,3,link_clayton(theta))
  
  f<- BiCopPDF(u,v,3,link_clayton(theta))
  
  FD <- (f_deriv / f) * (28*sigmoid(theta) * (1-sigmoid(theta)))
  
  return(sum(FD))
}

loglik_gumbel <- function(theta, u, v) {
  
  if(any(theta<=1)||any(theta>=17)) return(-Inf)
  
  log_density <- log(BiCopPDF(u,v,4,theta))
  
  return(sum(log_density))
}

FI_gumbel <- function(theta, u, v) {
  
  if(any(link_gumbel(theta)<=1)||any(link_gumbel(theta) >=17)) return(0)
  
  f_deriv <- BiCopDeriv(u,v,4,(link_gumbel(theta)))
  f_deriv2 <- BiCopDeriv2(u,v,4,(link_gumbel(theta)))
  
  f<- BiCopPDF(u,v,4,(link_gumbel(theta)))
  
  FI <- -(f_deriv2 / f - f_deriv^2/f^2) * (16*sigmoid(theta) * (1-sigmoid(theta)))^2 - (f_deriv/f) * (16*sigmoid(theta) * (1-sigmoid(theta)) * (1-2*sigmoid(theta)))
  
  return(sum(FI))
}

FD_log_gumbel <- function(theta, u, v) {
  
  if(any(link_gumbel(theta)<=1)||any(link_gumbel(theta) >=17)) return(0)
  
  f_deriv <- BiCopDeriv(u,v,4,(link_gumbel(theta)))
  
  f<- BiCopPDF(u,v,4,(link_gumbel(theta)))
  
  FD <- (f_deriv / f) * (16*sigmoid(theta) * (1-sigmoid(theta)))
  
  return(sum(FD))
}

log_prior_tbeta <- function(rho, bound = 1, theta_param_1, theta_param_2) {
  if (any(abs(rho) >= 1)) {
    return(rep(-Inf, length(rho)))
  }
  
  return((theta_param_1 - 1)*log(bound + rho) + (theta_param_2 - 1)*log(bound - rho))
}

logposterior <- function(rho, u, v, log_like_fun, log_prior_fun){
  return(log_like_fun(rho,u,v) + log_prior_fun(rho))
}

param_gauss <- function(tau) return(sin(tau*pi/2))
param_t <- function(tau) return(sin(tau*pi/2))
param_gumbel <- function(tau) return(1/(1-tau))
param_clayton <- function(tau) return((2*tau)/(1-tau))

conv_diag <- function(post_data, panel_name, n.burn, n.thin){
  
  post_mat <- (post_data %>% filter(panel.name == panel_name))$y
  
  post_mat_red <- na.omit((post_mat[n.burn+1:length(post_mat)])[c(rep(NA,(n.thin-1)), TRUE)])
  
  post_data_mcmc <- mcmc(post_mat_red, start = 1)
  
  auto_corr_sum <- autocorr.diag(post_data_mcmc)
  
  eff_size <- effectiveSize(post_data_mcmc)
  
  geweke_score <- geweke.diag(post_data_mcmc)
  
  
  summ_dat <- cbind(auto_corr_sum[2], 100*(eff_size/length(post_data_mcmc)), geweke_score$z)
  
  colnames(summ_dat) <- c("auto-correlation", "ess (%)", "geweke-score")
  rownames(summ_dat) <- panel_name
  
  return(summ_dat)
  
}

laplace_approx <- function(prop_mu, theta_param_2, u, v, fisher_fun_log, fd_fun_log, log_post){
  prior_fisher <- function(mu) (fisher_fun_log(mu, u, v) + 1/(theta_param_2)^2)
  prior_fd <- function(mu) (fd_fun_log(mu, u, v) - mu/(theta_param_2)^2)
  
  prop_mu <- optim(prop_mu, fn = log_post, gr = prior_fd, method = "BFGS", control = list(fnscale = -1))$par
  
  # n.iter <- 1
  # while((n.iter<21)&(abs(prior_fd(prop_mu)) > (sqrt(abs(prior_fisher(prop_mu)))/10))){
  #   prop_mu <- prop_mu + prior_fd(prop_mu)/(prior_fisher(prop_mu))
  #   n.iter <- n.iter + 1
  # }
  
  prop_sd <- 1/sqrt(abs(prior_fisher(prop_mu))) #+ 1e-4
  
  return(list("prop_mu" = prop_mu, "prop_sd" = prop_sd))
}

nterm_BART <- function(bart_tree_list){
  ll <- lapply(1:length(bart_tree_list$trees), \(idx) data.frame(idx = idx,
                                                                 nn = vapply(bart_tree_list$trees[[idx]], 
                                                                             \(x) get_num_terminal_nodes(x),
                                                                             0),
                                                                 trees = factor( 1:length(bart_tree_list$trees[[idx]]))))
  Reduce(rbind, ll)
}

depth_BART <- function(bart_tree_list){
  ll <- lapply(1:length(bart_tree_list$trees), \(idx) data.frame(idx = idx,
                                                                 nn = vapply(bart_tree_list$trees[[idx]], 
                                                                             \(x) get_depth(x),
                                                                             0),
                                                                 trees = factor( 1:length(bart_tree_list$trees[[idx]]))))
  Reduce(rbind, ll)
}

acc_BART <- function(bart_tree_list){
  ll <- lapply(1:length(bart_tree_list$df.res), \(idx) data.frame(idx = idx,
                                                                 nn = unlist(lapply(bart_tree_list$df.res[[idx]], 
                                                                                    \(x)as.factor(x$acceptance))),
                                                                 trees = factor( 1:length(bart_tree_list$df.res[[idx]]))))
  Reduce(rbind, ll)
}

BART_calculate_pred <- function(tree_list, X){
  gx <- lapply(tree_list, \(x) get_value_tree(x, X))
  Reduce('+', gx)
}

link_gauss <- function(rho) {
  return((2*sigmoid(rho)-1))
}

link_t <- function(rho) {
  return((2*sigmoid(rho)-1))
}

link_clayton <- function(x) {
  return(28 / (1 + exp(-x)))
}

link_gumbel <- function(x) {
  return(16 / (1 + exp(-x)) + 1)
}