library(parallel)
library(mc2d)
library(pracma)
library(MCMCpack)

multichain_MCMC_copula <- function(n.iter = 6000, 
                                   n.burn = 1000, 
                                   n.chain = 1,
                                   n.tree = 10,
                                   X, 
                                   U1, U2,
                                   prior_list,
                                   moves.prob = NULL, 
                                   starting.tree = NULL,
                                   n.cores = min(n.chain, 5),
                                   cont.unif = TRUE,
                                   include.split,
                                   prop_mu, prop_sigma,
                                   theta_param_1, theta_param_2,
                                   var_param_1, var_param_2,
                                   prior_type, cop_type, adapt = T){
  
  st_time = Sys.time()
  
  chain.list <- mclapply(1:n.chain,
                         \(x) MCMC_copula(n.iter = n.iter, 
                                          n.burn = n.burn,
                                          n.tree = n.tree,
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
                                          var_param_1 = var_param_1, var_param_2 = var_param_2,
                                          prior_type = prior_type, cop_type = cop_type,
                                          adapt = adapt),
                         mc.cores = n.cores)
  tree_list <- lapply(chain.list, function(x) x$trees)
  tree_list_comb <- Reduce(c, tree_list)
  
  df_list <- lapply(chain.list,
                    function(x) x$df_res)
  df.res <- Reduce(c, df_list)
  
  time_diff <- Sys.time() - st_time
  print(time_diff)
  return(list(trees = tree_list_comb, df_res = df.res))
}

MCMC_copula <- function(n.iter, 
                        n.burn, 
                        n.tree = 10, 
                        X, 
                        U1, U2, 
                        prior_list,
                        moves.prob = NULL, 
                        starting.tree = NULL, 
                        diag = FALSE, 
                        cont.unif = TRUE, 
                        include.split,
                        prop_mu, prop_sigma, 
                        theta_param_1, theta_param_2, 
                        var_param_1, var_param_2, 
                        prior_type, cop_type, adapt = T){
  
  ##############################################################################
  
  log_like_fun <- switch(cop_type,
                         gauss = function(rho, u, v) loglik_gauss(link_gauss(rho), u, v),
                         t = function(rho, u, v) loglik_t(link_t(rho), u, v, df = 3),
                         clayton = function(rho, u, v) loglik_clayton(link_clayton(rho), u, v),
                         gumbel = function(rho, u, v) loglik_gumbel(link_gumbel(rho), u, v),
                         frank = function(rho, u, v) loglik_frank(link_frank(rho), u, v),
                         stop("Unsupported copula type")  # default case
  )
  
  log_prior_fun <- switch(prior_type,
                          B = function(rho, theta_param_1, theta_param_2)
                            log_prior_tbeta(rho, bound = 1, theta_param_1, theta_param_2),
                          
                          LN = function(rho, theta_param_1, theta_param_2)
                            dlnorm(rho, meanlog = theta_param_1, sdlog = theta_param_2, log = TRUE),
                          
                          G = function(rho, theta_param_1, theta_param_2)
                            dgamma(rho, shape = theta_param_1, rate = theta_param_2, log = TRUE),
                          
                          N = function(rho, theta_param_1, theta_param_2)
                            dnorm(rho, mean = theta_param_1, sd = theta_param_2, log = TRUE),
                          
                          stop("Unsupported prior type")
  )
  ##
  
  # initialise output
  starting.tree.list <- vector("list", n.tree)
  last.tree.list <- starting.tree.list
  df_res.list <- vector("list", n.tree)
  idx.tree.vec <- seq_len(n.tree)
  res_theta_list <- matrix(0, nrow = nrow(X), ncol = n.tree)
  
  sig_2 <- rinvgamma(n.tree, shape = var_param_1, scale = var_param_2)
  
  # params for computation
  
  trees.for.iter <- vector("list", n.iter)
  df_res.for.iter <- vector("list", n.iter)
  theta_pred_list <- array(NA_real_, dim = c(nrow(X), n.tree, n.iter))
  
  st_time = Sys.time()
  
  # assign("sample.cond.mu.copula", sample.cond.mu.copula_before_burn)
  
  for(idx.iter in 1:n.burn){
    cat('Iteration: ', idx.iter, '\n')
    
    for(idx.tree in idx.tree.vec){
      
      res_theta <- rowSums(res_theta_list[,idx.tree.vec != idx.tree, drop = FALSE])
      
      # cat('Tree: ', idx.tree, '\n')
      old_tree <- last.tree.list[[idx.tree]]
      
      new_tree_single <- BART_single_step(X = X, 
                                          res_theta = res_theta, 
                                          U1 = U1, U2 = U2, 
                                          prop_mu = prop_mu, prop_sigma = prop_sigma, 
                                          prior_list = prior_list,
                                          moves.prob = moves.prob, 
                                          starting.tree = old_tree, 
                                          diag = diag, 
                                          cont.unif = cont.unif, 
                                          include.split = include.split,
                                          obs.per.term = 1, 
                                          empty.count.lim = 10, 
                                          log_like_fun = log_like_fun, 
                                          log_prior_fun = log_prior_fun,
                                          theta_param_1 = theta_param_1, theta_param_2 =  sqrt(sig_2[idx.tree]),
                                          var_param_1 = var_param_1, var_param_2 = var_param_2,
                                          theta_var_list = NULL,
                                          idx.tree = idx.tree)
      
      new_tree <- new_tree_single$tree
      
      sig_2[idx.tree] <- as.numeric(new_tree_single$df_res$sig2)
      
      last.tree.list[[idx.tree]] <- new_tree
      df_res.list[[idx.tree]] <- new_tree_single$df_res
      res_theta_list[,idx.tree] <- get_value_tree(new_tree, X)
      
    }
    trees.for.iter[[idx.iter]] <- last.tree.list
    df_res.for.iter[[idx.iter]] <- df_res.list
    theta_pred_list[, , idx.iter] <- res_theta_list
    
  }
  
  var_list = lapply(idx.tree.vec, function(x)var_adapt(theta_pred_list[,x,1:n.burn]))
  
  for(idx.iter in (n.burn+1):n.iter){
    cat('Iteration: ', idx.iter, '\n')
    
    for(idx.tree in idx.tree.vec){
      
      if(adapt){
        theta_var_list_in <- var_list[[idx.tree]]$cov
      }else {
        theta_var_list_in <- NULL
      }
      
      res_theta <- rowSums(res_theta_list[,idx.tree.vec != idx.tree, drop = FALSE])
      
      # cat('Tree: ', idx.tree, '\n')
      old_tree <- last.tree.list[[idx.tree]]
      
      new_tree_single <- BART_single_step(X = X, 
                                          res_theta = res_theta, 
                                          U1 = U1, U2 = U2, 
                                          prop_mu = prop_mu, prop_sigma = prop_sigma, 
                                          prior_list = prior_list,
                                          moves.prob = moves.prob, 
                                          starting.tree = old_tree, 
                                          diag = diag, 
                                          cont.unif = cont.unif, 
                                          include.split = include.split,
                                          obs.per.term = 1, 
                                          empty.count.lim = 10, 
                                          log_like_fun = log_like_fun, 
                                          log_prior_fun = log_prior_fun,
                                          theta_param_1 = theta_param_1, theta_param_2 =  sqrt(sig_2[idx.tree]),
                                          var_param_1 = var_param_1, var_param_2 = var_param_2,
                                          theta_var_list = theta_var_list_in,
                                          idx.tree = idx.tree)
      
      new_tree <- new_tree_single$tree
      
      sig_2[idx.tree] <- as.numeric(new_tree_single$df_res$sig2)
      
      last.tree.list[[idx.tree]] <- new_tree
      df_res.list[[idx.tree]] <- new_tree_single$df_res
      res_theta_list[,idx.tree] <- get_value_tree(new_tree, X)
      
    }
    trees.for.iter[[idx.iter]] <- last.tree.list
    df_res.for.iter[[idx.iter]] <- df_res.list
    theta_pred_list[, , idx.iter] <- res_theta_list
    
    var_list <- lapply(idx.tree.vec, function(x)iter_var(var_list[[x]]$cov, var_list[[x]]$mean, theta_pred_list[, x, idx.iter],idx.iter))
  }
  
  time_diff <- Sys.time() - st_time
  print(time_diff)
  return(list(trees = trees.for.iter, 
              df_res = df_res.for.iter))
}

BART_single_step <- function(X, 
                             res_theta, 
                             U1, U2, 
                             prop_mu, prop_sigma, 
                             prior_list,
                             moves.prob = NULL, 
                             starting.tree = NULL, 
                             diag = FALSE, 
                             cont.unif = TRUE, 
                             include.split,
                             obs.per.term = 1, 
                             empty.count.lim = 10,
                             log_like_fun, 
                             log_prior_fun,
                             theta_param_1, theta_param_2,
                             var_param_1, var_param_2,
                             theta_var_list,
                             idx.tree){
  if(is.null(starting.tree)){
    rt_old <- generate_random_binary_tree_depth_free(1)
    rt_old <- assign_node_idx(rt_old)
    rt_old <- assign_split_rules(rt_old, X)
    rt_old <- assign_term_node_values_copula(rt_old) # check input
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
                                      old_tree = rt_old, 
                                      res_theta = res_theta, 
                                      X = X, 
                                      U1 = U1, U2 = U2, 
                                      prior_list = prior_list, 
                                      cont.unif = cont.unif, 
                                      include.split = include.split,
                                      prop_mu = prop_mu, prop_sigma = prop_sigma, 
                                      log_like_fun = log_like_fun, 
                                      log_prior_fun = log_prior_fun,
                                      theta_param_1 = theta_param_1, theta_param_2 = theta_param_2,
                                      theta_var_list= NULL,
                                      idx.tree = idx.tree)
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
  
  # new.tree <- new.tree.list$tree
  new.tree_mu.old <- new.tree.list$tree
  term_node_old <- get_terminal_nodes_idx(new.tree_mu.old)
  term_val_old <- sapply(term_node_old, function(i)unique(get_value_tree(new.tree_mu.old, get_obs_at_node(i, X, new.tree_mu.old,1,X))))
  
  
  new.tree_mu.new <- assign_term_node_values_cond_copula(tree_top = new.tree_mu.old,
                                                         X = X, 
                                                         prop_mu = prop_mu,
                                                         prop_sigma = prop_sigma,
                                                         theta_var_list = theta_var_list,
                                                         idx.tree = idx.tree)
  
  term_node_new <- get_terminal_nodes_idx(new.tree_mu.new)
  term_val_new <- sapply(term_node_new, function(i)unique(get_value_tree(new.tree_mu.new, get_obs_at_node(i, X, new.tree_mu.new,1,X))))
  
  term_val_prior_old <- (sum(sapply(term_val_old, function(x)log_prior_fun(x, theta_param_1, theta_param_2))))
  term_val_prior_new <- (sum(sapply(term_val_new, function(x)log_prior_fun(x, theta_param_1, theta_param_2))))
  
  ll_new <- cart_log_lik_copula(new.tree_mu.new, U1, U2, res_theta, X, log_like_fun = log_like_fun)
  ll_old <- cart_log_lik_copula(new.tree_mu.old, U1, U2, res_theta, X, log_like_fun = log_like_fun)
  
  logHR <- (ll_new - ll_old) + (term_val_prior_new - term_val_prior_old)
  
  if (!is.finite(logHR)) {
    logHR <- -Inf
  }
  
  if(runif(1) < exp(logHR)){
    new.tree <- new.tree_mu.new
  }else{
    new.tree <- new.tree_mu.old
  }
  
  term_nodes <- get_terminal_nodes_idx(new.tree)
  
  term_vals <- sapply(term_nodes, function(i)unique(get_value_tree(new.tree, get_obs_at_node(i, X, new.tree,1,X))))
  
  sig_2 <- rinvgamma(1, shape = (var_param_1 + length(term_vals)/2), scale = (var_param_2 + sum(term_vals^2)/2))
  
  matrix.res <- matrix(c(move.type, get_depth(rt_old), get_num_terminal_nodes(rt_old),
                         new.tree.list$prior.ratio, new.tree.list$lik.ratio, 
                         new.tree.list$trans.ratio, new.tree.list$acc.prob, 
                         new.tree.list$accepted, sig_2),nrow = 1)
  
  colnames(matrix.res) <- c('move', 'old.depth', 'old.nterm', 'prior.ratio', 
                            'lik.ratio', 'trans.ratio', 'acc.prob', 'acceptance','sig2')
  
  return(list("tree" = new.tree, "df_res" = data.frame(matrix.res)))
  
}


tree_step_copula <- function(move.type, 
                             old_tree, 
                             X, 
                             U1, U2, 
                             res_theta, 
                             prior_list, 
                             cont.unif = TRUE, 
                             include.split,
                             obs.per.term = 1, 
                             empty.count.lim = 10,
                             prop_mu, prop_sigma,  
                             log_like_fun, 
                             log_prior_fun,
                             theta_param_1, theta_param_2,
                             theta_var_list, idx.tree){
  empty.flag = TRUE
  empty.count = 0
  if (move.type %in% c("swap", "change", "grow", "prune")) {
    while (empty.flag & (empty.count <= empty.count.lim)) {
      
      move.list <- switch(move.type,
                          swap = swap_move_copula(old_tree),
                          
                          change = change_move_copula(
                            tree_top = old_tree, 
                            X = X, 
                            obs.per.term = obs.per.term,
                            cont.unif = cont.unif),
                          
                          grow = grow_move_copula(
                            tree_top = old_tree, 
                            X = X, 
                            obs.per.term = obs.per.term, 
                            cont.unif = cont.unif,
                            prop_sigma = prop_sigma,
                            theta_var_list = theta_var_list,
                            idx.tree = idx.tree
                          ),
                          
                          prune = prune_move_copula(
                            tree_top = old_tree, 
                            X = X, 
                            prop_mu = prop_mu, 
                            prop_sigma = prop_sigma,
                            theta_var_list = theta_var_list,
                            idx.tree = idx.tree
                          )
      )
      
      # Check if resulting tree has empty terminal nodes
      nobs.per.term <- vapply(
        get_terminal_nodes_idx(move.list$tree),
        \(x) nrow(get_obs_at_node(node.idx = x, X = X, tree_top = move.list$tree, X.orig = X)),
        0
      )
      
      if (sum(nobs.per.term == 0) > 0) {
        empty.count <- empty.count + 1
      } else {
        empty.flag <- FALSE
      }
      
      # Optional: exit loop early for 'prune' if count limit reached
      if (move.type == "prune" & empty.count > empty.count.lim) break
    }
  }
  
  
  if(empty.count > empty.count.lim){
    print('empty count exceeded')
    return('empty count exceeded')
  }
  
  acc.prob <- acceptance.prob.list_copula(move_list = move.list, 
                                          old.tree = old_tree,
                                          X = X,
                                          U1 = U1, U2 = U2, 
                                          res_theta,
                                          prior_input_list = prior_list,
                                          include.split = include.split,
                                          prop_mu = prop_mu, prop_sigma = prop_sigma,  
                                          log_like_fun = log_like_fun, 
                                          log_prior_fun = log_prior_fun,
                                          theta_param_1 = theta_param_1, theta_param_2 = theta_param_2,
                                          theta_var_list = theta_var_list,
                                          idx.tree = idx.tree)
  
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

assign_term_node_values_copula <- function(tree_top){ 
  
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  for(node.idx in term.node.idx){
    tree_top <- set_term_node_value_copula(node.idx, tree_top)
  }
  return(tree_top)
}

assign_term_node_values_cond_copula <- function(tree_top, 
                                                X, 
                                                prop_mu, prop_sigma, 
                                                theta_var_list, 
                                                idx.tree){
  
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  for(node.idx in term.node.idx){
    obs.at.node <- get_obs_at_node(node.idx = node.idx, X = X, tree_top = tree_top, X.orig = X)
    
    tree_top <- set_term_node_value_cond_copula(node.idx = node.idx, 
                                                tree_top = tree_top, 
                                                X = X, 
                                                prop_mu = mean(get_value_tree(tree_top, obs.at.node)),
                                                prop_sigma = prop_sigma,
                                                theta_var_list = theta_var_list,
                                                obs.id = as.numeric(rownames(obs.at.node)),
                                                idx.tree = idx.tree)
  }
  return(tree_top)
}

set_term_node_value_cond_copula <- function(node.idx, 
                                            tree_top, 
                                            X, 
                                            prop_mu, prop_sigma, 
                                            theta_var_list, 
                                            obs.id = NULL, 
                                            idx.tree){
  
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      
      tree_top$value = sample.cond.mu.copula(tree_top = tree_top, 
                                             node.idx = node.idx,
                                             X = X, 
                                             prop_mu = prop_mu, 
                                             prop_sigma = prop_sigma,
                                             theta_var_list = theta_var_list,
                                             obs.id = obs.id, 
                                             idx.tree = idx.tree)$val
      return(tree_top)
    } else {
      return(tree_top)  
    }
  }  else {
    tree.left <- set_term_node_value_cond_copula(node.idx, 
                                                 tree_top$left, 
                                                 X, 
                                                 prop_mu = prop_mu,
                                                 prop_sigma = prop_sigma,
                                                 theta_var_list = theta_var_list,
                                                 obs.id = obs.id, 
                                                 idx.tree = idx.tree)
    
    tree.right <- set_term_node_value_cond_copula(node.idx, 
                                                  tree_top$right, 
                                                  X, 
                                                  prop_mu = prop_mu, 
                                                  prop_sigma = prop_sigma,
                                                  theta_var_list = theta_var_list,
                                                  obs.id = obs.id, 
                                                  idx.tree = idx.tree)
    
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}

sample.cond.mu.copula <- function(tree_top = NULL, 
                                  node.idx = NULL, 
                                  X = NULL, 
                                  prop_mu, prop_sigma,
                                  theta_var_list,  
                                  obs.id = NULL, 
                                  idx.tree){
  
  if(is.null(obs.id)){
    obs.at.node <- get_obs_at_node(node.idx = node.idx, X = X, tree_top = tree_top, X.orig = X)
    obs.id <- as.numeric(rownames(obs.at.node))
    prop_mu = mean(get_value_tree(tree_top, obs.at.node))
  } 
  
  ##############################################################################
  if(!is.null(theta_var_list)){
    prop_cov_us <- theta_var_list[obs.id,obs.id, drop = F]
    
    n.vars <- ncol(prop_cov_us)
    
    prop_cov <- 2.4^2 / n.vars * (prop_cov_us + diag(0.005, nrow = n.vars))

    prop_sigma <- sqrt(sum(prop_cov)) / n.vars

  }
  
  proposal = rnorm(1, prop_mu, prop_sigma)
  
  return(list("val"= proposal, "mu" = prop_mu, "sigma" = prop_sigma))
}

grow_move_copula <- function(tree_top, 
                             X, 
                             cont.unif = TRUE, 
                             obs.per.term = 1, 
                             prop_sigma,
                             theta_var_list, 
                             idx.tree){ 
  
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  
  nobs.at.nodes <- vapply(term.node.idx, \(x) nrow(get_obs_at_node(node.idx = x, X = X, tree_top = tree_top, X.orig = X)),0)
  term.node.idx <- term.node.idx[nobs.at.nodes > 1]
  check.valid.idx <- vapply(term.node.idx, \(x) check_grow_valid(x, tree_top, X), TRUE)
  term.node.idx <- term.node.idx[check.valid.idx]
  
  if(length(term.node.idx) == 0){
    return(list(tree = tree_top, node.idx = NULL, move = 'grow', 
                valid.pred = NULL, 
                valid.split = NULL))
  } else if(length(term.node.idx) > 1){ 
    grow.idx <- sample(term.node.idx, 1, replace = FALSE)
  } else if(length(term.node.idx) == 1){
    grow.idx <- term.node.idx
  }
  
  new.cond.list <- gen_node_condition(grow.idx, tree_top, X, obs.per.term, cont.unif = cont.unif, for.grow = TRUE)
  obs.at.node <- new.cond.list$obs.at.node
  
  if(is.numeric(new.cond.list$cond$x.val)){
    idx.left <- obs.at.node[,new.cond.list$cond$x.idx] <= new.cond.list$cond$x.val
  } else {
    idx.left <- obs.at.node[, new.cond.list$cond$x.idx] %in% new.cond.list$cond$x.val
  }
  
  idx.right <- !idx.left
  
  term.node.value.left <- sample.cond.mu.copula(prop_mu = mean(get_value_tree(tree_top, obs.at.node[idx.left, new.cond.list$cond$x.idx, drop = FALSE])),
                                                prop_sigma = prop_sigma,
                                                theta_var_list = theta_var_list, 
                                                obs.id = idx.left, 
                                                idx.tree = idx.tree)
  
  term.node.value.right <- sample.cond.mu.copula(prop_mu = mean(get_value_tree(tree_top, obs.at.node[idx.right, new.cond.list$cond$x.idx, drop = FALSE])),
                                                 prop_sigma = prop_sigma,
                                                 theta_var_list = theta_var_list, 
                                                 obs.id = idx.right, 
                                                 idx.tree = idx.tree)
  
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

prune_move_copula <- function(tree_top, 
                              X, 
                              prop_mu, prop_sigma,
                              theta_var_list, 
                              idx.tree){ 
  
  prune.node.idx <- get_prune_idx(tree_top)
  
  if(length(prune.node.idx) == 1){
    prune.idx <- prune.node.idx
  } else{
    prune.idx <- sample(prune.node.idx, 1, replace = FALSE)
  }
  
  prune.value <- sample.cond.mu.copula(tree_top = tree_top, 
                                       node.idx = prune.idx,
                                       X = X, 
                                       prop_mu = prop_mu, 
                                       prop_sigma = prop_sigma,
                                       theta_var_list = theta_var_list, 
                                       obs.id = NULL, 
                                       idx.tree = idx.tree)
  
  tree_top_prune <- prune_terminal(node.idx = prune.idx, 
                                   tree_top = tree_top, 
                                   value.node = prune.value$val)
  
  tree_top_prune <- assign_node_idx(tree_top_prune)
  
  return(list(tree = tree_top_prune, move = 'prune', node.idx = prune.idx,
              prune.node.idx = prune.node.idx,
              list_prop = prune.value))
}

# perform swap move
swap_move_copula <- function(tree_top){
  
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
change_move_copula <- function(tree_top, 
                               X, 
                               obs.per.term = 2, 
                               cont.unif = TRUE){
  
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


acceptance.prob.list_copula <- function(move_list, 
                                        old.tree, 
                                        X, 
                                        U1, U2, 
                                        res_theta, 
                                        prior_input_list,
                                        include.split, 
                                        cont.unif = TRUE, 
                                        prop_mu, prop_sigma, 
                                        log_like_fun, 
                                        log_prior_fun,
                                        theta_param_1, theta_param_2,
                                        theta_var_list, 
                                        idx.tree){ 
  
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
    
    prob.old.to.new <- 1
    prob.new.to.old <- 1
    
    term_val_prior_new <- 1
    term_val_prior_old <- 1
    
    prop_dist_old_val <- 1
    prop_dist_new_val <- 1
    
  } else if(move_list$move == 'change'){
    
    prob.old.to.new <- 1
    prob.new.to.old <- 1
    
    term_val_prior_new <- 1
    term_val_prior_old <- 1
    
    prop_dist_old_val <- 1
    prop_dist_new_val <- 1
    
  } else if(move_list$move == 'grow'){
    
    prob.old.to.new <- transition.p.grow(move_list, include.split = include.split)
    prob.new.to.old <- transition.p.grow.reverse(move_list)
    
    ############################################################################
    term_val_old <- sample.cond.mu.copula(tree_top = old.tree, 
                                          node.idx = move_list$node.idx,
                                          X = X, 
                                          prop_mu = prop_mu, 
                                          theta_var_list = theta_var_list, 
                                          prop_sigma = prop_sigma,
                                          obs.id = NULL, 
                                          idx.tree = idx.tree)
    
    term_val_new_right <- move_list$prop_right
    term_val_new_left <- move_list$prop_left
    
    term_val_old_val <- mean(get_value_tree(old.tree, get_obs_at_node(node.idx = move_list$node.idx, X, old.tree, 1,X)))
    ############################################################################
    
    term_val_prior_old <- (sum(sapply(term_val_old_val, function(x)log_prior_fun(x, theta_param_1, theta_param_2))))
    term_val_prior_new <- (sum(sapply(c(term_val_new_left$val,term_val_new_right$val), function(x)log_prior_fun(x, theta_param_1, theta_param_2))))
    
    prop_dist_old_val <- sum(dnorm(term_val_old_val, mean = term_val_old$mu, sd = term_val_old$sigma, log = T))
    prop_dist_new_val <- sum(dnorm(c(term_val_new_left$val,term_val_new_right$val), mean = c(term_val_new_left$mu,term_val_new_right$mu), sd = c(term_val_new_left$sigma,term_val_new_right$sigma), log = T))
    
  } else if(move_list$move == 'prune'){
    
    prob.old.to.new <- transition.p.prune(move_list)
    prob.new.to.old <- transition.p.prune.reverse(move_list, old.tree, X, obs.per.term = 1,
                                                  include.split = include.split)
    
    ############################################################################
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
    
    term_val_old_left <- sample.cond.mu.copula(prop_mu = mean(get_value_tree(old.tree, obs.at.node[idx.left, new.cond.list$x.idx, drop = FALSE])),
                                               prop_sigma = prop_sigma,
                                               theta_var_list = theta_var_list,
                                               obs.id = idx.left, 
                                               idx.tree = idx.tree)
    
    term_val_old_right <- sample.cond.mu.copula(prop_mu = mean(get_value_tree(old.tree, obs.at.node[idx.right, new.cond.list$x.idx, drop = FALSE])),
                                                prop_sigma = prop_sigma,
                                                theta_var_list = theta_var_list,
                                                obs.id = idx.right, 
                                                idx.tree = idx.tree)
    
    term_val_new <- move_list$list_prop
    
    term_val_old_left_val <- mean(get_value_tree(old.tree, obs.at.node[idx.left, new.cond.list$x.idx, drop = FALSE]))
    term_val_old_right_val <- mean(get_value_tree(old.tree, obs.at.node[idx.right, new.cond.list$x.idx, drop = FALSE]))
    ############################################################################
    
    term_val_prior_old <- (sum(sapply(c(term_val_old_right_val,term_val_old_left_val), function(x)log_prior_fun(x, theta_param_1, theta_param_2))))
    term_val_prior_new <- (sum(sapply(term_val_new$val, function(x)log_prior_fun(x, theta_param_1, theta_param_2))))
    
    prop_dist_old_val <- sum(dnorm(c(term_val_old_right_val,term_val_old_left_val), mean = c(term_val_old_right$mu,term_val_old_left$mu), sd = c(term_val_old_right$sigma,term_val_old_left$sigma), log = T))
    prop_dist_new_val <- sum(dnorm(term_val_new$val, mean = term_val_new$mu, sd = term_val_new$sigma, log = T))
    
  } else {
    stop('Unknown move')
  }
  
  prior.ratio <- (prior.new/prior.old) * exp(term_val_prior_new - term_val_prior_old)
  lik.ratio <- exp(log.lik.new - log.lik.old)
  trans.ratio <- (prob.new.to.old/prob.old.to.new) * exp(prop_dist_old_val-prop_dist_new_val)
  
  acc.prob <- prior.ratio * lik.ratio * trans.ratio
  
  return(list(prior.ratio = prior.ratio, lik.ratio = lik.ratio, 
              trans.ratio = trans.ratio, alpha = min(1, acc.prob)))
}

cart_log_lik_copula <- function(tree_top, U1, U2, res_theta, X, log_like_fun = log_like_fun){
  
  tree.at.obs <- get_value_tree(tree_top, X)
  
  rho <- tree.at.obs + res_theta
  
  log.prob.obs <- log_like_fun(rho, U1, U2)
  
  return(log.prob.obs)
}

set_term_node_value_copula <- function(node.idx, tree_top){
  
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      tree_top$value = 0
      return(tree_top)
    } else {
      return(tree_top)  
    }
  }  else {
    
    tree.left <- set_term_node_value_copula(node.idx, tree_top$left)
    
    tree.right <- set_term_node_value_copula(node.idx, tree_top$right)
    
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}

# manual doesn't work for pseudo observations (u or v = 0 or 1)
# gauss copula 

loglik_gauss <- function(rho, u, v) {
  
  if(any(abs(rho)>=1)) return(-Inf)
  
  log_lik <- sum(log(BiCopPDF(u,v,1,rho)))
  return(log_lik)
}

# t copula

loglik_t <- function(rho, u, v, df = 3) {
  
  if(any(abs(rho)>=1)) return(-Inf)
  
  log_lik <- sum(log(BiCopPDF(u,v,2,rho, par2 = df)))
  
  return(log_lik)
}

# clayton copula

loglik_clayton <- function(theta, u, v) {
  if(any(abs(theta - 14)>=14)) return(-Inf)
  
  log_density <- log(BiCopPDF(u,v,3,theta))
  
  return(sum(log_density))
}

# gumbel copula

loglik_gumbel <- function(theta, u, v) {
  if(any(abs(theta - 9)>8)) return(-Inf)
  
  log_density <- log(BiCopPDF(u,v,4,theta))
  
  return(sum(log_density))
}

# gumbel copula

loglik_frank <- function(theta, u, v) {
  if(any(abs(theta - 0)>35) || any(theta == 0)) return(-Inf)
  
  log_density <- log(BiCopPDF(u,v,5,theta))
  
  return(sum(log_density))
}

logposterior <- function(rho, u, v, log_like_fun, log_prior_fun){
  return(log_like_fun(rho,u,v) + log_prior_fun(rho))
}

param_gauss <- function(tau) return(sin(tau*pi/2))
param_t <- function(tau) return(sin(tau*pi/2))
param_gumbel <- function(tau) return(1/(1-tau))
param_clayton <- function(tau) return((2*tau)/(1-tau))
param_frank <- function(tau) return(BiCopTau2Par(5,tau))

conv_diag <- function(post_data, n.burn, n.thin){
  
  post_mat <- (post_data)$nn
  
  post_mat_red <- na.omit((post_mat[n.burn+1:length(post_mat)])[c(rep(NA,(n.thin-1)), TRUE)])
  
  post_data_mcmc <- mcmc(post_mat_red, start = 1)
  
  auto_corr_sum <- autocorr.diag(post_data_mcmc)
  
  eff_size <- effectiveSize(post_data_mcmc)
  
  geweke_score <- geweke.diag(post_data_mcmc)
  
  
  summ_dat <- cbind(auto_corr_sum[2], 100*(eff_size/length(post_data_mcmc)), geweke_score$z)
  
  colnames(summ_dat) <- c("auto-correlation", "ess (%)", "geweke-score")
  
  return(summ_dat)
  
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
  ll <- lapply(1:length(bart_tree_list$df_res), \(idx) data.frame(idx = idx,
                                                                  nn = unlist(lapply(bart_tree_list$df_res[[idx]], 
                                                                                     \(x)as.factor(x$acceptance))),
                                                                  trees = factor( 1:length(bart_tree_list$df_res[[idx]]))))
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
  return(exp(x) + 0.0001)
}

link_gumbel <- function(x) {
  return(exp(x) + 1.0001)
}

link_frank <- function(x) {
  return(x)
}

var_adapt <- function(theta_var_list) {
  
  pred_mat <- t(theta_var_list)
  
  return(list("cov" = var(pred_mat), "mean" = colMeans(pred_mat)))
}

iter_var <- function(cov_mat, mean_vec_prev, theta_pred_iter, idx.iter) {
  
  last_row <- theta_pred_iter
  
  col_means <- mean_vec_prev + (last_row - mean_vec_prev) / idx.iter
  
  tcross_prev_mean <- tcrossprod(mean_vec_prev)
  tcross_col_means <- tcrossprod(col_means)
  tcross_last_row <- tcrossprod(last_row)
  
  coef_old <- (idx.iter - 2) / (idx.iter - 1)
  coef_update <- 1 / ((idx.iter - 1))
  
  var_at_burn_new <- coef_old * cov_mat +
    coef_update * ((idx.iter - 1) * tcross_prev_mean - 
                     idx.iter * tcross_col_means + 
                     tcross_last_row)
  
  return(list("cov" = var_at_burn_new, "mean" = col_means))
}