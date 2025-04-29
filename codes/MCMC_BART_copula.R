library(parallel)
library(mc2d)
multichain_MCMC_copula <- function(n.iter, n.chain,
                                   X, U1, U2, 
                                   mu, sigma, alpha_val, beta_val,
                                   log_nor_mu, log_nor_sigma, prior_type, cop_type,
                                   prior_list, 
                                   moves.prob = NULL, starting.tree = NULL,
                                   n.cores = 5,
                                   cont.unif = TRUE,
                                   include.split){
  chain.list <- mclapply(1:n.chain, 
                         \(x) MCMC_copula(n.iter = n.iter, 
                                          X = X, 
                                          U1 = U1,
                                          U2 = U2, 
                                          mu = mu, cop_type = cop_type, 
                                          sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, 
                                          log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type,
                                          prior_list = prior_list, 
                                          moves.prob = moves.prob, 
                                          starting.tree = starting.tree,
                                          cont.unif = cont.unif,
                                          include.split = include.split),
                         mc.cores = n.cores)
  tree_list <- lapply(chain.list, function(x) x$trees)
  tree_list_comb <- Reduce(c, tree_list)
  
  df_list <- lapply(1:n.chain, 
                    function(x) chain.list[[x]]$df.res %>% mutate(chain = x))
  df.res <- Reduce(rbind, df_list)
  return(list(trees = tree_list_comb, df.res = df.res))
}

MCMC_copula <- function(n.iter, X, U1, U2, mu, sigma, alpha_val, beta_val, prior_list,
                        log_nor_mu, log_nor_sigma, prior_type, cop_type,
                        moves.prob = NULL, starting.tree = NULL, 
                        diag = FALSE, cont.unif = TRUE, include.split){
  if(is.null(starting.tree)){
    rt_old <- generate_random_binary_tree_depth_free(1)
    rt_old <- assign_node_idx(rt_old)
    rt_old <- assign_split_rules(rt_old, X)
    rt_old <- assign_term_node_values_copula(rt_old, mu, sigma) # check input
  } else {
    rt_old <- starting.tree
  }
  tree_list <- list()
  matrix.res <- matrix(NA, nrow = n.iter, ncol = 8)
  colnames(matrix.res) <- c('move', 'old.depth', 'old.nterm', 'prior.ratio', 
                            'lik.ratio', 'trans.ratio', 'acc.prob', 'acceptance')
  for(i in 1:n.iter){
    cat(i,'\n')
    #debugonce(tree_step)
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
    cat(move.type, '\n')
    move.flag = TRUE
    while(move.flag){
      #calculate obs at nodes with last tree
      if(diag == TRUE){
        nobs_diag <- vapply(get_terminal_nodes_idx(rt_old),
                            \(x) nrow(get_obs_at_node(node.idx = x, X = X, 
                                                      tree_top = rt_old, X.orig = X)), 0)
        cat('Before move nobs at nodes:', nobs_diag,'\n')
      }
      
      new.tree.list <- tree_step_copula(move.type = move.type, # check input
                                        old_tree = rt_old, 
                                        X = X, U1 = U1, U2 = U2, 
                                        mu = mu, cop_type = cop_type, 
                                        sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, 
                                        log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type, 
                                        prior_list = prior_list, 
                                        cont.unif = cont.unif, 
                                        include.split = include.split)
      cat('is tree char?' , is.character(new.tree.list), '\n')
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
          cat('new move: ', move.type, '\n')  
        }
    }
    
    
    new.tree <- new.tree.list$tree
    new.tree <- assign_term_node_values_cond_copula(tree_top = new.tree, 
                                                    mu = mu,
                                                    sigma = sigma, alpha_val = alpha_val, beta_val = beta_val,
                                                    log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type,cop_type = cop_type,
                                                    X = X, U1 = U1, U2 = U2)
    
    
    tree_list[[i]] <- new.tree
    matrix.res[i,] <- c(move.type, get_depth(rt_old), get_num_terminal_nodes(rt_old),
                        new.tree.list$prior.ratio, new.tree.list$lik.ratio, 
                        new.tree.list$trans.ratio, new.tree.list$acc.prob, 
                        new.tree.list$accepted)
    rt_old <- new.tree
  }
  return(list(trees = tree_list, df.res = data.frame(matrix.res)))
}

tree_step_copula <- function(move.type, old_tree, X, U1, U2, mu, sigma, alpha_val, beta_val, 
                             log_nor_mu, log_nor_sigma, prior_type, cop_type, # check input
                             prior_list, cont.unif = TRUE, include.split,
                             obs.per.term = 1, empty.count.lim = 10){
  empty.flag = TRUE
  empty.count = 0
  if(move.type == 'swap'){
    while(empty.flag & (empty.count <= empty.count.lim)){
      move.list <- swap_move(old_tree)
      cat('move=',move.list$move, ',idx = ',move.list$node.idx[1,1],',',move.list$node.idx[1,2],'\n')
      # calculate the number of obs per term node
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
      
      move.list <- change_move(tree_top = old_tree, X = X, obs.per.term = obs.per.term, 
                               cont.unif = cont.unif)
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
      move.list <- grow_move_copula(tree_top = old_tree, X = X, U1 = U1, U2 = U2, # check input
                                    mu = mu, cop_type = cop_type,
                                    sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, 
                                    log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type,
                                    obs.per.term = obs.per.term,
                                    cont.unif = cont.unif)
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
      move.list <- prune_move_copula(tree_top = old_tree, # check input 
                                     mu = mu, cop_type = cop_type,
                                     sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, 
                                     log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type, 
                                     X = X, U1 = U1, U2 = U2)
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
                                          U2 = U2, cop_type = cop_type,
                                          prior_input_list = prior_list,
                                          include.split = include.split)
  acceptance <- runif(1) <= acc.prob$alpha
  print(acceptance)
  
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

assign_term_node_values_copula <- function(tree_top, mu, sigma){ # check input
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  for(node.idx in term.node.idx){
    tree_top <- set_term_node_value_copula(node.idx, tree_top, mu, sigma) # check input
  }
  return(tree_top)
}

assign_term_node_values_cond_copula <- function(tree_top, mu, sigma, alpha_val, beta_val, 
                                                log_nor_mu, log_nor_sigma, prior_type,cop_type, # check input
                                                X, U1, U2){
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  for(node.idx in term.node.idx){
    obs.at.node <- get_obs_at_node(node.idx = node.idx, X = X, tree_top = tree_top, X.orig = X)
    U1.at.node <- U1[as.numeric(rownames(obs.at.node))]
    U2.at.node <- U2[as.numeric(rownames(obs.at.node))]
    tree_top <- set_term_node_value_cond_copula(node.idx = node.idx, tree_top = tree_top, # check input
                                                mu = mean(get_value_tree(tree_top, obs.at.node)), cop_type = cop_type,
                                                sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, 
                                                log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type,
                                                X = X, U1 =U1, U2=U2, U1.at.node = U1.at.node, 
                                                U2.at.node = U2.at.node)
  }
  return(tree_top)
}

set_term_node_value_cond_copula <- function(node.idx, tree_top, mu, sigma, alpha_val, beta_val, 
                                            log_nor_mu, log_nor_sigma, prior_type, cop_type, # check input
                                            X, U1, U2, 
                                            U1.at.node=NULL, U2.at.node=NULL, binary = FALSE){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      tree_top$value = sample.cond.mu.copula(tree_top = tree_top, # check input
                                             node.idx = node.idx, 
                                             mu = mu, cop_type = cop_type,
                                             sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, 
                                             log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type,
                                             X = X, 
                                             U1 = U1, U2 = U2, U1.at.node=U1.at.node, 
                                             U2.at.node=U2.at.node)
      return(tree_top)
    } else {
      return(tree_top)  
    }
  }  else {
    tree.left <- set_term_node_value_cond_copula(node.idx, tree_top$left, mu, sigma, alpha_val = alpha_val, beta_val = beta_val, cop_type = cop_type,
                                                 log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type, X, U1, U2, 
                                                 U1.at.node, U2.at.node)
    tree.right <- set_term_node_value_cond_copula(node.idx, tree_top$right, mu, sigma, alpha_val = alpha_val, beta_val = beta_val, cop_type = cop_type,
                                                  log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type, X, U1, U2, 
                                                  U1.at.node, U2.at.node)
    
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}

# sample conditional mu assuming binary observations and beta prior
sample.cond.mu.copula <- function(tree_top = NULL, 
                                  node.idx = NULL, 
                                  mu, 
                                  sigma, cop_type, alpha_val, beta_val, 
                                  log_nor_mu, log_nor_sigma, prior_type, 
                                  X = NULL, 
                                  U1 = NULL, 
                                  U2 = NULL, 
                                  U1.at.node = NULL,
                                  U2.at.node = NULL){
  if(is.null(U1.at.node)){
    obs.at.node <- get_obs_at_node(node.idx = node.idx, X = X, tree_top = tree_top, X.orig = X)
    U1.at.node <- U1[as.numeric(rownames(obs.at.node))]
    U2.at.node <- U2[as.numeric(rownames(obs.at.node))]
    mu = mean(get_value_tree(tree_top, obs.at.node))
  } 
  
  ##############################################################################
  proposal = rnorm(1, mu, sigma)
  print(proposal)
  if(cop_type == "Gauss"){
    if((proposal <= -1)||(proposal >= 1)) HR = 0 else # Q needs to be >0
      # Hastings ratio of the proposal
      HR = exp(logposterior(U1.at.node, U2.at.node, rho = proposal, alpha_val = alpha_val, cop_type = cop_type, 
                            beta_val = beta_val, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type) -
                 logposterior(U1.at.node, U2.at.node, rho = mu, alpha_val = alpha_val, cop_type = cop_type, 
                              beta_val = beta_val, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type))
  } else if(cop_type == "Frank"){
    HR = exp(logposterior(U1.at.node, U2.at.node, rho = proposal, alpha_val = alpha_val, cop_type = cop_type, 
                          beta_val = beta_val, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type) -
               logposterior(U1.at.node, U2.at.node, rho = mu, alpha_val = alpha_val, cop_type = cop_type, 
                            beta_val = beta_val, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type))
  }else if(cop_type == "Clayton"){
    if(proposal<=0) HR = 0 else{
      HR = exp(logposterior(U1.at.node, U2.at.node, rho = proposal, alpha_val = alpha_val, cop_type = cop_type, 
                            beta_val = beta_val, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type) -
                 logposterior(U1.at.node, U2.at.node, rho = mu, alpha_val = alpha_val, cop_type = cop_type, 
                              beta_val = beta_val, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type))
    }
  }
  
  
  
  
  if (runif(1) < HR){ 
    new_mu = proposal
    # if proposal is rejected, keep the values from the previous iteration
  }else{
    new_mu = mu
  }
  
  return(new_mu)
}

grow_move_copula <- function(tree_top, X, U1, U2, mu, sigma, alpha_val, beta_val, log_nor_mu, log_nor_sigma, prior_type, cop_type, cont.unif = TRUE, obs.per.term = 1){ # check input
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  # this is because if a terminal node has only 1 obs associated cannot be grow
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
  term.node.value.left <- sample.cond.mu.copula(U1.at.node = U1.at.left,
                                                U2.at.node = U2.at.left, cop_type = cop_type,
                                                mu = mean(get_value_tree(tree_top, obs.at.node[idx.left, new.cond.list$cond$x.idx, drop = FALSE])),
                                                sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type)
  term.node.value.right <- sample.cond.mu.copula(U1.at.node = U1.at.right,
                                                 U2.at.node = U2.at.right, cop_type = cop_type,
                                                 mu = mean(get_value_tree(tree_top, obs.at.node[idx.right, new.cond.list$cond$x.idx, drop = FALSE])),
                                                 sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type)
  
  tree_top_grow <- grow_terminal(node.idx = grow.idx, 
                                 tree_top = tree_top, 
                                 new.cond = new.cond.list$cond, 
                                 value.l = term.node.value.left,
                                 value.r = term.node.value.right)
  tree_top_grow <- assign_node_idx(tree_top_grow)
  return(list(tree = tree_top_grow, node.idx = grow.idx, move = 'grow', 
              valid.pred = new.cond.list$valid.pred, 
              valid.split = new.cond.list$valid.split))
}

prune_move_copula <- function(tree_top, mu, sigma, alpha_val, beta_val, log_nor_mu, log_nor_sigma, prior_type, cop_type, X, U1, U2){ # check input
  prune.node.idx <- get_prune_idx(tree_top)
  if(length(prune.node.idx) == 1){
    prune.idx <- prune.node.idx
  } else{
    prune.idx <- sample(prune.node.idx, 1, replace = FALSE)
  }
  prune.value <- sample.cond.mu.copula(tree_top = tree_top, # check input
                                       node.idx = prune.idx, cop_type = cop_type,
                                       mu = mu,
                                       sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, 
                                       log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type,
                                       X = X,
                                       U1 = U1,
                                       U2 = U2)
  tree_top_prune <- prune_terminal(node.idx = prune.idx, 
                                   tree_top = tree_top, 
                                   value.node = prune.value)
  tree_top_prune <- assign_node_idx(tree_top_prune)
  return(list(tree = tree_top_prune, move = 'prune', node.idx = prune.idx,
              prune.node.idx = prune.node.idx))
}

acceptance.prob.list_copula <- function(move_list, old.tree, X, U1, U2, prior_input_list, cop_type, 
                                        include.split, cont.unif = TRUE){ # check input
  if(identical(move_list$tree, old.tree)){
    return(list(prior.ratio = 1, lik.ratio = 1, 
                trans.ratio = 1, alpha = 1))
  }
  
  prior.tree.old <- prior_input_list$fun(old.tree, prior_input_list$param[1], prior_input_list$param[1])
  prior.old <- prior.tree.old*prior.split.rule(old.tree, X, cont.unif = cont.unif)
  
  prior.tree.new <- prior_input_list$fun(move_list$tree, prior_input_list$param[1], prior_input_list$param[1])
  prior.new <- prior.tree.new*prior.split.rule(move_list$tree, X, cont.unif = cont.unif)
  # likelihood 
  log.lik.old <- cart_log_lik_copula(tree_top = old.tree, U1 = U1, U2 = U2, X = X, cop_type = cop_type)
  log.lik.new <- cart_log_lik_copula(tree_top = move_list$tree, U1 = U1, U2 = U2, X = X, cop_type = cop_type)
  
  #transition probabilities
  if(move_list$move == 'swap'){
    prob.old.to.new <- 1#transition.p.swap(old.tree)
    prob.new.to.old <- 1#transition.p.swap(move_list$tree)
  } else if(move_list$move == 'change'){
    prob.old.to.new <- 1#transition.p.change(move_list)
    prob.new.to.old <- 1#transition.p.change.reverse(move_list, X, old.tree)
  } else if(move_list$move == 'grow'){
    prob.old.to.new <- transition.p.grow(move_list, include.split = include.split)
    prob.new.to.old <- transition.p.grow.reverse(move_list)
  } else if(move_list$move == 'prune'){
    prob.old.to.new <- transition.p.prune(move_list)
    prob.new.to.old <- transition.p.prune.reverse(move_list, old.tree, X, obs.per.term = 1,
                                                  include.split = include.split)
  } else {
    stop('Unknown move')
  }
  prior.ratio <- (prior.new/prior.old)
  lik.ratio <- exp(log.lik.new - log.lik.old)
  trans.ratio <- (prob.new.to.old/prob.old.to.new)
  acc.prob <- prior.ratio*lik.ratio*trans.ratio 
  print(lik.ratio)
  return(list(prior.ratio = prior.ratio, lik.ratio = lik.ratio, 
              trans.ratio = trans.ratio, alpha = min(1, acc.prob)))
}

cart_log_lik_copula <- function(tree_top, U1, U2, X, cop_type){ # check input
  tree.at.obs <- get_value_tree(tree_top, X)
  print(length(tree.at.obs))
  if(cop_type == "Gauss"){
    log.prob.obs <- gaussian_copula_loglik(tree.at.obs, U1, U2)
  }else if(cop_type == "Frank"){
    log.prob.obs <- loglik_frank(tree.at.obs, U1, U2)
  }else if(cop_type == "Clayton"){
    log.prob.obs <- loglik_clayton(tree.at.obs, U1, U2)
  }
  
  return(sum(log.prob.obs))
}

set_term_node_value_copula <- function(node.idx, tree_top, mu, sigma){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      tree_top$value = rnorm(1, mean = mu, sd = sigma)
      return(tree_top)
    } else {
      return(tree_top)  
    }
  }  else {
    tree.left <- set_term_node_value_copula(node.idx, tree_top$left, mu, sigma)
    tree.right <- set_term_node_value_copula(node.idx, tree_top$right, mu, sigma)
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}

gaussian_copula_loglik <- function(rho, u, v) {
  z1 <- qnorm(u)
  z2 <- qnorm(v)
  
  # Number of observations
  n <- length(z1)
  
  # Compute log-likelihood
  log_det = log(1 - rho^2)
  quad_form = (z1^2 - 2 * rho * z1 * z2 + z2^2) / (1 - rho^2)
  
  log_lik = -0.5 * sum(log_det + quad_form)
  return(log_lik)
}

loglik_frank <- function(theta, u, v) {
  
  theta = pmin(theta, 35)
  theta = pmax(theta, -35)
  densty = BiCopPDF(u, v, par = theta, family = 5)
  return(sum(log(densty)))
}

loglik_clayton <- function(theta, u, v) {
  
  theta = pmin(theta, 28)
  theta = pmax(theta, 1e-10)
  densty = BiCopPDF(u, v, par = theta, family = 3)
  return(sum(log(densty)))
}


logprior_unif <- function(rho, bound = 1, alpha_val, beta_val) {
  return((alpha_val - 1)*log(bound + rho) + (beta_val - 1)*log(bound - rho))
}

logprior_inv_gamma <- function(p, alpha_val, beta_val) {
  # Ensure p is within (-1, 1)
  if (any(abs(p) >= 1)) {
    return(rep(0, length(p)))
  }
  
  y <- -log(1 - p^2)
  
  dens <- (alpha_val * log(beta_val) + log(abs(p)) + (-alpha_val - 1) * log(y) - 
             beta_val / y - log(1 - p^2) - log(gamma(alpha_val))) 
  
  return(dens)
}

logprior_log_normal <- function(p, log_nor_mu, log_nor_sigma) {
  if (any(abs(p) >= 1)) {
    return(rep(0, length(p)))
  }
  
  y <- -log(1 - p^2)
  log_y <- log(y)
  
  # Density of p
  density <- log(2 * abs(p)) - (log(1 - p^2) + log(y)) +(-((log_y - log_nor_mu)^2) / (2 * log_nor_sigma^2))
  
  return(density)
}

logposterior <- function(rho, u, v, alpha_val, beta_val, log_nor_mu, log_nor_sigma, prior_type = "B", cop_type){
  if(cop_type == "Gauss"){
    if(prior_type == "B"){
      return(gaussian_copula_loglik(rho, u, v) + logprior_unif(rho, alpha_val = alpha_val, beta_val = beta_val))
    }else{
      if(prior_type == "LN"){
        return(gaussian_copula_loglik(rho, u, v) + logprior_log_normal(rho, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma))
      }else{
        return(gaussian_copula_loglik(rho, u, v) + logprior_inv_gamma(rho, alpha_val = alpha_val, beta_val = beta_val))
      }
    }
  }else if(cop_type == "Frank"){
    if(prior_type == "Normal"){
      return(loglik_frank(rho, u, v) + log(dnorm(rho, mean = log_nor_mu, sd = log_nor_sigma)))
    }else{
      if(prior_type == "T"){
        return(loglik_frank(rho, u, v) + log(dt(rho, ncp = log_nor_mu, df = 1/log_nor_sigma)))
      }else{
        if(prior_type == "Unif"){
          return(loglik_frank(rho, u, v) + log(dunif(rho, min = -35, max = 35)))
        }else{
          return(loglik_frank(rho, u, v) + logprior_unif(rho, bound = 35+1e-3, alpha_val, beta_val))
        }
      }
    }
  }else if(cop_type == "Clayton"){
    if(prior_type == "IG"){
      return(loglik_clayton(rho, u, v) + log(dgamma(rho, shape = alpha_val, rate = beta_val)))
    }else{
      if(prior_type == "LN"){
        return(loglik_clayton(rho, u, v) + log(dlnorm(rho, meanlog = log_nor_mu, sdlog = log_nor_sigma)))
      }
    }
  }
}