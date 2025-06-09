library(parallel)
library(mc2d)
library(pracma)
multichain_MCMC_copula <- function(n.iter, n.chain,
                                   X, U1, U2, 
                                   mu, sigma, alpha_val, beta_val,
                                   log_nor_mu, log_nor_sigma, prior_type, cop_type,
                                   prior_list, 
                                   moves.prob = NULL, starting.tree = NULL,
                                   n.cores = min(n.chain, 5),
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
    rt_old <- assign_term_node_values_copula(rt_old, mu, sigma, cop_type) # check input
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
    
    
    # sample mu from conditional 
    
    new.tree_mu.old <- new.tree.list$tree
    term_node_old <- get_terminal_nodes_idx(new.tree_mu.old)
    term_val_old <- sapply(term_node_old, function(i)unique(get_value_tree(new.tree_mu.old, get_obs_at_node(i, X, new.tree_mu.old,1,X))))
    
    
    new.tree_mu.new <- assign_term_node_values_cond_copula(tree_top = new.tree_mu.old,
                                                          mu = mu,
                                                          sigma = sigma, alpha_val = alpha_val, beta_val = beta_val,
                                                          log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type,cop_type = cop_type,
                                                          X = X, U1 = U1, U2 = U2)
    term_node_new <- get_terminal_nodes_idx(new.tree_mu.new)
    term_val_new <- sapply(term_node_new, function(i)unique(get_value_tree(new.tree_mu.new, get_obs_at_node(i, X, new.tree_mu.new,1,X))))
    
    if(cop_type == "Gauss"){
      if(prior_type == "B"){
        term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)log_prior_tbeta(x, alpha_val = alpha_val, beta_val = beta_val))))
        term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)log_prior_tbeta(x, alpha_val = alpha_val, beta_val = beta_val))))
      }else{
        if(prior_type == "LN"){
          term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)logprior_log_normal(x, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma))))
          term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)logprior_log_normal(x, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma))))
        }else{
          term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)logprior_inv_gamma(x, alpha_val = alpha_val, beta_val = beta_val))))
          term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)logprior_inv_gamma(x, alpha_val = alpha_val, beta_val = beta_val))))
        }
      }
    } else if(cop_type == "Clayton"){
      if(prior_type == "IG"){
        term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)log(dgamma(x, shape = alpha_val, rate = beta_val)))))
        term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)log(dgamma(x, shape = alpha_val, rate = beta_val)))))
      }else if(prior_type == "LN"){
        term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)log(dlnorm(x, meanlog = log_nor_mu, sdlog = log_nor_sigma)))))
        term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)log(dlnorm(x, meanlog = log_nor_mu, sdlog = log_nor_sigma)))))
      }
    }else if(cop_type == "gumbel"){
      if(prior_type == "IG"){
        term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)log(dgamma(x, shape = alpha_val, rate = beta_val)))))
        term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)log(dgamma(x, shape = alpha_val, rate = beta_val)))))
      }else{
        if(prior_type == "LN"){
          term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)log(dlnorm(x, meanlog = log_nor_mu, sdlog = log_nor_sigma)))))
          term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)log(dlnorm(x, meanlog = log_nor_mu, sdlog = log_nor_sigma)))))
        }
      }
    }else if(cop_type == "t"){
      if(prior_type == "B"){
        term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)log_prior_tbeta(x, alpha_val = alpha_val, beta_val = beta_val))))
        term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)log_prior_tbeta(x, alpha_val = alpha_val, beta_val = beta_val))))
      }else{
        if(prior_type == "LN"){
          term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)logprior_log_normal(x, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma))))
          term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)logprior_log_normal(x, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma))))
        }else{
          term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)logprior_inv_gamma(x, alpha_val = alpha_val, beta_val = beta_val))))
          term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)logprior_inv_gamma(x, alpha_val = alpha_val, beta_val = beta_val))))
        }
      }
    }
    
    HR <- exp(cart_log_lik_copula(new.tree_mu.new, U1, U2, X, cop_type) - cart_log_lik_copula(new.tree_mu.old, U1, U2, X, cop_type)) * (term_val_prior_new / term_val_prior_old)
    
    if(runif(1)<HR){
      new.tree <- new.tree_mu.new
    }else{
      new.tree <- new.tree_mu.old
    }
    
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
      move.list <- swap_move_copula(old_tree, X, U1, U2, mu = mean(get_value_tree(tree_top, obs.at.node_1)), cop_type = cop_type,
                                    sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, 
                                    log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type)
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
      
      move.list <- change_move_copula(tree_top = old_tree, X = X, obs.per.term = obs.per.term, 
                               cont.unif = cont.unif, U1, U2, mu = mu, cop_type = cop_type,
                               sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, 
                               log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type)
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
                                          mu = mu,
                                          sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, 
                                          log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type,
                                          prior_input_list = prior_list,
                                          include.split = include.split)
  acceptance <- runif(1) <= acc.prob$alpha
  
  print("tatha")
  print(acceptance)
  print("tatha")
  
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

assign_term_node_values_copula <- function(tree_top, mu, sigma, cop_type){ # check input
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  for(node.idx in term.node.idx){
    tree_top <- set_term_node_value_copula(node.idx, tree_top, mu, sigma, cop_type) # check input
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

assign_term_node_proposal_cond_copula <- function(tree_top, mu, sigma, alpha_val, beta_val, 
                                                log_nor_mu, log_nor_sigma, prior_type,cop_type, # check input
                                                X, U1, U2){
  term.node.idx <- get_terminal_nodes_idx(tree_top)
  
  prop_par <- c()
  for(node.idx in term.node.idx){
    obs.at.node <- get_obs_at_node(node.idx = node.idx, X = X, tree_top = tree_top, X.orig = X)
    U1.at.node <- U1[as.numeric(rownames(obs.at.node))]
    U2.at.node <- U2[as.numeric(rownames(obs.at.node))]
    c <- sample.cond.mu.copula(tree_top = tree_top, # check input
                               node.idx = node.idx, 
                               mu = mu, cop_type = cop_type,
                               sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, 
                               log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type,
                               X = X, 
                               U1 = U1, U2 = U2, U1.at.node=U1.at.node, 
                               U2.at.node=U2.at.node)
    
    prop_par <- rbind(prop_par, unlist(c)[-1])
  }
  return(prop_par)
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
                                             U2.at.node=U2.at.node)$val
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
  log_post <- function(x)logposterior(x, U1.at.node, U2.at.node, alpha_val = alpha_val, 
                                      beta_val= beta_val, log_nor_mu= log_nor_mu, log_nor_sigma= log_nor_sigma, 
                                      prior_type = prior_type, cop_type = cop_type)
  
  optim_mu = optim(0, function(x)-log_post(x))
  
  prop_mu = optim_mu$par
  
  prop_sd2_inv = fderiv(function(x)-log_post(x), prop_mu, n=2)
  
  prop_sd = sqrt(1/prop_sd2_inv)
  
  proposal = rnorm(1, prop_mu, prop_sd)
  
  print("tatha2")
  print(proposal)
  print("tatha2")
  # HR = exp(logposterior(U1.at.node, U2.at.node, rho = proposal, alpha_val = alpha_val, cop_type = cop_type, 
  #                                                 beta_val = beta_val, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type) -
  #                                      logposterior(U1.at.node, U2.at.node, rho = mu, alpha_val = alpha_val, cop_type = cop_type,
  #                                                   beta_val = beta_val, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type))
  # 
  # 
  # 
  # if (runif(1) < HR){
  #   new_mu = proposal
  #   # if proposal is rejected, keep the values from the previous iteration
  # }else{
  #   new_mu = mu
  # }
  return(list("val"= proposal, "mu" = prop_mu, "sigma" = prop_sd))
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
                                                sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, log_nor_mu = log_nor_mu, 
                                                log_nor_sigma = log_nor_sigma, prior_type = prior_type)$val
  term.node.value.right <- sample.cond.mu.copula(U1.at.node = U1.at.right,
                                                 U2.at.node = U2.at.right, cop_type = cop_type,
                                                 mu = mean(get_value_tree(tree_top, obs.at.node[idx.right, new.cond.list$cond$x.idx, drop = FALSE])),
                                                 sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, log_nor_mu = log_nor_mu, 
                                                 log_nor_sigma = log_nor_sigma, prior_type = prior_type)$val
  
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
                                       U2 = U2)$val
  tree_top_prune <- prune_terminal(node.idx = prune.idx, 
                                   tree_top = tree_top, 
                                   value.node = prune.value)
  tree_top_prune <- assign_node_idx(tree_top_prune)
  return(list(tree = tree_top_prune, move = 'prune', node.idx = prune.idx,
              prune.node.idx = prune.node.idx))
}

# perform swap move
swap_move_copula <- function(tree_top, X, U1, U2, mu, cop_type,
                             sigma, alpha_val, beta_val, 
                             log_nor_mu, log_nor_sigma, prior_type){
  swap.node.idx <- get_valid_swap_idx(tree_top)
  if(nrow(swap.node.idx) == 0){
    return(list(tree = tree_top, move = 'swap', valid.pred = NULL, valid.split = NULL))
  } else {
    swap.idx.row <- sample(1:nrow(swap.node.idx), 1)
    swap.idx <- swap.node.idx[swap.idx.row, ]
    swap.tree_top <- swap_node_condition(swap.idx[1,1], swap.idx[1,2], tree_top)
    
    # obs.at.node_1 <- get_obs_at_node(node.idx = swap.idx[1,1], X = X, tree_top = swap.tree_top, X.orig = X)
    # U1.at.node_1 <- U1[as.numeric(rownames(obs.at.node_1))]
    # U2.at.node_1 <- U2[as.numeric(rownames(obs.at.node_1))]
    # 
    # swap.tree_top <- set_term_node_value_cond_copula(node.idx = swap.idx[1,1], tree_top = swap.tree_top, # check input
    #                                 mu = mean(get_value_tree(swap.tree_top, obs.at.node_1)), cop_type = cop_type,
    #                                 sigma = sigma, alpha_val = alpha_val, beta_val = beta_val,
    #                                 log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type,
    #                                 X = X, U1 =U1, U2=U2, U1.at.node = U1.at.node_1,
    #                                 U2.at.node = U2.at.node_1)
    # 
    # obs.at.node_2 <- get_obs_at_node(node.idx = swap.idx[1,2], X = X, tree_top = swap.tree_top, X.orig = X)
    # U1.at.node_2 <- U1[as.numeric(rownames(obs.at.node_2))]
    # U2.at.node_2 <- U2[as.numeric(rownames(obs.at.node_2))]
    # 
    # swap.tree_top <- set_term_node_value_cond_copula(node.idx = swap.idx[1,2], tree_top = swap.tree_top, # check input
    #                                                  mu = mean(get_value_tree(swap.tree_top, obs.at.node_2)), cop_type = cop_type,
    #                                                  sigma = sigma, alpha_val = alpha_val, beta_val = beta_val,
    #                                                  log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type,
    #                                                  X = X, U1 =U1, U2=U2, U1.at.node = U1.at.node_2,
    #                                                  U2.at.node = U2.at.node_2)
    
    return(list(tree = swap.tree_top, move = 'swap', node.idx = swap.idx, valid.pred = NULL, valid.split = NULL))
  }
}

# perform change move
change_move_copula <- function(tree_top, X, obs.per.term = 2, cont.unif = TRUE, U1,U2,
                               mu, cop_type,
                               sigma, alpha_val, beta_val, 
                               log_nor_mu, log_nor_sigma, prior_type){
  internal.node.idx <- get_internal_nodes_idx(tree_top)
  change.idx <- sample(internal.node.idx, 1)
  new.cond <- gen_node_condition(change.idx, tree_top, X, obs.per.term, cont.unif)
  change.tree_top <- set_node_condition(change.idx, tree_top, new.cond$cond)
  
  change.tree_top <- set_term_node_value_cond_copula(node.idx = change.idx, tree_top = change.tree_top, # check input
                                                     mu = mu, cop_type = cop_type,
                                                     sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, 
                                                     log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type,
                                                     X = X, U1 =U1, U2=U2, U1.at.node = NULL, 
                                                     U2.at.node = NULL)
  
  return(list(tree = change.tree_top, move = 'change', 
              node.idx = change.idx, 
              valid.pred = new.cond$valid.pred, 
              valid.split = new.cond$valid.split,
              n.left = new.cond$n.left,
              n.right = new.cond$n.right))
}


acceptance.prob.list_copula <- function(move_list, old.tree, X, U1, U2, prior_input_list, cop_type, 
                                        mu,
                                        sigma, alpha_val, beta_val, 
                                        log_nor_mu, log_nor_sigma, prior_type,
                                        include.split, cont.unif = TRUE){ # check input
  if(identical(move_list$tree, old.tree)){
    return(list(prior.ratio = 1, lik.ratio = 1, 
                trans.ratio = 1, alpha = 1))
  }
  
  term_node_old <- get_terminal_nodes_idx(old.tree)
  term_node_new <- get_terminal_nodes_idx(move_list$tree)
  
  term_val_old <- sapply(term_node_old, function(i)unique(get_value_tree(old.tree, get_obs_at_node(i, X, old.tree,1,X))))
  term_val_new <- sapply(term_node_new, function(i)unique(get_value_tree(move_list$tree, get_obs_at_node(i, X, move_list$tree,1,X))))
  
  prop_dist_old <- assign_term_node_proposal_cond_copula(old.tree, mu = mu, sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, 
                                                         log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type, cop_type = cop_type, # check input
                                                         X, U1, U2)
  
  prop_dist_new <- assign_term_node_proposal_cond_copula(move_list$tree, mu = mu, sigma = sigma, alpha_val = alpha_val, beta_val = beta_val, 
                                                         log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma, prior_type = prior_type, cop_type = cop_type, # check input
                                                         X, U1, U2)
  
  prop_dist_old_val <- prod(dnorm(term_val_old, mean = prop_dist_old[,1], sd = prop_dist_old[,2]))
  prop_dist_new_val <- prod(dnorm(term_val_new, mean = prop_dist_new[,1], sd = prop_dist_new[,2]))
  
  
  if(move_list$move == "grow"){
    list2env(list("old_term" = prop_dist_old, "new_term" = prop_dist_new), envir = .GlobalEnv)
  }
  
  if(move_list$move == "grow"){
    list2env(list("old_term_val" = term_val_old, "new_term_val" = term_val_new), envir = .GlobalEnv)
  }
  
  if(cop_type == "Gauss"){
    log.lik.new <- cart_log_lik_copula(tree_top = move_list$tree, U1 = U1, U2 = U2, X = X, cop_type = cop_type)
    if(prior_type == "B"){
      term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)log_prior_tbeta(x, alpha_val = alpha_val, beta_val = beta_val))))
      term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)log_prior_tbeta(x, alpha_val = alpha_val, beta_val = beta_val))))
    }else{
      if(prior_type == "LN"){
        term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)logprior_log_normal(x, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma))))
        term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)logprior_log_normal(x, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma))))
      }else{
        term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)logprior_inv_gamma(x, alpha_val = alpha_val, beta_val = beta_val))))
        term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)logprior_inv_gamma(x, alpha_val = alpha_val, beta_val = beta_val))))
      }
    }
  } else if(cop_type == "Clayton"){
    log.lik.new <- cart_log_lik_copula(tree_top = move_list$tree, U1 = U1, U2 = U2, X = X, cop_type = cop_type)
    if(prior_type == "IG"){
      term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)log(dgamma(x, shape = alpha_val, rate = beta_val)))))
      term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)log(dgamma(x, shape = alpha_val, rate = beta_val)))))
    }else if(prior_type == "LN"){
      term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)log(dlnorm(x, meanlog = log_nor_mu, sdlog = log_nor_sigma)))))
      term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)log(dlnorm(x, meanlog = log_nor_mu, sdlog = log_nor_sigma)))))
    }
  }else if(cop_type == "gumbel"){
    log.lik.new <- cart_log_lik_copula(tree_top = move_list$tree, U1 = U1, U2 = U2, X = X, cop_type = cop_type)
    if(prior_type == "IG"){
      term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)log(dgamma(x, shape = alpha_val, rate = beta_val)))))
      term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)log(dgamma(x, shape = alpha_val, rate = beta_val)))))
    }else{
      if(prior_type == "LN"){
        term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)log(dlnorm(x, meanlog = log_nor_mu, sdlog = log_nor_sigma)))))
        term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)log(dlnorm(x, meanlog = log_nor_mu, sdlog = log_nor_sigma)))))
      }
    }
  }else if(cop_type == "t"){
    log.lik.new <- cart_log_lik_copula(tree_top = move_list$tree, U1 = U1, U2 = U2, X = X, cop_type = cop_type)
    if(prior_type == "B"){
      term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)log_prior_tbeta(x, alpha_val = alpha_val, beta_val = beta_val))))
      term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)log_prior_tbeta(x, alpha_val = alpha_val, beta_val = beta_val))))
    }else{
      if(prior_type == "LN"){
        term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)logprior_log_normal(x, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma))))
        term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)logprior_log_normal(x, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma))))
      }else{
        term_val_prior_old <- exp(sum(sapply(term_val_old, function(x)logprior_inv_gamma(x, alpha_val = alpha_val, beta_val = beta_val))))
        term_val_prior_new <- exp(sum(sapply(term_val_new, function(x)logprior_inv_gamma(x, alpha_val = alpha_val, beta_val = beta_val))))
      }
    }
  }
  
  prior.tree.old <- prior_input_list$fun(old.tree, prior_input_list$param[1], prior_input_list$param[2])
  prior.old <- prior.tree.old*prior.split.rule(old.tree, X, cont.unif = cont.unif) * term_val_prior_old
  
  prior.tree.new <- prior_input_list$fun(move_list$tree, prior_input_list$param[1], prior_input_list$param[2])
  prior.new <- prior.tree.new*prior.split.rule(move_list$tree, X, cont.unif = cont.unif) * term_val_prior_new
  # likelihood 
  log.lik.old <- cart_log_lik_copula(tree_top = old.tree, U1 = U1, U2 = U2, X = X, cop_type = cop_type)
  
  
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
  trans.ratio <- (prob.new.to.old/prob.old.to.new)*(prop_dist_old_val/prop_dist_new_val)
  acc.prob <- prior.ratio*lik.ratio*trans.ratio 
  
  return(list(prior.ratio = prior.ratio, lik.ratio = lik.ratio, 
              trans.ratio = trans.ratio, alpha = min(1, acc.prob)))
}

cart_log_lik_copula <- function(tree_top, U1, U2, X, cop_type){ # check input
  tree.at.obs <- get_value_tree(tree_top, X)
  if(cop_type == "Gauss"){
    log.prob.obs <- gaussian_copula_loglik(tree.at.obs, U1, U2)
  }else if(cop_type == "Frank"){
    log.prob.obs <- loglik_frank(tree.at.obs, U1, U2)
  }else if(cop_type == "Clayton"){
    log.prob.obs <- loglik_clayton(tree.at.obs, U1, U2)
  }else if(cop_type == "gumbel"){
    log.prob.obs <- loglik_gumbel(tree.at.obs, U1, U2)
  }else if(cop_type == "t"){
    log.prob.obs <- loglik_t(tree.at.obs, U1, U2)
  }
  
  return(sum(log.prob.obs))
}

set_term_node_value_copula <- function(node.idx, tree_top, mu, sigma, cop_type){
  if(is.null(tree_top$left) & is.null(tree_top$right)){
    if(tree_top$node.idx == node.idx){
      if(cop_type == "Gauss" || cop_type == "t"){
        tree_top$value = rnorm(1, mean = mu, sd = sigma)
      } else if(cop_type == "gumbel") {
        tree_top$value = rlnorm(1, meanlog = mu, sdlog = sigma) + 1
      } else if(cop_type == "Clayton") {
        tree_top$value = rlnorm(1, meanlog = mu, sdlog = sigma)
      }
      return(tree_top)
    } else {
      return(tree_top)  
    }
  }  else {
    tree.left <- set_term_node_value_copula(node.idx, tree_top$left, mu, sigma, cop_type)
    tree.right <- set_term_node_value_copula(node.idx, tree_top$right, mu, sigma, cop_type)
    return(list(left = tree.left, right = tree.right, node.idx = tree_top$node.idx,
                cond = tree_top$cond))
  }
}

gaussian_copula_loglik <- function(rho, u, v) {
  
  if(any(abs(rho)>=1)) return(-Inf)
  
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

loglik_t <- function(rho, u, v, df = 3) {
  
  if(any(abs(rho)>=1)) return(-Inf)
  
  # Inverse t CDF (quantiles)
  x <- qt(u, df)
  y <- qt(v, df)
  
  # Univariate t densities
  fx <- dt(x, df)
  fy <- dt(y, df)
  
  # Bivariate t density (manual formula)
  denom <- sqrt(1 - rho^2)
  z <- (x^2 - 2 * rho * x * y + y^2) / (df * (1 - rho^2))
  
  A <- gamma((df + 2) / 2) / (gamma(df / 2) * df * pi * denom)
  B <- (1 + z)^(-(df + 2) / 2)
  f_biv <- A * B
  
  # Copula density
  copula_density <- f_biv / (fx * fy)
  
  # Log-likelihood
  log_lik <- sum(log(copula_density))
  return(log_lik)
}

loglik_clayton <- function(theta, u, v) {
  if(min(theta)<0) return(-Inf)
  
  log_density <- log(1+theta) + (-1-theta)*log(u*v) + (-1/theta-2) * log((u^(-theta)+v^(-theta)-1))
  return(sum(log_density))
}


loglik_gumbel <- function(theta, u, v) {
  
  if(any(theta)<1) return(-Inf)
  
  # Transformations
  A <- ((-log(u))^theta+(-log(v))^theta)^(1/theta)
  
  # Density
  log_density <- -A -(log(u) + log(v)) + (-2+2/theta)*log((-log(u))^theta+ (-log(v))^theta) + (theta-1) *log(log(u)*log(v)) + log(1+(theta-1)*((-log(u))^theta+(-log(v))^theta)^(-1/theta))
  
  return(sum(log_density))
}

log_prior_tbeta <- function(rho, bound = 1, alpha_val, beta_val) {
  if (any(abs(rho) >= 1)) {
    return(rep(-Inf, length(rho)))
  }
  
  return((alpha_val - 1)*log(bound + rho) + (beta_val - 1)*log(bound - rho))
}

logprior_inv_gamma <- function(p, alpha_val, beta_val) {
  # Ensure p is within (-1, 1)
  if (any(abs(p) >= 1)) {
    return(rep(-Inf, length(p)))
  }
  
  y <- -log(1 - p^2)
  
  dens <- (alpha_val * log(beta_val) + log(abs(p)) + (-alpha_val - 1) * log(y) - 
             beta_val / y - log(1 - p^2) - log(gamma(alpha_val))) 
  
  return(dens)
}

logprior_log_normal <- function(p, log_nor_mu, log_nor_sigma) {
  if (any(abs(p) >= 1)) {
    return(rep(-Inf, length(p)))
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
      return(gaussian_copula_loglik(rho, u, v) + log_prior_tbeta(rho, alpha_val = alpha_val, beta_val = beta_val))
    }else{
      if(prior_type == "LN"){
        return(gaussian_copula_loglik(rho, u, v) + logprior_log_normal(rho, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma))
      }else{
        return(gaussian_copula_loglik(rho, u, v) + logprior_inv_gamma(rho, alpha_val = alpha_val, beta_val = beta_val))
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
  }else if(cop_type == "gumbel"){
    if(prior_type == "IG"){
      return(loglik_gumbel(rho, u, v) + log(dgamma(rho, shape = alpha_val, rate = beta_val)))
    }else{
      if(prior_type == "LN"){
        return(loglik_gumbel(rho, u, v) + log(dlnorm(rho, meanlog = log_nor_mu, sdlog = log_nor_sigma)))
      }
    }
  }else if(cop_type == "t"){
    if(prior_type == "B"){
      return(loglik_t(rho, u, v) + log_prior_tbeta(rho, alpha_val = alpha_val, beta_val = beta_val))
    }else{
      if(prior_type == "LN"){
        return(loglik_t(rho, u, v) + logprior_log_normal(rho, log_nor_mu = log_nor_mu, log_nor_sigma = log_nor_sigma))
      }else{
        return(loglik_t(rho, u, v) + logprior_inv_gamma(rho, alpha_val = alpha_val, beta_val = beta_val))
      }
    }
  }
}

param_gauss = function(tau) return(sin(tau*pi/2))
param_gumbel = function(tau) return(1/(1-tau))
param_clayton = function(tau) return((2*tau)/(1-tau))

conv_diag = function(post_data, panel_name, n.burn, n.thin){
  
  post_mat = (post_data %>% filter(panel.name == panel_name))$y
  
  post_mat_red = na.omit((post_mat[n.burn+1:length(post_mat)])[c(rep(NA,(n.thin-1)), TRUE)])
  
  post_data_mcmc = mcmc(post_mat_red, start = 1)
  
  auto_corr_sum = autocorr.diag(post_data_mcmc)
  
  eff_size = effectiveSize(post_data_mcmc)
  
  geweke_score = geweke.diag(post_data_mcmc)
  
  
  summ_dat = cbind(auto_corr_sum[2], 100*(eff_size/length(post_data_mcmc)), geweke_score$z)
  
  colnames(summ_dat) = c("auto-correlation", "ess (%)", "geweke-score")
  rownames(summ_dat) = panel_name
  
  return(summ_dat)
  
}

slice_sampler <- function(log_post, x0, w = 1, N = 1) {
  samples <- numeric(N)
  x <- x0
  
  for (i in 1:N) {
    # Evaluate log density at current position
    log_y <- log_post(x) - rexp(1)
    
    # Stepping out procedure
    L <- x - runif(1, 0, w)
    R <- L + w
    while (log_post(L) > log_y) L <- L - w
    while (log_post(R) > log_y) R <- R + w
    
    # Sample new x from the slice
    repeat {
      x_new <- runif(1, L, R)
      if (log_post(x_new) >= log_y) break
      # Shrink the interval
      if (x_new < x) {
        L <- x_new
      } else {
        R <- x_new
      }
    }
    
    x <- x_new
    samples[i] <- x
  }
  
  return(samples)
}