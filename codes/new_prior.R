f.odd <- function(n) {
  as.integer(n %% 2 == 1)
}

f.even <- function(n) {
  as.integer(n %% 2 == 0)
}

## if the number of terminal nodes is even
n.cost.even <- function(gam, n){
  num_ = exp(gam*n) - 1
  den_ = exp(gam*n)*(1 - exp(-2*gam))
  num_/den_
}

n.cost.odd <- function(gam, n){
  num_ <- exp(-gam)*(1 - exp(-gam*(n-1)))
  den_ <- (1 - exp(-2*gam))
  num_/den_
}

prior.nterm <- function(n, omeg){
  exp(-omeg*n)*(exp(omeg) - 1)
}

count.tree <- function(n, delt) {
  
  nl <- (n + delt) / 2
  nr <- (n - delt) / 2
  
  # Use log-factorials to avoid overflow
  log_cnl <- lfactorial(2 * (nl - 1)) - lfactorial(nl) - lfactorial(nl - 1)
  log_cnr <- lfactorial(2 * (nr - 1)) - lfactorial(nr) - lfactorial(nr - 1)
  
  count <- exp(log_cnl + log_cnr)
  
  if (delt == 0) {
    return(count)
  } else {
    return(2 * count)
  }
}

prior.delta <- function(delt, gam, n) {
  
  if (n == 1) {
    return(ifelse(delt == 1, 1, 0))
  }
  
  if (gam > 0) {
    n.cost <- if (f.odd(n)) {
      n.cost.odd(gam, n)
    } else {
      n.cost.even(gam, n)
    }
    return(exp(-gam * delt) / n.cost)
  } else {
    poss.delt <- if (f.odd(n)) {
      seq(1, n - 2, by = 2)
    } else {
      seq(0, n - 2, by = 2)
    }
    
    prior_val <- numeric(length(delt))
    is_valid <- delt %in% poss.delt
    prior_val[is_valid] <- 1 / length(poss.delt)
    return(prior_val)
  }
}


joint.prior.new <- function(n, delt, omeg, gam) {
  
  if (n == 1) {
    return(prior.nterm(n, omeg))
  }
  
  n.tree <- count.tree(n, delt)
  delta.prior <- prior.delta(delt, gam, n)
  nterm.prior <- prior.nterm(n, omeg)
  
  nterm.prior * delta.prior / n.tree
}


joint.prior.new.tree <- function(tree_top, omeg, gam){
  
  nterm_ <- get_num_terminal_nodes(tree_top)
  delta_ <- abs(get_num_terminal_nodes(tree_top$left) - get_num_terminal_nodes(tree_top$right))
  joint.prior.new(n = nterm_, delt = delta_, omeg = omeg, gam = gam)
}

# CDF of number of terminal nodes
F_nterm <- function(n, omeg){
  1 - exp(-omeg*n)
}

# inverse of CDF of number of terminal nodes
F_nterm_inv <- function(p, omeg){
  -log(1 - p)/omeg
}

# sample the number of terminal nodes prior distribution
rnterm <- function(size, omeg){
  uu <- runif(size)
  ceiling(F_nterm_inv(uu, omeg))
}



# sample from the conditional distribution of delta
rdelta.cond <- function(size, gamm, nl) {
  
  if (nl == 1) {
    return(rep(0, size))
  }
  
  delta.v <- if (f.odd(nl)) {
    seq(1, nl - 2, by = 2)
  } else {
    seq(0, nl - 2, by = 2)
  }
  
  if (length(delta.v) > 1) {
    prob.v <- prior.delta(delta.v, gamm, nl)
    sample(delta.v, size = size, replace = TRUE, prob = prob.v)
  } else {
    rep(delta.v, size)
  }
}

# sample from the marginal of rdelta
rdelta <- function(size, omeg, gamm) {
  
  nl.v <- rnterm(size, omeg)
  output <- mapply(function(nl) rdelta.cond(1, gamm, nl), nl.v)
  return(as.integer(output))
}


rdepth <- function(size, omeg, gamm){
  
  nl.v <- rnterm(size, omeg)
  delt.v <- vapply(nl.v, \(x) rdelta.cond(1, gamm, x), 0)
  depth.v <- vapply(1:size, \(x) get_depth(generate_random_binary_tree_n_delta(n = nl.v[x],
                                                                               delt = delt.v[x])), 0)
  return(depth.v)
}


